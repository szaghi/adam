!< ADAM, PRISM charge-conserving PIC current (londrillo, totorica, esirkepov_Q<N>), FNL backend.

#include "fundal.H"

module adam_prism_fnl_pic_conserving_object
!< ADAM, PRISM charge-conserving PIC current (londrillo, totorica, esirkepov_Q<N>), FNL backend.
!<
!< Device port of the CPU charge-conserving current chain of `prism_cpu_object` (integrate_rk_ssp_pic_charge_conserving and
!< the helpers it calls). Host-side twins of the CPU matrix builders/factorization live here (private) so the CPU backend is
!< left untouched; the arithmetic of every kernel follows the CPU routine named in its header line.
!<
!< Multi-block and multi-rank (single refinement level): every solver works on a *virtual block*, the whole domain seen as
!< one `NX x NY x NZ` block (NX = nbx*ni, ...). The conserving operators are global along grid lines (the inverse of the
!< londrillo line matrix is dense, the matched face flux is a prefix sum from the domain boundary, the P/Q
!< filters reach up to 7/5 cells), so a single-block solve of the whole domain is the reference result: on the virtual block the
!< phase-1 single-block kernels run unchanged and the multi-block result equals the single-block one by construction.
!< The real blocks are touched only at the two interfaces: (i) inputs, Esirkepov directly in global cell indices, the
!< matched rho through the regular block deposit + `reduce_ghost_local_gpu` + gather of the block interiors; (ii) outputs, J (and the filtered rho) scattered from the virtual block to every real cell, ghosts
!< included, so interface ghosts get the neighbour values and physical ghosts the single-block odd/zero values.
!< With one block the virtual block is the block itself (same origin, spacing and cell centers): phase-1 arithmetic.
!< Several ranks: the virtual block is replicated on every rank. Each rank adds the sources of its own particles on the
!< virtual cells (block interiors and in-domain ghosts, so that no ghost fold is needed), the complete sources are summed
!< over the ranks (`reduce_virtual`, MPI_Allreduce staged through the host), every rank runs the same virtual solve and
!< scatters the result on its own blocks. The per-particle cleanup sums the particle currents of each rank, then reduces
!< the virtual current.
!<
!< Device layouts (FNL convention, block index first):
!<+ real q-like fields: `[nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, nv]`;
!<+ virtual fields: the same with nb=1 and NX, NY, NZ (H, scratches, line buffers, the virtual current `jv_gpu`);
!<+ directional charge changes H: `[1, NX, NY, NZ, 3]` (interior only, the CPU solvers never read H ghosts);
!<+ particles: `[particle, variable]` (coalesced across particles).
!<
!< Like the CPU reference, a particle whose support touches the domain boundary stops the run (zero-flux closure, no
!< particle boundary conditions).

! ADAM classes, libraries, parameters
use :: adam_common_library
! ADAM FNL classes, libraries, parameters
use :: adam_fnl_library
! PRISM common classes, libraries, parameters
use :: adam_prism_common_library
! PRISM FNL classes
use :: adam_prism_fnl_pic_object
! third party modules
use :: fundal
use :: mpi
use :: penf

implicit none
private
public :: prism_fnl_pic_conserving_object
public :: CONSERVING_SOLVER_LONDRILLO
public :: CONSERVING_SOLVER_MODIFIED

integer(I4P), parameter :: CONSERVING_SOLVER_LONDRILLO = 2_I4P !< Esirkepov face fluxes + inverse reconstruction.
integer(I4P), parameter :: CONSERVING_SOLVER_MODIFIED  = 3_I4P !< Matched current, totorica or esirkepov_Q<N> (P/Q kernels).

integer(I4P), parameter :: SHAPE_NGP      = 0_I4P !< NGP particle shape.
integer(I4P), parameter :: SHAPE_CIC      = 1_I4P !< CIC particle shape.
integer(I4P), parameter :: SHAPE_TSC      = 2_I4P !< TSC particle shape.
integer(I4P), parameter :: SHAPE_CUBIC    = 3_I4P !< Cubic particle shape.
integer(I4P), parameter :: SHAPE_QUARTIC  = 4_I4P !< Quartic particle shape.
integer(I4P), parameter :: SHAPE_QUINTIC  = 5_I4P !< Quintic particle shape.
integer(I4P), parameter :: SHAPE_SEXTIC   = 6_I4P !< Sextic particle shape.
integer(I4P), parameter :: SHAPE_GAUSSIAN = 7_I4P !< Gaussian particle shape.

integer(I4P), parameter :: MAX_SPAN = 32_I4P !< Max 1D cells touched by one particle trajectory (gang-private weights).
integer(I4P), parameter :: LQ_SPAN = MAX_SPAN + 10_I4P !< Box cells widened by the Q taps (box current, nq <= 10).
integer(I4P), parameter :: LP_SPAN = MAX_SPAN + 14_I4P !< Box cells widened by the P taps (box current, radius <= 7).

integer(I4P), parameter :: ERR_OUTSIDE  = 1_I4P !< Particle outside the domain.
integer(I4P), parameter :: ERR_CHARGE   = 3_I4P !< Particle charge changed during the trajectory.
integer(I4P), parameter :: ERR_BOUNDARY = 4_I4P !< Particle shape touches the domain boundary.
integer(I4P), parameter :: ERR_FILTER   = 5_I4P !< esirkepov_Q<N> P filter support reaches the domain boundary.
integer(I4P), parameter :: ERR_SPAN     = 6_I4P !< Trajectory support wider than MAX_SPAN (FNL-only limit).
integer(I4P), parameter :: ERR_GAUSS    = 7_I4P !< Empty Gaussian deposition support.
integer(I4P), parameter :: ERR_PU       = 8_I4P !< B-spline weights do not sum to one on the deposition stencil.

real(R8P), parameter :: TAU_TAIL = 1.e-11_R8P !< Relative tail threshold of the particle current cleanup (CPU twin).
real(R8P), parameter :: TAU_PU   = 1.e-10_R8P !< Max |raw B-spline weight sum - 1| before the closure (CPU twin).

type :: prism_fnl_pic_conserving_object
   !< Charge-conserving PIC current, device state and kernels.
   logical      :: enabled      = .false. !< Conserving current active.
   logical      :: prepared     = .false. !< Geometry-dependent caches built.
   integer(I4P) :: solver       = 0_I4P   !< Selected solver (CONSERVING_SOLVER_*).
   integer(I4P) :: shape        = 0_I4P   !< Particle shape (SHAPE_*).
   logical      :: filter_deposition = .false. !< Binomial filter on deposition.
   logical      :: tail_cleanup = .false. !< Matched-current per-particle tail cleanup.
   logical      :: is_periodic(3) = .false. !< Periodic particle closure along each axis (host pic%bc).
   logical      :: box_current = .false.  !< Matched current and charge built per particle on its box (periodic faces).
   real(R8P)    :: sigma        = 0._R8P  !< Gaussian width.
   real(R8P)    :: cutoff_sigma = 0._R8P  !< Gaussian cutoff in sigma units.
   integer(I4P) :: gaussian_support_cells = -1_I4P !< Gaussian support radius in cells.
   integer(I4P) :: fdv_order    = 0_I4P   !< FD order.
   integer(I4P) :: hs           = 0_I4P   !< FD half stencil.
   integer(I4P) :: ni = 0_I4P, nj = 0_I4P, nk = 0_I4P, ngc = 0_I4P !< Virtual block sizes (whole domain), ghost cells.
   integer(I4P) :: nb = 1_I4P            !< Virtual blocks (always 1).
   integer(I4P) :: blocks_number = 1_I4P !< Virtual blocks (always 1).
   integer(I4P) :: bni = 0_I4P, bnj = 0_I4P, bnk = 0_I4P !< Real block sizes.
   integer(I4P) :: bnb = 0_I4P           !< Real allocated blocks.
   integer(I4P) :: bblocks = 0_I4P       !< Real actual blocks (cache key).
   real(R8P)    :: vemin(3) = 0._R8P     !< Virtual block origin (domain minimum).
   real(R8P)    :: vdx(3) = 0._R8P       !< Cell size (single refinement level).
   integer(I4P) :: nv = 0_I4P            !< Field variables number.
   integer(I4P) :: particle_number = 0_I4P !< Particles number owned by this rank.
   integer(I4P) :: particle_number_global = 0_I4P !< Particles number, all ranks.
   integer(I4P) :: nrk = 0_I4P           !< RK stages number.
   integer(I4P) :: esirkepov_q = 0_I4P   !< Q filter coefficients of esirkepov_Q<N>.
   character(len=:), allocatable :: solver_name !< Input name of the solver (messages).
   integer(I4P) :: radius = 0_I4P, qfirst = 0_I4P, nq = 0_I4P !< Modified P/Q kernels geometry.
   real(R8P)    :: shift(3) = 0._R8P     !< Base-shape position shift: dx/2 for the odd Q (dual nodes), 0 otherwise.
   real(R8P)    :: p(-7:7) = 0._R8P      !< Modified P kernel (cell centers / transverse lines).
   real(R8P)    :: qw(0:9) = 0._R8P      !< Modified Q kernel (normal faces to centers).
   real(R8P)    :: residual_max = 0._R8P !< Max current-solver residual of the current step.
   logical      :: coherent = .false.    !< Coherent matched gather (field_weighting_model totorica or esirkepov_Q<N>).
   integer(I4P) :: gradius = 0_I4P       !< Coherent gather: taps of its P kernel.
   real(R8P)    :: gshift(3) = 0._R8P    !< Coherent gather: dx/2 for the odd esirkepov_Q<N> (dual nodes), 0 otherwise.
   ! device data
   integer(I4P), pointer :: off_gpu(:,:)              => null() !< Real block origin in the virtual block, cells [bnb,3].
   real(R8P),    pointer :: vx_gpu(:)                 => null() !< Virtual x cell centers [NX].
   real(R8P),    pointer :: vy_gpu(:)                 => null() !< Virtual y cell centers [NY].
   real(R8P),    pointer :: vz_gpu(:)                 => null() !< Virtual z cell centers [NZ].
   real(R8P),    pointer :: jv_gpu(:,:,:,:,:)         => null() !< Virtual current [1,ghosts...,3].
   real(R8P),    pointer :: p_gpu(:)                  => null() !< P kernel [-7:7].
   real(R8P),    pointer :: qw_gpu(:)                 => null() !< Q kernel [0:9].
   real(R8P),    pointer :: fv1_gpu(:)                => null() !< FV1_CC(1:hs,hs), face reconstruction coefficients.
   real(R8P),    pointer :: h_gpu(:,:,:,:,:,:)        => null() !< Stage directional charge changes [nb,ni,nj,nk,3,nrk].
   real(R8P),    pointer :: q_next_gpu(:,:)           => null() !< Particle stage target [np,8].
   real(R8P),    pointer :: q_ref_gpu(:,:)            => null() !< Particle reference (step start / t0 minus) [np,8].
   real(R8P),    pointer :: q_hist_gpu(:,:,:)         => null() !< Particle stage targets history [np,8,nrk].
   real(R8P),    pointer :: q_mixed_gpu(:,:)          => null() !< Shifted positions, odd esirkepov_Q<N> rho [np,8].
   real(R8P),    pointer :: work_a_gpu(:,:,:,:)       => null() !< Interior scratch [nb,ni,nj,nk].
   real(R8P),    pointer :: work_b_gpu(:,:,:,:)       => null() !< Interior scratch [nb,ni,nj,nk].
   real(R8P),    pointer :: line_gpu(:)               => null() !< Line-major scratch for line solves/scans.
   real(R8P),    pointer :: gather_gpu(:,:,:,:,:)     => null() !< Real-block field copy, external-field gather [bnb,...,nv].
   real(R8P),    pointer :: gp_gpu(:)                 => null() !< Coherent gather P kernel [-7:7].
   real(R8P),    pointer :: gq_gpu(:,:,:,:,:)         => null() !< Coherent gather: E, B filtered by P^T [bnb,...,6].
   real(R8P),    pointer :: q_dual_gpu(:,:)           => null() !< Coherent gather: stage particles at x-dx/2 [np,8].
   real(R8P),    pointer :: jp_gpu(:,:,:,:,:)         => null() !< Single-particle virtual current, cleanup [1,ghosts...,3].
   real(R8P),    pointer :: src_gpu(:,:,:,:,:)        => null() !< Single-particle source, cleanup [nb,ni,nj,nk,3].
   real(R8P),    pointer :: src_stage_gpu(:,:,:,:,:,:)=> null() !< Single-particle stage sources [nb,ni,nj,nk,3,nrk].
   real(R8P),    pointer :: mat_x_gpu(:,:,:)          => null() !< Factored x-line matrices [ni,ni,nb].
   real(R8P),    pointer :: mat_y_gpu(:,:,:)          => null() !< Factored y-line matrices [nj,nj,nb].
   real(R8P),    pointer :: mat_z_gpu(:,:,:)          => null() !< Factored z-line matrices [nk,nk,nb].
   integer(I4P), pointer :: piv_x_gpu(:,:)            => null() !< x-line pivots [ni,nb].
   integer(I4P), pointer :: piv_y_gpu(:,:)            => null() !< y-line pivots [nj,nb].
   integer(I4P), pointer :: piv_z_gpu(:,:)            => null() !< z-line pivots [nk,nb].
   contains
      procedure, pass(self) :: destroy                     !< Free device data.
      procedure, pass(self) :: initialize                  !< Initialize from host configuration.
      procedure, pass(self) :: prepare                     !< Build geometry-dependent caches.
      procedure, pass(self) :: allocate_gather             !< Allocate the external-field gather copy.
      procedure, pass(self) :: initialize_coherent_gather  !< Kernels and buffers of the coherent matched gather.
      procedure, pass(self) :: coherent_field_weighting    !< Gather E, B with P^T and the base shape.
      procedure, pass(self) :: copy_field                  !< Copy a q-like device field.
      procedure, pass(self) :: copy_particles              !< Copy a device particle array.
      procedure, pass(self) :: build_stage_target          !< q_next = q_pic + dt sum_r c_r k_r.
      procedure, pass(self) :: build_centered_trajectory   !< t=0 virtual trajectory q_pic -/+ dt/2 v.
      procedure, pass(self) :: esirkepov_displacement      !< Esirkepov directional charge change.
      procedure, pass(self) :: combine_stage_source        !< H_s = (H_s - sum_r c_r H_r)/c_s.
      procedure, pass(self) :: solve_dispatch              !< Select the configured current solver.
      procedure, pass(self) :: solve_lines                 !< Londrillo line solves.
      procedure, pass(self) :: solve_modified              !< Matched current (totorica, esirkepov_Q<N>).
      procedure, pass(self) :: solve_modified_with_cleanup !< Per-particle modified current with tail cleanup.
      procedure, pass(self) :: deposit_modified_charge     !< Unfiltered deposit + P_c filter of rho.
      procedure, pass(self) :: deposit_matched_current_box !< Matched current per particle on its box, folded (periodic).
      procedure, pass(self) :: store_box_current           !< Real q(J) = virtual stage current of the box solver.
      procedure, pass(self) :: impose_current_ghosts       !< Odd/zero current ghosts on the physical faces of real blocks.
      procedure, pass(self) :: reduce_virtual4             !< Sum a rank-4 virtual-block array over the ranks.
      procedure, pass(self) :: reduce_virtual5             !< Sum a rank-5 virtual-block array over the ranks.
      generic               :: reduce_virtual => reduce_virtual4, reduce_virtual5 !< Sum a virtual-block array over the ranks.
      ! private methods
      procedure, pass(self), private :: build_virtual_geometry  !< Virtual block of the whole domain.
      procedure, pass(self), private :: scatter_current         !< Real q(J) = virtual J, ghosts included.
      procedure, pass(self), private :: impose_virtual_ghosts   !< Odd/zero current ghosts of the virtual block.
      procedure, pass(self), private :: impose_virtual_ghosts_periodic !< Virtual current ghosts with periodic axes.
      procedure, pass(self), private :: deposit_matched_charge_box !< Matched rho per particle, folded (periodic).
      procedure, pass(self), private :: check_source_support    !< Particle source inside its trajectory support box.
      procedure, pass(self), private :: transposed_filter_axis  !< P^T convolution of the coherent gather along one axis.
endtype prism_fnl_pic_conserving_object

contains
   ! public methods
   subroutine destroy(self)
   !< Free device data owned by this object.
   class(prism_fnl_pic_conserving_object), intent(inout) :: self !< Conserving current object.

   call free_r1(self%p_gpu)
   call free_r1(self%qw_gpu)
   call free_r1(self%gp_gpu)
   call free_r2(self%q_dual_gpu)
   call free_r5(self%gq_gpu)
   call free_r1(self%fv1_gpu)
   call free_r1(self%line_gpu)
   call free_r1(self%vx_gpu)
   call free_r1(self%vy_gpu)
   call free_r1(self%vz_gpu)
   call free_i2(self%off_gpu)
   call free_r2(self%q_next_gpu)
   call free_r2(self%q_ref_gpu)
   call free_r2(self%q_mixed_gpu)
   call free_r3(self%q_hist_gpu)
   call free_r3(self%mat_x_gpu)
   call free_r3(self%mat_y_gpu)
   call free_r3(self%mat_z_gpu)
   call free_r4(self%work_a_gpu)
   call free_r4(self%work_b_gpu)
   call free_r5(self%gather_gpu)
   call free_r5(self%jp_gpu)
   call free_r5(self%jv_gpu)
   call free_r5(self%src_gpu)
   call free_r6(self%h_gpu)
   call free_r6(self%src_stage_gpu)
   call free_i2(self%piv_x_gpu)
   call free_i2(self%piv_y_gpu)
   call free_i2(self%piv_z_gpu)
   self%enabled  = .false.
   self%prepared = .false.
   self%coherent = .false.
   contains
      subroutine free_r1(a)
      real(R8P), pointer, intent(inout) :: a(:)
      if (associated(a)) then ; call dev_free(a, mydev) ; nullify(a) ; endif
      endsubroutine free_r1
      subroutine free_r2(a)
      real(R8P), pointer, intent(inout) :: a(:,:)
      if (associated(a)) then ; call dev_free(a, mydev) ; nullify(a) ; endif
      endsubroutine free_r2
      subroutine free_r3(a)
      real(R8P), pointer, intent(inout) :: a(:,:,:)
      if (associated(a)) then ; call dev_free(a, mydev) ; nullify(a) ; endif
      endsubroutine free_r3
      subroutine free_r4(a)
      real(R8P), pointer, intent(inout) :: a(:,:,:,:)
      if (associated(a)) then ; call dev_free(a, mydev) ; nullify(a) ; endif
      endsubroutine free_r4
      subroutine free_r5(a)
      real(R8P), pointer, intent(inout) :: a(:,:,:,:,:)
      if (associated(a)) then ; call dev_free(a, mydev) ; nullify(a) ; endif
      endsubroutine free_r5
      subroutine free_r6(a)
      real(R8P), pointer, intent(inout) :: a(:,:,:,:,:,:)
      if (associated(a)) then ; call dev_free(a, mydev) ; nullify(a) ; endif
      endsubroutine free_r6
      subroutine free_i2(a)
      integer(I4P), pointer, intent(inout) :: a(:,:)
      if (associated(a)) then ; call dev_free(a, mydev) ; nullify(a) ; endif
      endsubroutine free_i2
   endsubroutine destroy

   subroutine initialize(self, pic, grid, field, nb, nv, nrk, fdv_order, hs)
   !< Initialize from the host configuration and allocate the step-invariant device buffers.
   !< Twin of the CPU validation in prism_cpu_object%initialize plus the allocations of integrate_rk_ssp_pic_charge_conserving.
   !< The virtual block (whole domain) is built here: the blocks must already be in place (uniform refinement included).
   class(prism_fnl_pic_conserving_object), intent(inout) :: self      !< Conserving current object.
   type(prism_pic_object),                 intent(in)    :: pic       !< Host PIC object.
   type(grid_object),                      intent(in)    :: grid      !< Grid.
   type(field_object),                     intent(in)    :: field     !< Host field (block geometry).
   integer(I4P),                           intent(in)    :: nb        !< Allocated real blocks.
   integer(I4P),                           intent(in)    :: nv        !< Field variables number.
   integer(I4P),                           intent(in)    :: nrk       !< RK stages number.
   integer(I4P),                           intent(in)    :: fdv_order !< FD order.
   integer(I4P),                           intent(in)    :: hs        !< FD half stencil.
   integer(I4P)                                          :: np, ni, nj, nk, ngc, nl, ierr

   call self%destroy
   self%enabled = .true.
   select case(trim(pic%current_conserving_solver))
   case(LONDRILLO_CURRENT_CONSERVING_SOLVER)                                     ; self%solver = CONSERVING_SOLVER_LONDRILLO
   case(TOTORICA_CURRENT_CONSERVING_SOLVER, ESIRKEPOV_MODIFIED_CURRENT_CONSERVING_SOLVER); self%solver = CONSERVING_SOLVER_MODIFIED
   case default
      call mpih%error_stop(msg=': unsupported current_conserving_solver '//trim(pic%current_conserving_solver))
   endselect
   select case(trim(pic%particle_weighting_model))
   case('NGP')      ; self%shape = SHAPE_NGP
   case('CIC')      ; self%shape = SHAPE_CIC
   case('TSC')      ; self%shape = SHAPE_TSC
   case('cubic')    ; self%shape = SHAPE_CUBIC
   case('quartic')  ; self%shape = SHAPE_QUARTIC
   case('quintic')  ; self%shape = SHAPE_QUINTIC
   case('sextic')   ; self%shape = SHAPE_SEXTIC
   case('Gaussian') ; self%shape = SHAPE_GAUSSIAN
   case default
      call mpih%error_stop(msg=': unsupported particle weighting in Esirkepov current')
   endselect
   self%esirkepov_q            = pic%esirkepov_q
   self%solver_name            = pic%conserving_solver_name()
   self%filter_deposition      = pic%filter_deposition
   self%tail_cleanup           = pic%esirkepov_tail_cleanup
   self%is_periodic            = pic%bc%is_periodic
   self%sigma                  = pic%sigma
   self%cutoff_sigma           = pic%cutoff_sigma
   self%gaussian_support_cells = pic%gaussian_support_cells
   self%fdv_order              = fdv_order
   self%hs                     = hs
   self%bni = grid%ni ; self%bnj = grid%nj ; self%bnk = grid%nk ; self%ngc = grid%ngc
   self%bnb = nb ; self%nv = nv ; self%nrk = nrk
   self%nb = 1_I4P ; self%blocks_number = 1_I4P
   call self%build_virtual_geometry(field=field, grid=grid)
   self%particle_number = pic%particle_number
   self%particle_number_global = pic%particle_number_global
   np = self%particle_number ; ni = self%ni ; nj = self%nj ; nk = self%nk ; ngc = self%ngc

   if (self%tail_cleanup .and. self%solver /= CONSERVING_SOLVER_MODIFIED) &
      call mpih%error_stop(msg=': esirkepov_tail_cleanup requires totorica or an esirkepov_Q<N> solver')
   ! periodic particle closure: the matched currents built on the trajectory boxes only (CPU twin:
   ! check_conserving_single_block); without periodic faces the virtual line solves are kept bitwise
   if (any(self%is_periodic)) then
      if (self%solver /= CONSERVING_SOLVER_MODIFIED) &
         call mpih%error_stop(msg=': '//self%solver_name//' does not support periodic faces yet (totorica and '// &
                                  'esirkepov_Q<N> do)')
      if (self%tail_cleanup) &
         call mpih%error_stop(msg=': esirkepov_tail_cleanup does not support periodic faces (and the matched current '// &
                                  'built on the trajectory boxes has no tails: disable it)')
      self%box_current = .true.
   endif
   if (self%shape == SHAPE_GAUSSIAN) then
      if (self%sigma <= 0._R8P .or. self%cutoff_sigma <= 0._R8P .or. self%gaussian_support_cells < 0_I4P) &
         call mpih%error_stop(msg=': invalid Gaussian width, cutoff or support in Esirkepov current')
   endif
   if (ngc > min(ni,nj,nk)) &
      call mpih%error_stop(msg=': esirkepov needs at least ngc interior cells per direction')
   self%shift = 0._R8P
   if (self%solver == CONSERVING_SOLVER_MODIFIED) then
      call pic%check_esirkepov_q_order(fdv_order=fdv_order)
      if (trim(pic%current_conserving_solver) == TOTORICA_CURRENT_CONSERVING_SOLVER) then
         call totorica_kernels(hs=hs, p=self%p, qw=self%qw, radius=self%radius, qfirst=self%qfirst, nq=self%nq)
      else
         call modified_kernels(esirkepov_q=self%esirkepov_q, p=self%p, qw=self%qw, radius=self%radius, &
                               qfirst=self%qfirst, nq=self%nq)
      endif
      ! odd Q: the base shape lives on the dual nodes x_c+dx/2, i.e. the regular shape of the particle at x-dx/2
      if (mod(self%esirkepov_q, 2_I4P) == 1_I4P) self%shift = 0.5_R8P*self%vdx
      ! uploaded 1-based: p_gpu(s+8) = p(s), qw_gpu(t+1) = qw(t)
      call dev_assign_to_device(src=[self%p],  dst=self%p_gpu)
      call dev_assign_to_device(src=[self%qw], dst=self%qw_gpu)
   endif
   ! Stencil coefficients uploaded once: the kernels do not read the module tables FV1_CC/FD1_CC, whose device copy
   ! (`!$acc declare copyin` of a parameter) is not reliably initialized when no other kernel references it.
   call dev_assign_to_device(src=FV1_CC(1:hs,hs), dst=self%fv1_gpu)
   ! FNL-only limit: the Esirkepov weights of one trajectory live in gang-private arrays of MAX_SPAN cells
   if (2_I4P*(support_radius(self)+1_I4P)+2_I4P > MAX_SPAN) &
      call mpih%error_stop(msg=': particle support too wide for the FNL Esirkepov current (MAX_SPAN='// &
                               trim(str(MAX_SPAN,.true.))//')')

   call alloc_r6(self%h_gpu, [1,ni,nj,nk,3,nrk], [1,1,1,1,1,1], 'h_gpu')
   call alloc_r4(self%work_a_gpu, [1,ni,nj,nk], 'work_a_gpu')
   call alloc_r4(self%work_b_gpu, [1,ni,nj,nk], 'work_b_gpu')
   nl = (ni+1)*(nj+1)*(nk+1)
   call dev_alloc(fptr_dev=self%line_gpu, ubounds=[nl], lbounds=[1], init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih%error_stop(msg=': failed to allocate line_gpu in prism_fnl_pic_conserving_object')
   if (np > 0_I4P) then
      call alloc_r2(self%q_next_gpu,  [np,8], 'q_next_gpu')
      call alloc_r2(self%q_ref_gpu,   [np,8], 'q_ref_gpu')
      call alloc_r3(self%q_hist_gpu,  [np,8,nrk], 'q_hist_gpu')
      ! odd esirkepov_Q: shifted positions of the rho deposit
      if (any(self%shift /= 0._R8P)) call alloc_r2(self%q_mixed_gpu, [np,8], 'q_mixed_gpu')
   endif
   call alloc_r5(self%jv_gpu, [1,ni+ngc,nj+ngc,nk+ngc,3], [1,1-ngc,1-ngc,1-ngc,1], 'jv_gpu')
   if (self%tail_cleanup) then
      call alloc_r5(self%jp_gpu, [1,ni+ngc,nj+ngc,nk+ngc,3], [1,1-ngc,1-ngc,1-ngc,1], 'jp_gpu')
      if (np > 1_I4P) then
         call alloc_r5(self%src_gpu, [1,ni,nj,nk,3], [1,1,1,1,1], 'src_gpu')
         call alloc_r6(self%src_stage_gpu, [1,ni,nj,nk,3,nrk], [1,1,1,1,1,1], 'src_stage_gpu')
      endif
   endif
   contains
      subroutine alloc_r2(a, ub, label)
      real(R8P), pointer, intent(inout) :: a(:,:)
      integer(I4P),       intent(in)    :: ub(2)
      character(*),       intent(in)    :: label
      call dev_alloc(fptr_dev=a, ubounds=ub, lbounds=[1,1], init_value=0._R8P, ierr=ierr)
      if (ierr /= 0_I4P) call mpih%error_stop(msg=': failed to allocate '//label//' in prism_fnl_pic_conserving_object')
      endsubroutine alloc_r2
      subroutine alloc_r3(a, ub, label)
      real(R8P), pointer, intent(inout) :: a(:,:,:)
      integer(I4P),       intent(in)    :: ub(3)
      character(*),       intent(in)    :: label
      call dev_alloc(fptr_dev=a, ubounds=ub, lbounds=[1,1,1], init_value=0._R8P, ierr=ierr)
      if (ierr /= 0_I4P) call mpih%error_stop(msg=': failed to allocate '//label//' in prism_fnl_pic_conserving_object')
      endsubroutine alloc_r3
      subroutine alloc_r4(a, ub, label)
      real(R8P), pointer, intent(inout) :: a(:,:,:,:)
      integer(I4P),       intent(in)    :: ub(4)
      character(*),       intent(in)    :: label
      call dev_alloc(fptr_dev=a, ubounds=ub, lbounds=[1,1,1,1], init_value=0._R8P, ierr=ierr)
      if (ierr /= 0_I4P) call mpih%error_stop(msg=': failed to allocate '//label//' in prism_fnl_pic_conserving_object')
      endsubroutine alloc_r4
      subroutine alloc_r5(a, ub, lb, label)
      real(R8P), pointer, intent(inout) :: a(:,:,:,:,:)
      integer(I4P),       intent(in)    :: ub(5), lb(5)
      character(*),       intent(in)    :: label
      call dev_alloc(fptr_dev=a, ubounds=ub, lbounds=lb, init_value=0._R8P, ierr=ierr)
      if (ierr /= 0_I4P) call mpih%error_stop(msg=': failed to allocate '//label//' in prism_fnl_pic_conserving_object')
      endsubroutine alloc_r5
      subroutine alloc_r6(a, ub, lb, label)
      real(R8P), pointer, intent(inout) :: a(:,:,:,:,:,:)
      integer(I4P),       intent(in)    :: ub(6), lb(6)
      character(*),       intent(in)    :: label
      call dev_alloc(fptr_dev=a, ubounds=ub, lbounds=lb, init_value=0._R8P, ierr=ierr)
      if (ierr /= 0_I4P) call mpih%error_stop(msg=': failed to allocate '//label//' in prism_fnl_pic_conserving_object')
      endsubroutine alloc_r6
   endsubroutine initialize

   subroutine prepare(self, field, grid, force)
   !< Build the geometry-dependent caches: virtual block geometry and the factored line matrices.
   !< Twin of CPU ensure_londrillo_current_solver_cache (host factorization, one upload).
   !< The matrices span whole virtual lines (the domain extent), one per direction (single refinement level).
   class(prism_fnl_pic_conserving_object), intent(inout)        :: self  !< Conserving current object.
   type(field_object),                     intent(in)           :: field !< Host field.
   type(grid_object),                      intent(in)           :: grid  !< Grid.
   logical,                                intent(in), optional :: force !< Rebuild even if already prepared.
   real(R8P),    allocatable                                    :: mat(:,:,:), face_map(:,:)
   integer(I4P), allocatable                                    :: piv(:,:)
   integer(I4P)                                                 :: dir, n, b, face, col, nvirt(3)
   logical                                                      :: force_

   force_ = .false. ; if (present(force)) force_ = force
   if (self%prepared .and. (.not. force_) .and. self%bblocks == field%blocks_number) return
   nvirt = [self%ni, self%nj, self%nk]
   call self%build_virtual_geometry(field=field, grid=grid)
   if (any([self%ni, self%nj, self%nk] /= nvirt)) &
      call mpih%error_stop(msg=': the FNL conserving PIC virtual block changed size (refinement level change is unsupported)')
   if (self%solver /= CONSERVING_SOLVER_MODIFIED) then
      do dir=1, 3
         select case(dir)
         case(1) ; n = self%ni
         case(2) ; n = self%nj
         case(3) ; n = self%nk
         endselect
         allocate(mat(1:n,1:n,1:self%nb), piv(1:n,1:self%nb))
         mat = 0._R8P ; piv = 0_I4P
         do b=1, self%blocks_number
            allocate(face_map(0:n,1:n))
            call build_centered_reconstruction_matrix(hs=self%hs, n=n, a=face_map)
            do face=1, n-1
               mat(face,1:n,b) = face_map(face,1:n)
            enddo
            do col=1, n
               mat(n,col,b) = (-1._R8P)**col
            enddo
            deallocate(face_map)
            call factorize_matrix_pivot(a=mat(:,:,b), pivot=piv(:,b), n=n)
         enddo
         select case(dir)
         case(1)
            call dev_assign_to_device(src=mat, dst=self%mat_x_gpu) ; call dev_assign_to_device(src=piv, dst=self%piv_x_gpu)
         case(2)
            call dev_assign_to_device(src=mat, dst=self%mat_y_gpu) ; call dev_assign_to_device(src=piv, dst=self%piv_y_gpu)
         case(3)
            call dev_assign_to_device(src=mat, dst=self%mat_z_gpu) ; call dev_assign_to_device(src=piv, dst=self%piv_z_gpu)
         endselect
         deallocate(mat, piv)
      enddo
   endif
   self%prepared = .true.
   endsubroutine prepare

   subroutine allocate_gather(self)
   !< Allocate the field copy used to gather with external fields (CPU twin: q_gather).
   class(prism_fnl_pic_conserving_object), intent(inout) :: self !< Conserving current object.
   integer(I4P)                                          :: ierr !< Error status.

   if (associated(self%gather_gpu)) return
   call dev_alloc(fptr_dev=self%gather_gpu, ubounds=[self%bnb,self%bni+self%ngc,self%bnj+self%ngc,self%bnk+self%ngc,self%nv], &
                  lbounds=[1,1-self%ngc,1-self%ngc,1-self%ngc,1], init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih%error_stop(msg=': failed to allocate gather_gpu in prism_fnl_pic_conserving_object')
   endsubroutine allocate_gather

   subroutine initialize_coherent_gather(self, pic)
   !< Kernels and buffers of the coherent gather of totorica and esirkepov_Q<N> (pic%set_coherent_gather already ran). CPU
   !< twin: the gather_* state of prism_pic_object.
   class(prism_fnl_pic_conserving_object), intent(inout) :: self !< Conserving current object.
   type(prism_pic_object),                 intent(in)    :: pic  !< Host PIC object.
   integer(I4P)                                          :: ierr !< Error status.

   if (.not.pic%matched_gather()) return
   self%coherent = .true.
   self%gradius  = pic%gather_radius
   self%gshift   = 0._R8P
   if (pic%gather_dual) self%gshift = 0.5_R8P*self%vdx
   ! uploaded 1-based: gp_gpu(s+8) = gather_p(s)
   call dev_assign_to_device(src=[pic%gather_p], dst=self%gp_gpu)
   call dev_alloc(fptr_dev=self%gq_gpu, ubounds=[self%bnb,self%bni+self%ngc,self%bnj+self%ngc,self%bnk+self%ngc,6], &
                  lbounds=[1,1-self%ngc,1-self%ngc,1-self%ngc,1], init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih%error_stop(msg=': failed to allocate gq_gpu in prism_fnl_pic_conserving_object')
   if (self%particle_number > 0_I4P .and. any(self%gshift /= 0._R8P)) then
      call dev_alloc(fptr_dev=self%q_dual_gpu, ubounds=[self%particle_number,8], lbounds=[1,1], init_value=0._R8P, ierr=ierr)
      if (ierr /= 0_I4P) call mpih%error_stop(msg=': failed to allocate q_dual_gpu in prism_fnl_pic_conserving_object')
   endif
   endsubroutine initialize_coherent_gather

   subroutine deposit_matched_charge_box(self, pic_fnl, q, nv)
   !< Matched charge rho = P S rho^P built particle by particle and folded into the virtual block along the periodic axes,
   !< summed over the ranks, then scattered to every real cell (ghosts beyond a periodic face take the periodic image,
   !< physical ghosts 0). CPU twin: deposit_matched_charge_box. Used with periodic faces only.
   !<
   !< Parallel mapping: one gang per particle, the three 1D factors P w computed sequentially (one column per axis), the
   !< deposit one vector lane per cell, atomic.
   class(prism_fnl_pic_conserving_object), intent(inout) :: self    !< Conserving current object.
   type(prism_fnl_pic_object),             intent(inout) :: pic_fnl !< Device PIC helper.
   real(R8P),                              intent(inout) :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Real field.
   integer(I4P),                           intent(in)    :: nv      !< Charge density slot.
   real(R8P),    pointer                                 :: pk(:), wa(:,:,:,:), q_pic(:,:)
   integer(I4P), pointer                                 :: off(:,:)
   real(R8P)                                             :: w(MAX_SPAN,3), xc(MAX_SPAN), pw(LP_SPAN,3)
   real(R8P)                                             :: dxyz(3), x(3), emin(3), cutoff_limit, charge_density, sigma, acc
   real(R8P)                                             :: value, e1, e2, e3, d1, d2, d3, s1, s2, s3
   integer(I4P)                                          :: part, d, m, l, sm, li, lj, lk, perr, ierr, rad, radius, shp
   integer(I4P)                                          :: ni, nj, nk, nx, ny, nz, nbl, ngc, np, b, i, j, k, gi, gj, gk
   integer(I4P)                                          :: ii, jj, kk, ic(3), lo(3), rr(3), nn(3)
   logical                                               :: px, py, pz, per(3)

   nx = self%ni ; ny = self%nj ; nz = self%nk ; np = self%particle_number
   ni = self%bni ; nj = self%bnj ; nk = self%bnk ; nbl = self%bblocks ; ngc = self%ngc
   pk => self%p_gpu ; wa => self%work_a_gpu ; off => self%off_gpu
   e1 = self%vemin(1) ; e2 = self%vemin(2) ; e3 = self%vemin(3)
   d1 = self%vdx(1)   ; d2 = self%vdx(2)   ; d3 = self%vdx(3)
   s1 = self%shift(1) ; s2 = self%shift(2) ; s3 = self%shift(3)
   px = self%is_periodic(1) ; py = self%is_periodic(2) ; pz = self%is_periodic(3)
   rad = support_radius(self) ; radius = self%radius ; shp = self%shape ; sigma = self%sigma
   cutoff_limit = self%cutoff_sigma + 64._R8P*epsilon(self%cutoff_sigma)*max(1._R8P, abs(self%cutoff_sigma))
   !$acc parallel loop independent gang vector collapse(3) DEVICEVAR(wa)
   !$omp OMPLOOP collapse(3) DEVICEPTR(wa)
   do k=1, nz
   do j=1, ny
   do i=1, nx
      wa(1,i,j,k) = 0._R8P
   enddo
   enddo
   enddo
   ierr = 0_I4P
   if (np > 0_I4P) then
      q_pic => pic_fnl%q_pic_gpu
      !$acc parallel loop gang DEVICEVAR(q_pic,wa,pk) &
      !$acc& firstprivate(nx,ny,nz,rad,radius,shp,sigma,cutoff_limit,e1,e2,e3,d1,d2,d3,s1,s2,s3,px,py,pz) &
      !$acc& private(nn,per,w,xc,pw,dxyz,x,emin,ic,lo,rr,perr,charge_density) reduction(max:ierr)
      !$omp OMPLOOP DEVICEPTR(q_pic,wa,pk) &
      !$omp& firstprivate(nx,ny,nz,rad,radius,shp,sigma,cutoff_limit,e1,e2,e3,d1,d2,d3,s1,s2,s3,px,py,pz) &
      !$omp& private(nn,per,w,xc,pw,dxyz,x,emin,ic,lo,rr,perr,charge_density) reduction(max:ierr)
      do part=1, np
         perr = 0_I4P
         nn(1) = nx ; nn(2) = ny ; nn(3) = nz
         per(1) = px ; per(2) = py ; per(3) = pz
         dxyz(1) = d1 ; dxyz(2) = d2 ; dxyz(3) = d3
         emin(1) = e1 ; emin(2) = e2 ; emin(3) = e3
         x(1) = q_pic(part,1) - s1 ; x(2) = q_pic(part,2) - s2 ; x(3) = q_pic(part,3) - s3
         !$acc loop seq
         do d=1, 3
            ic(d) = ceiling((x(d) - emin(d)) / dxyz(d), kind=I4P)
            lo(d) = ic(d) - rad
            rr(d) = 2_I4P*rad + 1_I4P
            if (.not.per(d)) then
               if (lo(d) <= radius .or. nn(d)-(ic(d)+rad) <= radius) perr = ERR_FILTER
            endif
         enddo
         if (perr == 0_I4P) then
            !$acc loop seq
            do d=1, 3
               !$acc loop seq
               do m=1, rr(d)
                  xc(m) = emin(d) + (real(lo(d)+m-1_I4P, R8P) - 0.5_R8P)*dxyz(d)
               enddo
               call shape_weights_dev(shp=shp, filter=.false., sigma=sigma, cutoff_limit=cutoff_limit, &
                                      radius=rad, x=x(d), ic=ic(d)-lo(d)+1_I4P, span=rr(d), dx=dxyz(d), xc=xc, &
                                      w=w(:,d), perr=perr)
               ! P on the weights, cell lo-radius+l-1
               !$acc loop seq
               do l=1, rr(d)+2_I4P*radius
                  acc = 0._R8P
                  !$acc loop seq
                  do sm=-radius, radius
                     m = l - radius + sm
                     if (m < 1_I4P .or. m > rr(d)) cycle
                     acc = acc + pk(sm+8)*w(m,d)
                  enddo
                  pw(l,d) = acc
               enddo
            enddo
         endif
         if (perr /= 0_I4P) then
            ierr = max(ierr, perr)
         else
            charge_density = q_pic(part,7) / (d1*d2*d3)
            !$acc loop vector collapse(3) private(value,ii,jj,kk)
            do lk=1, rr(3)+2*radius
            do lj=1, rr(2)+2*radius
            do li=1, rr(1)+2*radius
               value = charge_density*pw(li,1)*pw(lj,2)*pw(lk,3)
               if (value /= 0._R8P) then
                  ii = lo(1) - radius + li - 1_I4P
                  jj = lo(2) - radius + lj - 1_I4P
                  kk = lo(3) - radius + lk - 1_I4P
                  if (px) ii = modulo(ii-1_I4P, nx) + 1_I4P
                  if (py) jj = modulo(jj-1_I4P, ny) + 1_I4P
                  if (pz) kk = modulo(kk-1_I4P, nz) + 1_I4P
                  if (ii >= 1_I4P .and. ii <= nx .and. jj >= 1_I4P .and. jj <= ny .and. kk >= 1_I4P .and. kk <= nz) then
                     !$acc atomic update
                     !$omp atomic update
                     wa(1,ii,jj,kk) = wa(1,ii,jj,kk) + value
                  endif
               endif
            enddo
            enddo
            enddo
         endif
      enddo
   endif
   call check_particle_error(ierr)
   call self%reduce_virtual(a=self%work_a_gpu)
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q,wa,off) &
   !$acc& firstprivate(nv,nx,ny,nz,px,py,pz) private(gi,gj,gk)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q,wa,off) firstprivate(nv,nx,ny,nz,px,py,pz) private(gi,gj,gk)
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do b=1, nbl
      gi = off(b,1) + i ; gj = off(b,2) + j ; gk = off(b,3) + k
      if (px) gi = modulo(gi-1_I4P, nx) + 1_I4P
      if (py) gj = modulo(gj-1_I4P, ny) + 1_I4P
      if (pz) gk = modulo(gk-1_I4P, nz) + 1_I4P
      if (gi >= 1 .and. gi <= nx .and. gj >= 1 .and. gj <= ny .and. gk >= 1 .and. gk <= nz) then
         q(b,i,j,k,nv) = wa(1,gi,gj,gk)
      else
         q(b,i,j,k,nv) = 0._R8P
      endif
   enddo
   enddo
   enddo
   enddo
   endsubroutine deposit_matched_charge_box

   subroutine coherent_field_weighting(self, pic_fnl, field_fnl, field, grid, pic_fields_gpu, q, q_pic_gpu, nv)
   !< Gather E and B with the transpose of the matched charge operator: P^T on the fields, then the unfiltered base shape
   !< (on the dual nodes x_c+dx/2 for the odd esirkepov_Q<N>), the adjoint of deposit_modified_charge. CPU twin:
   !< coherent_field_weighting of prism_pic_object.
   !<
   !< Multi-block, multi-rank: as the charge deposit, P^T runs on the virtual block. Each field component is copied from
   !< the real interiors (every virtual cell has one owner block on one rank, the others add zero), summed over the ranks,
   !< filtered, and scattered back to every real cell, ghosts included (interface ghosts get the neighbour values, physical
   !< ghosts 0, the zero extension of the single-block filter).
   class(prism_fnl_pic_conserving_object), intent(inout) :: self              !< Conserving current object.
   type(prism_fnl_pic_object),             intent(inout) :: pic_fnl           !< Device PIC helper.
   type(field_fnl_object),                 intent(in)    :: field_fnl         !< Device field helper.
   type(field_object),                     intent(in)    :: field             !< Host field.
   type(grid_object),                      intent(in)    :: grid              !< Grid.
   real(R8P),                              intent(inout) :: pic_fields_gpu(:,:) !< Fields at the particles [np,6].
   real(R8P),                              intent(in)    :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Real field.
   real(R8P),                              intent(in)    :: q_pic_gpu(:,:)    !< Gather particles [np,8].
   integer(I4P),                           intent(in)    :: nv                !< Field variables number.
   integer(I4P), pointer                                 :: off(:,:)
   real(R8P),    pointer                                 :: wa(:,:,:,:), wb(:,:,:,:), gq(:,:,:,:,:), q_dual(:,:)
   integer(I4P)                                          :: b, i, j, k, v, n, ni, nj, nk, nbl, ngc, nx, ny, nz, gi, gj, gk, np
   real(R8P)                                             :: s1, s2, s3
   logical                                               :: old_filter

   if (.not.self%coherent) call mpih%error_stop(msg=': coherent gather kernels not initialized in conserving_fnl')
   ni = self%bni ; nj = self%bnj ; nk = self%bnk ; nbl = self%bblocks ; ngc = self%ngc
   nx = self%ni ; ny = self%nj ; nz = self%nk ; np = self%particle_number
   off => self%off_gpu ; wa => self%work_a_gpu ; wb => self%work_b_gpu ; gq => self%gq_gpu
   do v=1, 6
      !$acc parallel loop independent gang vector collapse(3) DEVICEVAR(wa)
      !$omp OMPLOOP collapse(3) DEVICEPTR(wa)
      do k=1, nz
      do j=1, ny
      do i=1, nx
         wa(1,i,j,k) = 0._R8P
      enddo
      enddo
      enddo
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q,wa,off) firstprivate(v)
      !$omp OMPLOOP collapse(4) DEVICEPTR(q,wa,off) firstprivate(v)
      do k=1, nk
      do j=1, nj
      do i=1, ni
      do b=1, nbl
         wa(1,off(b,1)+i,off(b,2)+j,off(b,3)+k) = q(b,i,j,k,v)
      enddo
      enddo
      enddo
      enddo
      call self%reduce_virtual(a=self%work_a_gpu)
      call self%transposed_filter_axis(input=self%work_a_gpu, output=self%work_b_gpu, axis=1_I4P)
      call self%transposed_filter_axis(input=self%work_b_gpu, output=self%work_a_gpu, axis=2_I4P)
      call self%transposed_filter_axis(input=self%work_a_gpu, output=self%work_b_gpu, axis=3_I4P)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(gq,wb,off) firstprivate(v,nx,ny,nz) &
      !$acc& private(gi,gj,gk)
      !$omp OMPLOOP collapse(4) DEVICEPTR(gq,wb,off) firstprivate(v,nx,ny,nz) private(gi,gj,gk)
      do k=1-ngc, nk+ngc
      do j=1-ngc, nj+ngc
      do i=1-ngc, ni+ngc
      do b=1, nbl
         gi = off(b,1) + i ; gj = off(b,2) + j ; gk = off(b,3) + k
         if (gi >= 1 .and. gi <= nx .and. gj >= 1 .and. gj <= ny .and. gk >= 1 .and. gk <= nz) then
            gq(b,i,j,k,v) = wb(1,gi,gj,gk)
         else
            gq(b,i,j,k,v) = 0._R8P
         endif
      enddo
      enddo
      enddo
      enddo
   enddo
   old_filter = pic_fnl%filter_deposition
   pic_fnl%filter_deposition = .false.
   if (np > 0_I4P .and. any(self%gshift /= 0._R8P)) then
      ! odd filters: the base shape on the dual nodes, read as the particles at x-dx/2; then restore the real lookup
      q_dual => self%q_dual_gpu
      s1 = self%gshift(1) ; s2 = self%gshift(2) ; s3 = self%gshift(3)
      !$acc parallel loop independent gang vector collapse(2) DEVICEVAR(q_dual,q_pic_gpu) firstprivate(s1,s2,s3)
      !$omp OMPLOOP collapse(2) DEVICEPTR(q_dual,q_pic_gpu) firstprivate(s1,s2,s3)
      do v=1, 8
      do n=1, np
         select case(v)
         case(1)      ; q_dual(n,v) = q_pic_gpu(n,v) - s1
         case(2)      ; q_dual(n,v) = q_pic_gpu(n,v) - s2
         case(3)      ; q_dual(n,v) = q_pic_gpu(n,v) - s3
         case default ; q_dual(n,v) = q_pic_gpu(n,v)
         endselect
      enddo
      enddo
      call pic_fnl%particle_cartesian_grid_index_dev(field_fnl=field_fnl, field=field, grid=grid, q_pic_gpu=q_dual)
      call pic_fnl%field_weighting_dev(field_fnl=field_fnl, field=field, grid=grid, pic_fields_gpu=pic_fields_gpu, &
                                       q_gpu=gq, q_pic_gpu=q_dual, nv=nv)
      call pic_fnl%particle_cartesian_grid_index_dev(field_fnl=field_fnl, field=field, grid=grid, q_pic_gpu=q_pic_gpu)
   else
      call pic_fnl%field_weighting_dev(field_fnl=field_fnl, field=field, grid=grid, pic_fields_gpu=pic_fields_gpu, &
                                       q_gpu=gq, q_pic_gpu=q_pic_gpu, nv=nv)
   endif
   pic_fnl%filter_deposition = old_filter
   endsubroutine coherent_field_weighting

   subroutine copy_field(self, src, dst)
   !< dst = src, full real-block q-like device field (ghosts included).
   class(prism_fnl_pic_conserving_object), intent(in)    :: self !< Conserving current object.
   real(R8P),                              intent(in)    :: src(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Source.
   real(R8P),                              intent(inout) :: dst(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Destination.
   integer(I4P)                                          :: b, i, j, k, v, ngc, ni, nj, nk, nbl, nvv

   ngc = self%ngc ; ni = self%bni ; nj = self%bnj ; nk = self%bnk ; nbl = self%bblocks ; nvv = size(src, dim=5)
   !$acc parallel loop independent gang vector collapse(5) DEVICEVAR(src,dst)
   !$omp OMPLOOP collapse(5) DEVICEPTR(src,dst)
   do v=1, nvv
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do b=1, nbl
      dst(b,i,j,k,v) = src(b,i,j,k,v)
   enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine copy_field

   subroutine copy_particles(self, src, dst)
   !< dst = src, device particle arrays [np,8].
   class(prism_fnl_pic_conserving_object), intent(in)    :: self     !< Conserving current object.
   real(R8P),                              intent(in)    :: src(1:,1:) !< Source.
   real(R8P),                              intent(inout) :: dst(1:,1:) !< Destination.
   integer(I4P)                                          :: n, v, np

   np = self%particle_number
   if (np == 0_I4P) return
   !$acc parallel loop independent gang vector collapse(2) DEVICEVAR(src,dst)
   !$omp OMPLOOP collapse(2) DEVICEPTR(src,dst)
   do v=1, 8
   do n=1, np
      dst(n,v) = src(n,v)
   enddo
   enddo
   endsubroutine copy_particles

   subroutine build_stage_target(self, q_pic, q_pic_rk, s, coeff, dt)
   !< q_next(1:6) = q_pic(1:6) + dt sum_{r<=s} coeff(r) k_r(1:6), q_next(7:8) = q_pic(7:8).
   !< CPU twin: the q_pic_next loop of integrate_rk_ssp_pic_charge_conserving (same r order).
   class(prism_fnl_pic_conserving_object), intent(inout) :: self             !< Conserving current object.
   real(R8P),                              intent(in)    :: q_pic(1:,1:)     !< Particles at step start [np,8].
   real(R8P),                              intent(in)    :: q_pic_rk(1:,1:,1:) !< Stage right-hand sides [np,8,nrk].
   integer(I4P),                           intent(in)    :: s                !< Current stage.
   real(R8P),                              intent(in)    :: coeff(1:)        !< Stage coefficients (1:s used).
   real(R8P),                              intent(in)    :: dt               !< Time step.
   real(R8P), pointer                                    :: q_next(:,:)
   real(R8P)                                             :: c(1:8), val
   integer(I4P)                                          :: n, v, r, np

   np = self%particle_number
   if (np == 0_I4P) return
   q_next => self%q_next_gpu
   c = 0._R8P ; c(1:s) = coeff(1:s)
   !$acc parallel loop independent gang vector collapse(2) DEVICEVAR(q_pic,q_pic_rk,q_next) copyin(c) private(val)
   !$omp OMPLOOP collapse(2) DEVICEPTR(q_pic,q_pic_rk,q_next) map(to:c) private(val)
   do v=1, 8
   do n=1, np
      val = q_pic(n,v)
      if (v <= 6) then
         !$acc loop seq
         do r=1, s
            val = val + dt * c(r) * q_pic_rk(n,v,r)
         enddo
      endif
      q_next(n,v) = val
   enddo
   enddo
   endsubroutine build_stage_target

   subroutine build_centered_trajectory(self, q_pic, dt)
   !< q_ref = q_pic - dt/2 v, q_next = q_pic + dt/2 v (positions only). CPU twin: initialize_pic_conserving_current_time_zero.
   class(prism_fnl_pic_conserving_object), intent(inout) :: self         !< Conserving current object.
   real(R8P),                              intent(in)    :: q_pic(1:,1:) !< Particles [np,8].
   real(R8P),                              intent(in)    :: dt           !< Virtual time step.
   real(R8P), pointer                                    :: q_minus(:,:), q_plus(:,:)
   integer(I4P)                                          :: n, v, np

   np = self%particle_number
   if (np == 0_I4P) return
   q_minus => self%q_ref_gpu
   q_plus  => self%q_next_gpu
   !$acc parallel loop independent gang vector collapse(2) DEVICEVAR(q_pic,q_minus,q_plus)
   !$omp OMPLOOP collapse(2) DEVICEPTR(q_pic,q_minus,q_plus)
   do v=1, 8
   do n=1, np
      if (v <= 3) then
         q_minus(n,v) = q_pic(n,v) - 0.5_R8P*dt*q_pic(n,v+3)
         q_plus(n,v)  = q_pic(n,v) + 0.5_R8P*dt*q_pic(n,v+3)
      else
         q_minus(n,v) = q_pic(n,v)
         q_plus(n,v)  = q_pic(n,v)
      endif
   enddo
   enddo
   endsubroutine build_centered_trajectory

   subroutine esirkepov_displacement(self, field_fnl, q_ref, q_to, h, p_first, p_last)
   !< Symmetric Esirkepov density decomposition. CPU twin: compute_pic_esirkepov_displacement.
   !<
   !< Parallel mapping: one gang per particle, one vector lane per cell of the particle's trajectory box. The six 1D weight
   !< vectors (3 directions x 2 endpoints) are computed once per gang into gang-private arrays (CPU recomputes them inside
   !< the triple loop for the non-Gaussian/non-quartic shapes: same values, fewer evaluations). Different particles may
   !< overlap, hence the atomic updates.
   !<
   !< Cells are virtual-block (global) indices: a trajectory box may straddle any number of real blocks, only the domain
   !< boundary (and MAX_SPAN) limits it.
   class(prism_fnl_pic_conserving_object), intent(in)    :: self           !< Conserving current object.
   type(field_fnl_object),                 intent(in)    :: field_fnl      !< Device field helper (unused: virtual geometry).
   real(R8P),                              intent(in)    :: q_ref(1:,1:)   !< Trajectory start [np,8].
   real(R8P),                              intent(in)    :: q_to(1:,1:)    !< Trajectory end [np,8].
   real(R8P),                              intent(inout) :: h(1:,1:,1:,1:,1:) !< Directional charge change [1,NX,NY,NZ,3].
   integer(I4P),                           intent(in)    :: p_first        !< First particle.
   integer(I4P),                           intent(in)    :: p_last         !< Last particle.
   real(R8P), pointer                                    :: vx_gpu(:), vy_gpu(:), vz_gpu(:)
   real(R8P)                                             :: wref(MAX_SPAN,3), wto(MAX_SPAN,3), xc(MAX_SPAN)
   real(R8P)                                             :: dxyz(3), xr(3), xt(3), cutoff_limit, charge_density
   real(R8P)                                             :: a0, a1, da, b0, b1, db, c0, c1, dc, sigma
   real(R8P)                                             :: e1, e2, e3, d1, d2, d3, s1, s2, s3
   integer(I4P)                                          :: p, d, m, i, j, k, perr, ierr
   integer(I4P)                                          :: cr(3), ct(3), lo(3), hi(3), rr(3), rad, rmod
   integer(I4P)                                          :: ni, nj, nk, shp, nn(3)
   logical                                               :: modified, filter_esk

   call zero_h(self=self, h=h)
   if (p_last < p_first) return
   vx_gpu => self%vx_gpu ; vy_gpu => self%vy_gpu ; vz_gpu => self%vz_gpu
   ni = self%ni ; nj = self%nj ; nk = self%nk ; shp = self%shape
   e1 = self%vemin(1) ; e2 = self%vemin(2) ; e3 = self%vemin(3)
   d1 = self%vdx(1)   ; d2 = self%vdx(2)   ; d3 = self%vdx(3)
   s1 = self%shift(1) ; s2 = self%shift(2) ; s3 = self%shift(3)
   modified   = self%solver == CONSERVING_SOLVER_MODIFIED
   filter_esk = self%filter_deposition .and. (.not. modified)
   rad  = support_radius(self) ; if (filter_esk) rad = rad + 1_I4P
   rmod = 0_I4P ; if (modified) rmod = self%radius
   sigma = self%sigma
   cutoff_limit = self%cutoff_sigma + 64._R8P*epsilon(self%cutoff_sigma)*max(1._R8P, abs(self%cutoff_sigma))
   ierr = 0_I4P

   !$acc parallel loop gang DEVICEVAR(q_ref,q_to,h,vx_gpu,vy_gpu,vz_gpu) &
   !$acc& firstprivate(ni,nj,nk,rad,rmod,shp,modified,filter_esk,sigma,cutoff_limit,e1,e2,e3,d1,d2,d3,s1,s2,s3) &
   !$acc& private(nn,wref,wto,xc,dxyz,xr,xt,cr,ct,lo,hi,rr,perr,charge_density) reduction(max:ierr)
   !$omp OMPLOOP DEVICEPTR(q_ref,q_to,h,vx_gpu,vy_gpu,vz_gpu) &
   !$omp& firstprivate(ni,nj,nk,rad,rmod,shp,modified,filter_esk,sigma,cutoff_limit,e1,e2,e3,d1,d2,d3,s1,s2,s3) &
   !$omp& private(nn,wref,wto,xc,dxyz,xr,xt,cr,ct,lo,hi,rr,perr,charge_density) reduction(max:ierr)
   do p=p_first, p_last
      perr = 0_I4P
      nn(1) = ni ; nn(2) = nj ; nn(3) = nk
      dxyz(1) = d1 ; dxyz(2) = d2 ; dxyz(3) = d3
      ! base-shape positions (shifted by dx/2 on the dual nodes of the odd esirkepov_Q filters)
      xr(1) = q_ref(p,1) - s1 ; xr(2) = q_ref(p,2) - s2 ; xr(3) = q_ref(p,3) - s3
      xt(1) = q_to(p,1)  - s1 ; xt(2) = q_to(p,2)  - s2 ; xt(3) = q_to(p,3)  - s3
      ! find_pic_position_cell for both endpoints, virtual-block cells
      cr(1) = ceiling((xr(1) - e1) / d1, kind=I4P) ; ct(1) = ceiling((xt(1) - e1) / d1, kind=I4P)
      cr(2) = ceiling((xr(2) - e2) / d2, kind=I4P) ; ct(2) = ceiling((xt(2) - e2) / d2, kind=I4P)
      cr(3) = ceiling((xr(3) - e3) / d3, kind=I4P) ; ct(3) = ceiling((xt(3) - e3) / d3, kind=I4P)
      !$acc loop seq
      do d=1, 3
         if (cr(d) < 1_I4P .or. cr(d) > nn(d) .or. ct(d) < 1_I4P .or. ct(d) > nn(d)) perr = ERR_OUTSIDE
      enddo
      if (perr == 0_I4P .and. q_ref(p,7) /= q_to(p,7)) perr = ERR_CHARGE
      if (perr == 0_I4P) then
         !$acc loop seq
         do d=1, 3
            lo(d) = min(cr(d),ct(d)) - rad
            hi(d) = max(cr(d),ct(d)) + rad
            if (lo(d) < 1_I4P .or. hi(d) > nn(d)) perr = ERR_BOUNDARY
         enddo
      endif
      if (perr == 0_I4P .and. modified) then
         if (min(lo(1),lo(2),lo(3)) <= rmod .or. nn(1)-hi(1) <= rmod .or. nn(2)-hi(2) <= rmod .or. nn(3)-hi(3) <= rmod) &
            perr = ERR_FILTER
      endif
      if (perr == 0_I4P) then
         !$acc loop seq
         do d=1, 3
            rr(d) = hi(d) - lo(d) + 1_I4P
            if (rr(d) > MAX_SPAN) perr = ERR_SPAN
         enddo
      endif
      if (perr == 0_I4P) then
         ! 1D weights of both endpoints, gang-private
         !$acc loop seq
         do d=1, 3
            !$acc loop seq
            do m=1, rr(d)
               select case(d)
               case(1) ; xc(m) = vx_gpu(lo(d)+m-1)
               case(2) ; xc(m) = vy_gpu(lo(d)+m-1)
               case(3) ; xc(m) = vz_gpu(lo(d)+m-1)
               endselect
            enddo
            call shape_weights_dev(shp=shp, filter=filter_esk, sigma=sigma, cutoff_limit=cutoff_limit, &
                                   radius=rad, x=xr(d), ic=cr(d)-lo(d)+1_I4P, span=rr(d), dx=dxyz(d), xc=xc, &
                                   w=wref(:,d), perr=perr)
            call shape_weights_dev(shp=shp, filter=filter_esk, sigma=sigma, cutoff_limit=cutoff_limit, &
                                   radius=rad, x=xt(d), ic=ct(d)-lo(d)+1_I4P, span=rr(d), dx=dxyz(d), xc=xc, &
                                   w=wto(:,d), perr=perr)
         enddo
      endif
      if (perr /= 0_I4P) then
         ierr = max(ierr, perr)
      else
         charge_density = q_ref(p,7) / (dxyz(1)*dxyz(2)*dxyz(3))
         !$acc loop vector collapse(3) private(a0,a1,da,b0,b1,db,c0,c1,dc)
         do k=lo(3), hi(3)
         do j=lo(2), hi(2)
         do i=lo(1), hi(1)
            c0 = wref(k-lo(3)+1,3) ; c1 = wto(k-lo(3)+1,3) ; dc = c1-c0
            b0 = wref(j-lo(2)+1,2) ; b1 = wto(j-lo(2)+1,2) ; db = b1-b0
            a0 = wref(i-lo(1)+1,1) ; a1 = wto(i-lo(1)+1,1) ; da = a1-a0
            !$acc atomic update
            !$omp atomic update
            h(1,i,j,k,1) = h(1,i,j,k,1) + charge_density*da*(b0*c0 + 0.5_R8P*(db*c0+b0*dc) + db*dc/3._R8P)
            !$acc atomic update
            !$omp atomic update
            h(1,i,j,k,2) = h(1,i,j,k,2) + charge_density*db*(a0*c0 + 0.5_R8P*(da*c0+a0*dc) + da*dc/3._R8P)
            !$acc atomic update
            !$omp atomic update
            h(1,i,j,k,3) = h(1,i,j,k,3) + charge_density*dc*(a0*b0 + 0.5_R8P*(da*b0+a0*db) + da*db/3._R8P)
         enddo
         enddo
         enddo
      endif
   enddo
   call check_particle_error(ierr)
   endsubroutine esirkepov_displacement

   subroutine deposit_matched_current_box(self, q_ref, q_to, h, p_first, p_last, dt)
   !< Matched current (totorica, esirkepov_Q<N>) of the trajectories q_ref -> q_to, built particle by particle on the box of
   !< each trajectory and folded into the virtual block along the periodic axes: J_d = Q_d P_e1 P_e2 F_d^E. CPU twin:
   !< deposit_matched_current (the separable factors are explained there). Used with periodic faces only: h receives the
   !< stage current in place of the stage source, combined over the stages by combine_stage_source as the sources are (the
   !< solve is linear) and written to the real blocks by store_box_current.
   !<
   !< Parallel mapping: one gang per particle. The 1D factors (weights of both endpoints, face prefix sum + Q along d, P on
   !< the transverse weights) are computed sequentially into gang-private arrays, one column per direction and written
   !< once (the zero-flux closure past the last weight change is applied when the prefix sums are read: nvfortran 25.11
   !< -fast can lose a second store into a private array). The deposit is one vector lane per cell of the widened box,
   !< atomic (particles overlap). Cells are unwrapped, center emin+(i-1/2)dx, and folded modulo n along the periodic axes.
   class(prism_fnl_pic_conserving_object), intent(inout) :: self              !< Conserving current object.
   real(R8P),                              intent(in)    :: q_ref(1:,1:)      !< Trajectory start [np,8].
   real(R8P),                              intent(in)    :: q_to(1:,1:)       !< Trajectory end [np,8].
   real(R8P),                              intent(inout) :: h(1:,1:,1:,1:,1:) !< Current [1,NX,NY,NZ,3].
   integer(I4P),                           intent(in)    :: p_first           !< First particle.
   integer(I4P),                           intent(in)    :: p_last            !< Last particle.
   real(R8P),                              intent(in)    :: dt                !< Time step.
   real(R8P), pointer                                    :: pk(:), qw(:)
   real(R8P)                                             :: wref(MAX_SPAN,3), wto(MAX_SPAN,3), xc(MAX_SPAN)
   real(R8P)                                             :: g(0:MAX_SPAN,3), lq(LQ_SPAN,3)
   real(R8P)                                             :: pt0(LP_SPAN,3), pdt(LP_SPAN,3), pu(LP_SPAN,3), pv(LP_SPAN,3)
   real(R8P)                                             :: dxyz(3), xr(3), xt(3), emin(3)
   real(R8P)                                             :: cutoff_limit, charge_density, sigma, res, acc, acd, value, w0, dw
   real(R8P)                                             :: e1, e2, e3, d1, d2, d3, s1, s2, s3
   integer(I4P)                                          :: part, d, a, b, m, l, t, sm, li, lj, lk, perr, ierr
   integer(I4P)                                          :: rad, radius, qfirst, nq, shp, ni, nj, nk, nli, nlj, nlk
   integer(I4P)                                          :: ci, cj, ck, ii, jj, kk
   integer(I4P)                                          :: cr(3), ct(3), lo(3), hi(3), rr(3), nn(3), last(3)
   logical                                               :: px, py, pz, per(3)

   call zero_h(self=self, h=h)
   if (p_last < p_first) return
   pk => self%p_gpu ; qw => self%qw_gpu
   ni = self%ni ; nj = self%nj ; nk = self%nk ; shp = self%shape
   e1 = self%vemin(1) ; e2 = self%vemin(2) ; e3 = self%vemin(3)
   d1 = self%vdx(1)   ; d2 = self%vdx(2)   ; d3 = self%vdx(3)
   s1 = self%shift(1) ; s2 = self%shift(2) ; s3 = self%shift(3)
   px = self%is_periodic(1) ; py = self%is_periodic(2) ; pz = self%is_periodic(3)
   rad = support_radius(self) ; radius = self%radius ; qfirst = self%qfirst ; nq = self%nq
   sigma = self%sigma
   cutoff_limit = self%cutoff_sigma + 64._R8P*epsilon(self%cutoff_sigma)*max(1._R8P, abs(self%cutoff_sigma))
   ierr = 0_I4P ; res = 0._R8P

   !$acc parallel loop gang DEVICEVAR(q_ref,q_to,h,pk,qw) &
   !$acc& firstprivate(ni,nj,nk,rad,radius,qfirst,nq,shp,sigma,cutoff_limit,e1,e2,e3,d1,d2,d3,s1,s2,s3,px,py,pz,dt) &
   !$acc& private(nn,per,wref,wto,xc,g,lq,pt0,pdt,pu,pv,dxyz,xr,xt,emin,cr,ct,lo,hi,rr,last,perr,charge_density) &
   !$acc& private(a,b,nli,nlj,nlk) reduction(max:ierr,res)
   !$omp OMPLOOP DEVICEPTR(q_ref,q_to,h,pk,qw) &
   !$omp& firstprivate(ni,nj,nk,rad,radius,qfirst,nq,shp,sigma,cutoff_limit,e1,e2,e3,d1,d2,d3,s1,s2,s3,px,py,pz,dt) &
   !$omp& private(nn,per,wref,wto,xc,g,lq,pt0,pdt,pu,pv,dxyz,xr,xt,emin,cr,ct,lo,hi,rr,last,perr,charge_density) &
   !$omp& private(a,b,nli,nlj,nlk) reduction(max:ierr,res)
   do part=p_first, p_last
      perr = 0_I4P
      nn(1) = ni ; nn(2) = nj ; nn(3) = nk
      per(1) = px ; per(2) = py ; per(3) = pz
      dxyz(1) = d1 ; dxyz(2) = d2 ; dxyz(3) = d3
      emin(1) = e1 ; emin(2) = e2 ; emin(3) = e3
      xr(1) = q_ref(part,1) - s1 ; xr(2) = q_ref(part,2) - s2 ; xr(3) = q_ref(part,3) - s3
      xt(1) = q_to(part,1)  - s1 ; xt(2) = q_to(part,2)  - s2 ; xt(3) = q_to(part,3)  - s3
      !$acc loop seq
      do d=1, 3
         cr(d) = ceiling((xr(d) - emin(d)) / dxyz(d), kind=I4P)
         ct(d) = ceiling((xt(d) - emin(d)) / dxyz(d), kind=I4P)
         if (.not.per(d)) then
            if (cr(d) < 1_I4P .or. cr(d) > nn(d) .or. ct(d) < 1_I4P .or. ct(d) > nn(d)) perr = ERR_OUTSIDE
         endif
      enddo
      if (perr == 0_I4P .and. q_ref(part,7) /= q_to(part,7)) perr = ERR_CHARGE
      if (perr == 0_I4P) then
         !$acc loop seq
         do d=1, 3
            lo(d) = min(cr(d),ct(d)) - rad
            hi(d) = max(cr(d),ct(d)) + rad
            rr(d) = hi(d) - lo(d) + 1_I4P
            if (rr(d) > MAX_SPAN) perr = ERR_SPAN
            if (.not.per(d)) then
               if (lo(d) <= radius .or. nn(d)-hi(d) <= radius) perr = ERR_FILTER
            endif
         enddo
      endif
      if (perr == 0_I4P) then
         !$acc loop seq
         do d=1, 3
            !$acc loop seq
            do m=1, rr(d)
               xc(m) = emin(d) + (real(lo(d)+m-1_I4P, R8P) - 0.5_R8P)*dxyz(d)
            enddo
            call shape_weights_dev(shp=shp, filter=.false., sigma=sigma, cutoff_limit=cutoff_limit, &
                                   radius=rad, x=xr(d), ic=cr(d)-lo(d)+1_I4P, span=rr(d), dx=dxyz(d), xc=xc, &
                                   w=wref(:,d), perr=perr)
            call shape_weights_dev(shp=shp, filter=.false., sigma=sigma, cutoff_limit=cutoff_limit, &
                                   radius=rad, x=xt(d), ic=ct(d)-lo(d)+1_I4P, span=rr(d), dx=dxyz(d), xc=xc, &
                                   w=wto(:,d), perr=perr)
         enddo
      endif
      if (perr /= 0_I4P) then
         ierr = max(ierr, perr)
      else
         charge_density = q_ref(part,7) / (d1*d2*d3)
         ! longitudinal factors: face prefix sum of the weight change (faces lo-1..hi, index m = face-(lo-1)), then Q
         !$acc loop seq
         do d=1, 3
            g(0,d) = 0._R8P
            last(d) = 0_I4P
            !$acc loop seq
            do m=1, rr(d)
               g(m,d) = g(m-1,d) - dxyz(d)*charge_density*(wto(m,d) - wref(m,d))/dt
               if (wto(m,d) /= wref(m,d)) last(d) = m
            enddo
            res = max(res, abs(g(rr(d),d)))
            ! cell index of lq(l,d): lo(d)-1-qfirst-(nq-1)+l-1; its faces m = l-nq+t, zero past the last weight change
            !$acc loop seq
            do l=1, rr(d)+nq
               acc = 0._R8P
               !$acc loop seq
               do t=0, nq-1
                  m = l - nq + t
                  if (m < 0_I4P .or. m > rr(d)) cycle
                  if (last(d) > 0_I4P .and. m >= last(d)) cycle
                  acc = acc + qw(t+1)*g(m,d)
               enddo
               lq(l,d) = acc
            enddo
         enddo
         ! transverse factors of direction d: P on t0, dt (axis a) and on s0+ds/2, s0/2+ds/3 (axis b); cell lo-radius+l-1
         !$acc loop seq
         do d=1, 3
            a = 1_I4P + merge(1_I4P, 0_I4P, d == 1_I4P)
            b = 3_I4P - merge(1_I4P, 0_I4P, d == 3_I4P)
            !$acc loop seq
            do l=1, rr(a)+2_I4P*radius
               acc = 0._R8P ; acd = 0._R8P
               !$acc loop seq
               do sm=-radius, radius
                  m = l - radius + sm
                  if (m < 1_I4P .or. m > rr(a)) cycle
                  acc = acc + pk(sm+8)*wref(m,a)
                  acd = acd + pk(sm+8)*(wto(m,a) - wref(m,a))
               enddo
               pt0(l,d) = acc ; pdt(l,d) = acd
            enddo
            !$acc loop seq
            do l=1, rr(b)+2_I4P*radius
               acc = 0._R8P ; acd = 0._R8P
               !$acc loop seq
               do sm=-radius, radius
                  m = l - radius + sm
                  if (m < 1_I4P .or. m > rr(b)) cycle
                  w0 = wref(m,b) ; dw = wto(m,b) - wref(m,b)
                  acc = acc + pk(sm+8)*(w0 + 0.5_R8P*dw)
                  acd = acd + pk(sm+8)*(0.5_R8P*w0 + dw/3._R8P)
               enddo
               pu(l,d) = acc ; pv(l,d) = acd
            enddo
         enddo
         ! deposit, one vector lane per cell of the widened box, folded along the periodic axes
         !$acc loop seq
         do d=1, 3
            a = 1_I4P + merge(1_I4P, 0_I4P, d == 1_I4P)
            b = 3_I4P - merge(1_I4P, 0_I4P, d == 3_I4P)
            nli = rr(d) + nq ; nlj = rr(a) + 2_I4P*radius ; nlk = rr(b) + 2_I4P*radius
            !$acc loop vector collapse(3) private(value,ci,cj,ck,ii,jj,kk)
            do lk=1, nlk
            do lj=1, nlj
            do li=1, nli
               value = lq(li,d)*(pt0(lj,d)*pu(lk,d) + pdt(lj,d)*pv(lk,d))
               if (value /= 0._R8P) then
                  ci = lo(d) - 1_I4P - qfirst - (nq - 1_I4P) + li - 1_I4P
                  cj = lo(a) - radius + lj - 1_I4P
                  ck = lo(b) - radius + lk - 1_I4P
                  if (d == 1_I4P) then
                     ii = ci ; jj = cj ; kk = ck
                  elseif (d == 2_I4P) then
                     ii = cj ; jj = ci ; kk = ck
                  else
                     ii = cj ; jj = ck ; kk = ci
                  endif
                  if (px) ii = modulo(ii-1_I4P, ni) + 1_I4P
                  if (py) jj = modulo(jj-1_I4P, nj) + 1_I4P
                  if (pz) kk = modulo(kk-1_I4P, nk) + 1_I4P
                  if (ii >= 1_I4P .and. ii <= ni .and. jj >= 1_I4P .and. jj <= nj .and. kk >= 1_I4P .and. kk <= nk) then
                     !$acc atomic update
                     !$omp atomic update
                     h(1,ii,jj,kk,d) = h(1,ii,jj,kk,d) + value
                  endif
               endif
            enddo
            enddo
            enddo
         enddo
      endif
   enddo
   call check_particle_error(ierr)
   self%residual_max = max(self%residual_max, res)
   endsubroutine deposit_matched_current_box

   subroutine store_box_current(self, q, var_jx, hq)
   !< Real q(J) = virtual stage current of the box solver (deposit_matched_current_box + combine_stage_source), ghosts
   !< included: the virtual ghosts are the periodic images along the periodic axes and the odd mirror along the others
   !< (impose_virtual_ghosts_periodic), then scatter_current. CPU twin: store_matched_current.
   class(prism_fnl_pic_conserving_object), intent(inout) :: self               !< Conserving current object.
   real(R8P),                              intent(inout) :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Real field.
   integer(I4P),                           intent(in)    :: var_jx             !< Jx slot in q.
   real(R8P),                              intent(in)    :: hq(1:,1:,1:,1:,1:) !< Virtual current [1,NX,NY,NZ,3].
   real(R8P), pointer                                    :: jv(:,:,:,:,:)
   integer(I4P)                                          :: i, j, k, d, ni, nj, nk

   ni = self%ni ; nj = self%nj ; nk = self%nk
   jv => self%jv_gpu
   call zero_current(self=self, q=jv, var_jx=1_I4P)
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(jv,hq)
   !$omp OMPLOOP collapse(4) DEVICEPTR(jv,hq)
   do d=1, 3
   do k=1, nk
   do j=1, nj
   do i=1, ni
      jv(1,i,j,k,d) = hq(1,i,j,k,d)
   enddo
   enddo
   enddo
   enddo
   call self%impose_virtual_ghosts_periodic(q=jv, var_jx=1_I4P)
   call self%scatter_current(q=q, var_jx=var_jx)
   endsubroutine store_box_current

   subroutine combine_stage_source(self, h, s, coeff)
   !< H_s = (H_s - sum_{r<s} coeff(r) H_r) / coeff(s). CPU twin: the h_q SSP combination of the conserving integrator.
   class(prism_fnl_pic_conserving_object), intent(in)    :: self                 !< Conserving current object.
   real(R8P),                              intent(inout) :: h(1:,1:,1:,1:,1:,1:) !< Stage sources [nb,ni,nj,nk,3,nrk].
   integer(I4P),                           intent(in)    :: s                    !< Current stage.
   real(R8P),                              intent(in)    :: coeff(1:)            !< Stage coefficients (1:s used).
   real(R8P)                                             :: c(1:8), val
   integer(I4P)                                          :: b, i, j, k, d, r, ni, nj, nk, nbl

   if (abs(coeff(s)) <= tiny(1._R8P)) call mpih%error_stop(msg=': zero SSP coefficient in charge-conserving PIC current')
   ni = self%ni ; nj = self%nj ; nk = self%nk ; nbl = self%blocks_number
   c = 0._R8P ; c(1:s) = coeff(1:s)
   !$acc parallel loop independent gang vector collapse(5) DEVICEVAR(h) copyin(c) private(val)
   !$omp OMPLOOP collapse(5) DEVICEPTR(h) map(to:c) private(val)
   do d=1, 3
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, nbl
      val = h(b,i,j,k,d,s)
      !$acc loop seq
      do r=1, s-1
         val = val - c(r) * h(b,i,j,k,d,r)
      enddo
      h(b,i,j,k,d,s) = val / c(s)
   enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine combine_stage_source

   subroutine solve_dispatch(self, field_fnl, q, var_jx, hq, dt)
   !< Select the configured current reconstruction. CPU twin: solve_pic_charge_conserving_current_dispatch.
   !< The solve runs on the virtual block (jv_gpu), then J is scattered to the real blocks, ghosts included.
   class(prism_fnl_pic_conserving_object), intent(inout) :: self      !< Conserving current object.
   type(field_fnl_object),                 intent(in)    :: field_fnl !< Device field helper.
   real(R8P),                              intent(inout) :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Real target field.
   integer(I4P),                           intent(in)    :: var_jx    !< Jx slot in q (Jy, Jz follow).
   real(R8P),                              intent(in)    :: hq(1:,1:,1:,1:,1:) !< Directional source [1,NX,NY,NZ,3].
   real(R8P),                              intent(in)    :: dt        !< Time step.

   select case(self%solver)
   case(CONSERVING_SOLVER_LONDRILLO)
      call self%solve_lines(field_fnl=field_fnl, q=self%jv_gpu, var_jx=1_I4P, hq=hq, dt=dt)
   case(CONSERVING_SOLVER_MODIFIED)
      call self%solve_modified(field_fnl=field_fnl, q=self%jv_gpu, var_jx=1_I4P, hq=hq, dt=dt)
   endselect
   call self%scatter_current(q=q, var_jx=var_jx)
   endsubroutine solve_dispatch

   subroutine solve_lines(self, field_fnl, q, var_jx, hq, dt)
   !< Londrillo line solves along x, y, z: Esirkepov face prefix sum + inverse reconstruction.
   !< CPU twins: solve_pic_charge_conserving_current_londrillo, build_esirkepov_face_line, update_esirkepov_face_residual,
   !< solve_factored_matrix_pivot.
   !<
   !< Parallel mapping: one thread per grid line (b, t1, t2), sequential LU substitution along the line against the cached
   !< factored matrix of its block (all lanes of a block read the same matrix entry: broadcast). The line is staged in a
   !< line-major scratch `line(line + nl*(m-1))`, so at each substitution step adjacent lanes touch adjacent words
   !< (coalesced), instead of the strided q(b,i,j,k,v) layout (stride nb*(ni+2ngc) for x-lines). The CPU flux array is not
   !< stored: the residual pass rebuilds the prefix sum in the same order.
   class(prism_fnl_pic_conserving_object), intent(inout) :: self      !< Conserving current object.
   type(field_fnl_object),                 intent(in)    :: field_fnl !< Device field helper.
   real(R8P),                              intent(inout) :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Target field.
   integer(I4P),                           intent(in)    :: var_jx    !< Jx slot in q.
   real(R8P),                              intent(in)    :: hq(1:,1:,1:,1:,1:) !< Directional source [nb,ni,nj,nk,3].
   real(R8P),                              intent(in)    :: dt        !< Time step.
   real(R8P), pointer                                    :: mat(:,:,:), line(:), fv1(:)
   integer(I4P), pointer                                 :: piv(:,:)
   real(R8P)                                             :: res, flux, tmp, hdx, recon, coeff, sgn
   integer(I4P)                                          :: dir, n, n1, n2, nl, nbl, b, t1, t2, lid, m, c, pv
   integer(I4P)                                          :: i, j, k, face, mm, col, hs, var

   call zero_current(self=self, q=q, var_jx=var_jx)
   hs = self%hs ; nbl = self%blocks_number
   line => self%line_gpu ; fv1 => self%fv1_gpu
   do dir=1, 3
      hdx = self%vdx(dir)
      select case(dir)
      case(1) ; n = self%ni ; n1 = self%nj ; n2 = self%nk ; mat => self%mat_x_gpu ; piv => self%piv_x_gpu
      case(2) ; n = self%nj ; n1 = self%ni ; n2 = self%nk ; mat => self%mat_y_gpu ; piv => self%piv_y_gpu
      case(3) ; n = self%nk ; n1 = self%ni ; n2 = self%nj ; mat => self%mat_z_gpu ; piv => self%piv_z_gpu
      endselect
      nl = nbl*n1*n2
      var = var_jx + dir - 1_I4P
      res = 0._R8P
      !$acc parallel loop independent gang vector collapse(3) DEVICEVAR(q,hq,mat,piv,line,fv1) &
      !$acc& firstprivate(dir,n,n1,nl,hs,var,dt,hdx) &
      !$acc& private(lid,flux,tmp,recon,coeff,sgn,m,c,pv,i,j,k,face,mm,col) reduction(max:res)
      !$omp OMPLOOP collapse(3) DEVICEPTR(q,hq,mat,piv,line,fv1) &
      !$omp& firstprivate(dir,n,n1,nl,hs,var,dt,hdx) &
      !$omp& private(lid,flux,tmp,recon,coeff,sgn,m,c,pv,i,j,k,face,mm,col) reduction(max:res)
      do t2=1, n2
      do t1=1, n1
      do b=1, nbl
         lid = b + nbl*((t1-1) + n1*(t2-1))
         ! right-hand side
         flux = 0._R8P
         !$acc loop seq
         do m=1, n
            call line_cell(dir, m, t1, t2, i, j, k)
            flux = flux - hdx*hq(b,i,j,k,dir)/dt
            if (m < n) line(lid+nl*(m-1)) = flux
         enddo
         line(lid+nl*(n-1)) = 0._R8P
         res = max(res, abs(flux))
         ! solve_factored_matrix_pivot
         !$acc loop seq
         do m=1, n-1
            pv = piv(m,b)
            if (pv /= m) then
               tmp = line(lid+nl*(m-1)) ; line(lid+nl*(m-1)) = line(lid+nl*(pv-1)) ; line(lid+nl*(pv-1)) = tmp
            endif
         enddo
         !$acc loop seq
         do m=2, n
            tmp = line(lid+nl*(m-1))
            !$acc loop seq
            do c=1, m-1
               tmp = tmp - mat(m,c,b) * line(lid+nl*(c-1))
            enddo
            line(lid+nl*(m-1)) = tmp
         enddo
         !$acc loop seq
         do m=n, 1, -1
            tmp = line(lid+nl*(m-1))
            !$acc loop seq
            do c=m+1, n
               tmp = tmp - mat(m,c,b) * line(lid+nl*(c-1))
            enddo
            line(lid+nl*(m-1)) = tmp / mat(m,m,b)
         enddo
         ! residual
         flux = 0._R8P
         !$acc loop seq
         do face=0, n
            if (face > 0) then
               call line_cell(dir, face, t1, t2, i, j, k)
               flux = flux - hdx*hq(b,i,j,k,dir)/dt
            endif
            ! odd mirror inlined (CPU twin: odd_mirror_index): a routine-seq version returned wrong values on device
            recon = 0._R8P
            !$acc loop seq
            do mm=1, hs
               coeff = fv1(mm)
               col = face + mm ; sgn = 1._R8P
               if (col < 1_I4P) then
                  col = 1_I4P - col ; sgn = -1._R8P
               elseif (col > n) then
                  col = 2_I4P*n + 1_I4P - col ; sgn = -1._R8P
               endif
               recon = recon + coeff*sgn*line(lid+nl*(col-1))
               col = face + 1_I4P - mm ; sgn = 1._R8P
               if (col < 1_I4P) then
                  col = 1_I4P - col ; sgn = -1._R8P
               elseif (col > n) then
                  col = 2_I4P*n + 1_I4P - col ; sgn = -1._R8P
               endif
               recon = recon + coeff*sgn*line(lid+nl*(col-1))
            enddo
            res = max(res, abs(recon-flux))
         enddo
         ! scatter
         !$acc loop seq
         do m=1, n
            call line_cell(dir, m, t1, t2, i, j, k)
            q(b,i,j,k,var) = line(lid+nl*(m-1))
         enddo
      enddo
      enddo
      enddo
      self%residual_max = max(self%residual_max, res)
   enddo
   call self%impose_virtual_ghosts(q=q, var_jx=var_jx)
   endsubroutine solve_lines

   subroutine solve_modified(self, field_fnl, q, var_jx, hq, dt)
   !< J_d = Q_d P_transverse F_d^E. CPU twins: solve_pic_charge_conserving_current_modified, filter_modified_axis,
   !< modified_current_line, build_esirkepov_face_line.
   !<
   !< Parallel mapping: (1) the two transverse P filters are fully parallel separable convolutions (one thread per cell,
   !< zero extension outside the block interior, taps summed in the CPU order); (2) the face prefix sum F^E is a scan along
   !< each line, one thread per line, written line-major into the scratch; (3) the Q convolution F^E -> J is again one
   !< thread per cell, reading the line-major scratch with line-fastest (coalesced) indexing. Only steps (2) is sequential
   !< and it is O(n) per line.
   class(prism_fnl_pic_conserving_object), intent(inout) :: self      !< Conserving current object.
   type(field_fnl_object),                 intent(in)    :: field_fnl !< Device field helper.
   real(R8P),                              intent(inout) :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Target field.
   integer(I4P),                           intent(in)    :: var_jx    !< Jx slot in q.
   real(R8P),                              intent(in)    :: hq(1:,1:,1:,1:,1:) !< Directional source [nb,ni,nj,nk,3].
   real(R8P),                              intent(in)    :: dt        !< Time step.
   real(R8P), pointer                                    :: line(:), qw(:), wb(:,:,:,:)
   real(R8P)                                             :: res, flux, hdx, sol
   integer(I4P)                                          :: dir, axis_a, axis_b, n, n1, n2, nl, nbl, b, t1, t2, lid
   integer(I4P)                                          :: m, t, f, i, j, k, var, qfirst, nq, last

   call zero_current(self=self, q=q, var_jx=var_jx)
   nbl = self%blocks_number ; qfirst = self%qfirst ; nq = self%nq
   line => self%line_gpu ; qw => self%qw_gpu ; wb => self%work_b_gpu
   do dir=1, 3
      hdx = self%vdx(dir)
      select case(dir)
      case(1) ; axis_a = 2 ; axis_b = 3 ; n = self%ni ; n1 = self%nj ; n2 = self%nk
      case(2) ; axis_a = 1 ; axis_b = 3 ; n = self%nj ; n1 = self%ni ; n2 = self%nk
      case(3) ; axis_a = 1 ; axis_b = 2 ; n = self%nk ; n1 = self%ni ; n2 = self%nj
      endselect
      call filter_axis(self=self, input=hq(:,:,:,:,dir), output=self%work_a_gpu, axis=axis_a)
      call filter_axis(self=self, input=self%work_a_gpu, output=self%work_b_gpu, axis=axis_b)
      nl  = nbl*n1*n2
      var = var_jx + dir - 1_I4P
      ! face prefix sums, one thread per line. Exact closure: the line sum of the source vanishes analytically (zero-flux
      ! closure, supports kept off the boundary), so the faces past the last nonzero source cell carry zero flux; their
      ! roundoff residual (still reported in residual_max) would otherwise run to the domain boundary as a constant current.
      ! CPU twin: modified_current_line.
      res = 0._R8P
      !$acc parallel loop independent gang vector collapse(3) DEVICEVAR(wb,line) &
      !$acc& firstprivate(dir,n,n1,nl,dt,hdx) private(lid,flux,m,i,j,k,last) reduction(max:res)
      !$omp OMPLOOP collapse(3) DEVICEPTR(wb,line) &
      !$omp& firstprivate(dir,n,n1,nl,dt,hdx) private(lid,flux,m,i,j,k,last) reduction(max:res)
      do t2=1, n2
      do t1=1, n1
      do b=1, nbl
         lid = b + nbl*((t1-1) + n1*(t2-1))
         flux = 0._R8P
         line(lid) = flux
         last = 0
         !$acc loop seq
         do m=1, n
            call line_cell(dir, m, t1, t2, i, j, k)
            flux = flux - hdx*wb(b,i,j,k)/dt
            line(lid+nl*m) = flux
            if (wb(b,i,j,k) /= 0._R8P) last = m
         enddo
         res = max(res, abs(flux))
         if (last > 0) then
            !$acc loop seq
            do m=last, n
               line(lid+nl*m) = 0._R8P
            enddo
         endif
      enddo
      enddo
      enddo
      self%residual_max = max(self%residual_max, res)
      ! Q convolution, one thread per cell
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q,line,qw) &
      !$acc& firstprivate(dir,n,n1,nl,var,qfirst,nq) private(lid,sol,f,i,j,k)
      !$omp OMPLOOP collapse(4) DEVICEPTR(q,line,qw) &
      !$omp& firstprivate(dir,n,n1,nl,var,qfirst,nq) private(lid,sol,f,i,j,k)
      do m=1, n
      do t2=1, n2
      do t1=1, n1
      do b=1, nbl
         lid = b + nbl*((t1-1) + n1*(t2-1))
         sol = 0._R8P
         !$acc loop seq
         do t=0, nq-1
            f = m + qfirst + t
            if (f < 0_I4P .or. f > n) cycle
            sol = sol + qw(t+1)*line(lid+nl*f)
         enddo
         call line_cell(dir, m, t1, t2, i, j, k)
         q(b,i,j,k,var) = sol
      enddo
      enddo
      enddo
      enddo
   enddo
   call self%impose_virtual_ghosts(q=q, var_jx=var_jx)
   endsubroutine solve_modified

   subroutine solve_modified_with_cleanup(self, pic_fnl, field_fnl, q, var_jx, q_ref, q_hist, active_stage, alph, beta, dt, &
                                          mixed_source)
   !< Per-particle modified current with tail cleanup, then accumulation. CPU twins: solve_modified_current_with_cleanup,
   !< cleanup_modified_particle_current, report_particle_tail_error.
   !<
   !< The particle loop stays on the host, as on the CPU (each particle needs a whole-grid modified solve, so the work per
   !< particle is already a set of grid-parallel kernels). The support box of each particle is obtained with min/max
   !< reductions (no device-to-host array copy). The tail check reports the largest offending ratio instead of the first
   !< offending cell. Everything runs on the virtual block (support boxes in global cells, so a box may straddle real
   !< blocks); the accumulated current is scattered to the real blocks at the end.
   class(prism_fnl_pic_conserving_object), intent(inout)        :: self          !< Conserving current object.
   type(prism_fnl_pic_object),             intent(in)           :: pic_fnl       !< Device PIC helper.
   type(field_fnl_object),                 intent(in)           :: field_fnl     !< Device field helper.
   real(R8P),                              intent(inout)        :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Real target.
   integer(I4P),                           intent(in)           :: var_jx        !< Jx slot in q.
   real(R8P),                              intent(in)           :: q_ref(1:,1:)  !< Trajectory start [np,8].
   real(R8P),                              intent(in)           :: q_hist(1:,1:,1:) !< Stage targets [np,8,>=stages].
   integer(I4P),                           intent(in)           :: active_stage  !< Active stage (0 = t0 virtual step).
   real(R8P),                              intent(in)           :: alph(1:,1:)   !< SSP alpha coefficients (host).
   real(R8P),                              intent(in)           :: beta(1:)      !< SSP beta coefficients (host).
   real(R8P),                              intent(in)           :: dt            !< Time step.
   real(R8P),                              intent(in), optional :: mixed_source(1:,1:,1:,1:,1:) !< Mixed stage source.
   real(R8P)                                                    :: coeff(1:8), j_max(3), j_ref, bad
   integer(I4P)                                                 :: stage_count, particle, t, r, block, dir, np
   integer(I4P)                                                 :: slo(3), shi(3), lo(3,3), hi(3,3)
   integer(I4P)                                                 :: p_min, p_max, q_min, q_max, s
   logical                                                      :: single

   np = self%particle_number
   stage_count = max(1_I4P, active_stage)
   ! the mixed stage source is the whole source: it is the particle's own only with one particle in the whole run
   single = self%particle_number_global == 1_I4P .and. np == 1_I4P .and. present(mixed_source)
   call zero_current(self=self, q=self%jv_gpu, var_jx=1_I4P)
   ! P/Q support offsets (CPU: cleanup_modified_particle_current)
   p_min = huge(1_I4P) ; p_max = -huge(1_I4P)
   do s=-self%radius, self%radius
      if (self%p(s) == 0._R8P) cycle
      p_min = min(p_min,s) ; p_max = max(p_max,s)
   enddo
   q_min = huge(1_I4P) ; q_max = -huge(1_I4P)
   do t=0, self%nq-1
      if (self%qw(t) == 0._R8P) cycle
      q_min = min(q_min,self%qfirst+t) ; q_max = max(q_max,self%qfirst+t)
   enddo

   do particle=1, np
      if (.not. single) then
         do t=1, stage_count
            call self%esirkepov_displacement(field_fnl=field_fnl, q_ref=q_ref, q_to=q_hist(:,:,t), h=self%src_gpu, &
                                             p_first=particle, p_last=particle)
            if (active_stage == 0_I4P) then
               call copy_source(src=self%src_gpu, dst=self%src_stage_gpu(:,:,:,:,:,t))
            else
               coeff = 0._R8P
               if (t < self%nrk) then
                  coeff(1:t) = alph(t+1,1:t)
               else
                  coeff(1:t) = beta(1:t)
               endif
               if (abs(coeff(t)) <= tiny(1._R8P)) &
                  call mpih%error_stop(msg=': zero SSP coefficient in particle current cleanup')
               call combine_particle_source(t=t, coeff=coeff)
            endif
         enddo
      endif

      call particle_support(particle=particle, block=block, slo=slo, shi=shi)
      ! the box and the source must come from the same trajectory end points
      if (single) then
         call self%check_source_support(h=mixed_source, slo=slo, shi=shi, particle=particle)
      else
         call self%check_source_support(h=self%src_stage_gpu(:,:,:,:,:,stage_count), slo=slo, shi=shi, particle=particle)
      endif
      do dir=1, 3
         lo(:,dir) = slo - p_max
         hi(:,dir) = shi - p_min
         lo(dir,dir) = slo(dir) - q_max
         hi(dir,dir) = shi(dir) - 1_I4P - q_min
      enddo

      if (single) then
         call self%solve_modified(field_fnl=field_fnl, q=self%jp_gpu, var_jx=1_I4P, hq=mixed_source, dt=dt)
      else
         call self%solve_modified(field_fnl=field_fnl, q=self%jp_gpu, var_jx=1_I4P, &
                                  hq=self%src_stage_gpu(:,:,:,:,:,stage_count), dt=dt)
      endif
      do dir=1, 3
         j_max(dir) = box_max(dir=dir, block=block, lo=lo(:,dir), hi=hi(:,dir))
      enddo
      j_ref = maxval(j_max)
      do dir=1, 3
         bad = clean_tail(dir=dir, block=block, lo=lo(:,dir), hi=hi(:,dir), j_ref=j_ref)
         if (bad > 0._R8P) then
            write(*,'(a)') 'ERROR: non-negligible current outside theoretical particle support'
            write(*,'(a,i0)') 'particle = ',particle
            write(*,'(a,i0)') 'component = ',dir
            write(*,'(a,i0)') 'block = ',block
            write(*,'(a,es24.16e3)') 'J_ref = ',j_ref
            write(*,'(a,es24.16e3)') 'max ratio = ',bad
            write(*,'(a,3(1x,i0))') 'support_lo = ',lo(:,dir)
            write(*,'(a,3(1x,i0))') 'support_hi = ',hi(:,dir)
            call mpih%error_stop(msg=': particle current outside theoretical support exceeds tail threshold')
         endif
      enddo
      call self%impose_virtual_ghosts(q=self%jp_gpu, var_jx=1_I4P)
      call add_particle_current
   enddo
   ! the particle currents of all the ranks
   call self%reduce_virtual(a=self%jv_gpu)
   call self%scatter_current(q=q, var_jx=var_jx)
   contains
      subroutine copy_source(src, dst)
      real(R8P),    intent(in)    :: src(1:,1:,1:,1:,1:)
      real(R8P),    intent(inout) :: dst(1:,1:,1:,1:,1:)
      integer(I4P)                :: b, i, j, k, d, ni, nj, nk, nbl
      ni = self%ni ; nj = self%nj ; nk = self%nk ; nbl = self%blocks_number
      !$acc parallel loop independent gang vector collapse(5) DEVICEVAR(src,dst)
      !$omp OMPLOOP collapse(5) DEVICEPTR(src,dst)
      do d=1, 3
      do k=1, nk
      do j=1, nj
      do i=1, ni
      do b=1, nbl
         dst(b,i,j,k,d) = src(b,i,j,k,d)
      enddo
      enddo
      enddo
      enddo
      enddo
      endsubroutine copy_source

      subroutine combine_particle_source(t, coeff)
      !< src_stage(t) = (src - sum_{r<t} coeff(r) src_stage(r)) / coeff(t).
      integer(I4P), intent(in) :: t
      real(R8P),    intent(in) :: coeff(1:8)
      real(R8P), pointer       :: src(:,:,:,:,:), stg(:,:,:,:,:,:)
      real(R8P)                :: val
      integer(I4P)             :: b, i, j, k, d, ni, nj, nk, nbl
      src => self%src_gpu ; stg => self%src_stage_gpu
      ni = self%ni ; nj = self%nj ; nk = self%nk ; nbl = self%blocks_number
      !$acc parallel loop independent gang vector collapse(5) DEVICEVAR(src,stg) copyin(coeff) private(val)
      !$omp OMPLOOP collapse(5) DEVICEPTR(src,stg) map(to:coeff) private(val)
      do d=1, 3
      do k=1, nk
      do j=1, nj
      do i=1, ni
      do b=1, nbl
         val = src(b,i,j,k,d)
         !$acc loop seq
         do r=1, t-1
            val = val - coeff(r)*stg(b,i,j,k,d,r)
         enddo
         stg(b,i,j,k,d,t) = val/coeff(t)
      enddo
      enddo
      enddo
      enddo
      enddo
      endsubroutine combine_particle_source

      subroutine particle_support(particle, block, slo, shi)
      !< Virtual-block cell box spanned by the particle shape over its reference and stage positions (min/max reductions).
      integer(I4P), intent(in)  :: particle
      integer(I4P), intent(out) :: block, slo(3), shi(3)
      real(R8P)                 :: e1, e2, e3, d1, d2, d3, s1, s2, s3
      integer(I4P)              :: tt, c1, c2, c3, r1, r2, r3, rx, ni, nj, nk, ierr
      integer(I4P)              :: lo1, lo2, lo3, hi1, hi2, hi3
      ni = self%ni ; nj = self%nj ; nk = self%nk
      e1 = self%vemin(1) ; e2 = self%vemin(2) ; e3 = self%vemin(3)
      ! base-shape cells, same arithmetic as esirkepov_displacement (dual nodes of the odd esirkepov_Q filters)
      s1 = self%shift(1) ; s2 = self%shift(2) ; s3 = self%shift(3)
      d1 = self%vdx(1)   ; d2 = self%vdx(2)   ; d3 = self%vdx(3)
      rx = support_radius(self)
      lo1 = huge(1_I4P) ; lo2 = huge(1_I4P) ; lo3 = huge(1_I4P)
      hi1 = -huge(1_I4P) ; hi2 = -huge(1_I4P) ; hi3 = -huge(1_I4P)
      ierr = 0_I4P
      ! Scalars only: with -fast, nvfortran 25.11 miscompiles a private x(3) overwritten twice in the iteration (first
      ! with q_ref, then with q_hist): the second store is lost and every stage target read as q_ref, so the box missed the
      ! stage targets (NGP/TSC tail errors at the first face crossing). Correct at -O0, on the host and without -acc.
      !$acc parallel loop independent gang vector DEVICEVAR(q_ref,q_hist) &
      !$acc& firstprivate(ni,nj,nk,rx,particle,e1,e2,e3,d1,d2,d3,s1,s2,s3) private(c1,c2,c3,r1,r2,r3) &
      !$acc& reduction(min:lo1,lo2,lo3) reduction(max:hi1,hi2,hi3,ierr)
      !$omp OMPLOOP DEVICEPTR(q_ref,q_hist) &
      !$omp& firstprivate(ni,nj,nk,rx,particle,e1,e2,e3,d1,d2,d3,s1,s2,s3) private(c1,c2,c3,r1,r2,r3) &
      !$omp& reduction(min:lo1,lo2,lo3) reduction(max:hi1,hi2,hi3,ierr)
      do tt=1, stage_count
         r1 = ceiling((q_ref(particle,1) - s1 - e1) / d1, kind=I4P)
         r2 = ceiling((q_ref(particle,2) - s2 - e2) / d2, kind=I4P)
         r3 = ceiling((q_ref(particle,3) - s3 - e3) / d3, kind=I4P)
         c1 = ceiling((q_hist(particle,1,tt) - s1 - e1) / d1, kind=I4P)
         c2 = ceiling((q_hist(particle,2,tt) - s2 - e2) / d2, kind=I4P)
         c3 = ceiling((q_hist(particle,3,tt) - s3 - e3) / d3, kind=I4P)
         if (r1 < 1_I4P .or. r1 > ni .or. r2 < 1_I4P .or. r2 > nj .or. r3 < 1_I4P .or. r3 > nk .or. &
             c1 < 1_I4P .or. c1 > ni .or. c2 < 1_I4P .or. c2 > nj .or. c3 < 1_I4P .or. c3 > nk) then
            ierr = max(ierr, ERR_OUTSIDE)
         else
            lo1 = min(lo1, r1-rx, c1-rx) ; hi1 = max(hi1, r1+rx, c1+rx)
            lo2 = min(lo2, r2-rx, c2-rx) ; hi2 = max(hi2, r2+rx, c2+rx)
            lo3 = min(lo3, r3-rx, c3-rx) ; hi3 = max(hi3, r3+rx, c3+rx)
         endif
      enddo
      if (ierr == ERR_OUTSIDE) call mpih%error_stop(msg=': PIC particle outside the domain in Esirkepov current cleanup')
      block = 1_I4P
      slo = [lo1, lo2, lo3]
      shi = [hi1, hi2, hi3]
      endsubroutine particle_support

      function box_max(dir, block, lo, hi) result(vmax)
      !< max |J_dir| over the expanded support box clipped to the block interior.
      integer(I4P), intent(in) :: dir, block, lo(3), hi(3)
      real(R8P)                :: vmax
      real(R8P), pointer       :: jp(:,:,:,:,:)
      integer(I4P)             :: i, j, k, i0, i1, j0, j1, k0, k1
      jp => self%jp_gpu
      i0 = max(1_I4P,lo(1)) ; i1 = min(self%ni,hi(1))
      j0 = max(1_I4P,lo(2)) ; j1 = min(self%nj,hi(2))
      k0 = max(1_I4P,lo(3)) ; k1 = min(self%nk,hi(3))
      vmax = 0._R8P
      !$acc parallel loop independent gang vector collapse(3) DEVICEVAR(jp) reduction(max:vmax)
      !$omp OMPLOOP collapse(3) DEVICEPTR(jp) reduction(max:vmax)
      do k=k0, k1
      do j=j0, j1
      do i=i0, i1
         vmax = max(vmax, abs(jp(block,i,j,k,dir)))
      enddo
      enddo
      enddo
      endfunction box_max

      function clean_tail(dir, block, lo, hi, j_ref) result(bad)
      !< Zero sub-threshold leakage outside the box; return the largest offending ratio (0 if none).
      integer(I4P), intent(in) :: dir, block, lo(3), hi(3)
      real(R8P),    intent(in) :: j_ref
      real(R8P)                :: bad
      real(R8P), pointer       :: jp(:,:,:,:,:)
      real(R8P)                :: value
      integer(I4P)             :: b, i, j, k, ni, nj, nk, nbl, l1, l2, l3, h1, h2, h3
      logical                  :: inside
      jp => self%jp_gpu
      ni = self%ni ; nj = self%nj ; nk = self%nk ; nbl = self%blocks_number
      l1 = lo(1) ; l2 = lo(2) ; l3 = lo(3) ; h1 = hi(1) ; h2 = hi(2) ; h3 = hi(3)
      bad = 0._R8P
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(jp) firstprivate(l1,l2,l3,h1,h2,h3,j_ref,block,dir) &
      !$acc& private(value,inside) reduction(max:bad)
      !$omp OMPLOOP collapse(4) DEVICEPTR(jp) firstprivate(l1,l2,l3,h1,h2,h3,j_ref,block,dir) &
      !$omp& private(value,inside) reduction(max:bad)
      do k=1, nk
      do j=1, nj
      do i=1, ni
      do b=1, nbl
         inside = b == block .and. i >= l1 .and. i <= h1 .and. j >= l2 .and. j <= h2 .and. &
                  k >= l3 .and. k <= h3
         if (.not. inside) then
            value = jp(b,i,j,k,dir)
            if (j_ref == 0._R8P) then
               if (value /= 0._R8P) bad = max(bad, huge(1._R8P))
            elseif (abs(value) < TAU_TAIL*j_ref) then
               jp(b,i,j,k,dir) = 0._R8P
            else
               bad = max(bad, abs(value)/j_ref)
            endif
         endif
      enddo
      enddo
      enddo
      enddo
      endfunction clean_tail

      subroutine add_particle_current
      !< jv(J) += jp, full virtual extent (ghosts included, as the CPU accumulation).
      real(R8P), pointer :: jp(:,:,:,:,:), jv(:,:,:,:,:)
      integer(I4P)       :: b, i, j, k, d, ni, nj, nk, nbl, ngc
      jp => self%jp_gpu ; jv => self%jv_gpu
      ni = self%ni ; nj = self%nj ; nk = self%nk ; nbl = self%blocks_number ; ngc = self%ngc
      !$acc parallel loop independent gang vector collapse(5) DEVICEVAR(jv,jp)
      !$omp OMPLOOP collapse(5) DEVICEPTR(jv,jp)
      do d=1, 3
      do k=1-ngc, nk+ngc
      do j=1-ngc, nj+ngc
      do i=1-ngc, ni+ngc
      do b=1, nbl
         jv(b,i,j,k,d) = jv(b,i,j,k,d) + jp(b,i,j,k,d)
      enddo
      enddo
      enddo
      enddo
      enddo
      endsubroutine add_particle_current
   endsubroutine solve_modified_with_cleanup

   subroutine deposit_modified_charge(self, pic_fnl, field_fnl, field, grid, q, nv)
   !< Deposit the unfiltered particle shape, then apply P_c in all directions. CPU twin: deposit_pic_modified_charge.
   !< The neighbour list must be current (caller runs the grid index first, as on the CPU).
   !<
   !< Multi-block: the regular deposit lands in the real blocks (ghost spill folded back with reduce_ghost_local_gpu), the
   !< complete interiors are gathered into the virtual block, filtered there (P reaches 7 cells, beyond ngc), and the result
   !< is scattered back to every real cell: interior and interface ghosts get the filtered values, physical ghosts 0 (as the
   !< single-block CPU, which zeroes the rho ghosts).
   class(prism_fnl_pic_conserving_object), intent(inout) :: self      !< Conserving current object.
   type(prism_fnl_pic_object),             intent(inout) :: pic_fnl   !< Device PIC helper.
   type(field_fnl_object),                 intent(in)    :: field_fnl !< Device field helper.
   type(field_object),                     intent(in)    :: field     !< Host field.
   type(grid_object),                      intent(in)    :: grid      !< Grid.
   real(R8P),                              intent(inout) :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Real field.
   integer(I4P),                           intent(in)    :: nv        !< Charge density slot.
   integer(I4P), pointer                                 :: nl_gpu(:,:), off(:,:)
   real(R8P),    pointer                                 :: wa(:,:,:,:), wb(:,:,:,:)
   integer(I4P)                                          :: n, np, rx, margin, ierr, ni, nj, nk, nbl, ngc, b, i, j, k
   integer(I4P)                                          :: nx, ny, nz, gi, gj, gk, v
   real(R8P),    pointer                                 :: q_shift(:,:), q_pic(:,:)
   real(R8P)                                             :: s1, s2, s3
   logical                                               :: old_filter

   ! periodic faces: built particle by particle and folded (deposit_matched_charge_box); without them the deposit below,
   ! whose results the box deposit reproduces to roundoff only, is kept bitwise
   if (self%box_current) then
      call self%deposit_matched_charge_box(pic_fnl=pic_fnl, q=q, nv=nv)
      return
   endif
   ni = self%bni ; nj = self%bnj ; nk = self%bnk ; nbl = self%bblocks ; ngc = self%ngc
   nx = self%ni ; ny = self%nj ; nz = self%nk
   rx = support_radius(self) ; margin = self%radius ; np = self%particle_number
   off => self%off_gpu
   if (np > 0_I4P) then
      nl_gpu => pic_fnl%neighbour_list_gpu
      ierr = 0_I4P
      !$acc parallel loop independent gang vector DEVICEVAR(nl_gpu,off) firstprivate(rx,margin,nx,ny,nz) &
      !$acc& private(gi,gj,gk) reduction(max:ierr)
      !$omp OMPLOOP DEVICEPTR(nl_gpu,off) firstprivate(rx,margin,nx,ny,nz) private(gi,gj,gk) reduction(max:ierr)
      do n=1, np
         if (nl_gpu(n,1) <= 0_I4P) then
            ierr = max(ierr, ERR_OUTSIDE)
         else
            ! global (virtual-block) cell: the filter margin is measured from the domain boundary
            gi = off(nl_gpu(n,1),1) + nl_gpu(n,2)
            gj = off(nl_gpu(n,1),2) + nl_gpu(n,3)
            gk = off(nl_gpu(n,1),3) + nl_gpu(n,4)
            if (gi-rx <= margin .or. nx-gi-rx <= margin .or. &
                gj-rx <= margin .or. ny-gj-rx <= margin .or. &
                gk-rx <= margin .or. nz-gk-rx <= margin) ierr = max(ierr, ERR_FILTER)
         endif
      enddo
      if (ierr == ERR_OUTSIDE) call mpih%error_stop(msg=': PIC particle outside the domain in '//self%solver_name//' charge')
      if (ierr == ERR_FILTER) call mpih%error_stop(msg=': '//self%solver_name//' P filter support reaches a domain boundary')
   endif
   old_filter = pic_fnl%filter_deposition
   pic_fnl%filter_deposition = .false.
   if (np > 0_I4P .and. any(self%shift /= 0._R8P)) then
      ! odd filters: base shape on the dual nodes, deposited as the particles at x-dx/2; then restore the real lookup
      q_shift => self%q_mixed_gpu ; q_pic => pic_fnl%q_pic_gpu
      s1 = self%shift(1) ; s2 = self%shift(2) ; s3 = self%shift(3)
      !$acc parallel loop independent gang vector collapse(2) DEVICEVAR(q_shift,q_pic) firstprivate(s1,s2,s3)
      !$omp OMPLOOP collapse(2) DEVICEPTR(q_shift,q_pic) firstprivate(s1,s2,s3)
      do v=1, 8
      do n=1, np
         select case(v)
         case(1)      ; q_shift(n,v) = q_pic(n,v) - s1
         case(2)      ; q_shift(n,v) = q_pic(n,v) - s2
         case(3)      ; q_shift(n,v) = q_pic(n,v) - s3
         case default ; q_shift(n,v) = q_pic(n,v)
         endselect
      enddo
      enddo
      call pic_fnl%particle_cartesian_grid_index_dev(field_fnl=field_fnl, field=field, grid=grid, q_pic_gpu=q_shift)
      call pic_fnl%particle_weighting_dev(field_fnl=field_fnl, field=field, grid=grid, q_gpu=q, q_pic_gpu=q_shift, nv=nv)
      call pic_fnl%particle_cartesian_grid_index_dev(field_fnl=field_fnl, field=field, grid=grid, q_pic_gpu=pic_fnl%q_pic_gpu)
   else
      call pic_fnl%particle_weighting_dev(field_fnl=field_fnl, field=field, grid=grid, q_gpu=q, q_pic_gpu=pic_fnl%q_pic_gpu, &
                                          nv=nv)
   endif
   pic_fnl%filter_deposition = old_filter
   wa => self%work_a_gpu ; wb => self%work_b_gpu
   ! gather the real interiors and the in-domain ghosts (a ghost deposit is the share of the virtual cell it mirrors: the
   ! ghost fold in global indices), then sum over the ranks; physical ghosts are skipped (zero after the margin check)
   !$acc parallel loop independent gang vector collapse(3) DEVICEVAR(wa)
   !$omp OMPLOOP collapse(3) DEVICEPTR(wa)
   do k=1, nz
   do j=1, ny
   do i=1, nx
      wa(1,i,j,k) = 0._R8P
   enddo
   enddo
   enddo
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q,wa,off) firstprivate(nv,nx,ny,nz) private(gi,gj,gk)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q,wa,off) firstprivate(nv,nx,ny,nz) private(gi,gj,gk)
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do b=1, nbl
      gi = off(b,1) + i ; gj = off(b,2) + j ; gk = off(b,3) + k
      if (gi >= 1 .and. gi <= nx .and. gj >= 1 .and. gj <= ny .and. gk >= 1 .and. gk <= nz) then
         if (q(b,i,j,k,nv) /= 0._R8P) then
            !$acc atomic update
            !$omp atomic update
            wa(1,gi,gj,gk) = wa(1,gi,gj,gk) + q(b,i,j,k,nv)
         endif
      endif
   enddo
   enddo
   enddo
   enddo
   call self%reduce_virtual(a=self%work_a_gpu)
   call filter_axis(self=self, input=self%work_a_gpu, output=self%work_b_gpu, axis=1_I4P)
   call filter_axis(self=self, input=self%work_b_gpu, output=self%work_a_gpu, axis=2_I4P)
   call filter_axis(self=self, input=self%work_a_gpu, output=self%work_b_gpu, axis=3_I4P)
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q,wb,off) firstprivate(nv,nx,ny,nz) private(gi,gj,gk)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q,wb,off) firstprivate(nv,nx,ny,nz) private(gi,gj,gk)
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do b=1, nbl
      gi = off(b,1) + i ; gj = off(b,2) + j ; gk = off(b,3) + k
      if (gi >= 1 .and. gi <= nx .and. gj >= 1 .and. gj <= ny .and. gk >= 1 .and. gk <= nz) then
         q(b,i,j,k,nv) = wb(1,gi,gj,gk)
      else
         q(b,i,j,k,nv) = 0._R8P
      endif
   enddo
   enddo
   enddo
   enddo
   endsubroutine deposit_modified_charge

   subroutine impose_current_ghosts(self, q, var_jx)
   !< Odd mirror (zero normal face flux) current ghosts on the PHYSICAL faces of the real blocks; interface ghosts are left
   !< to the ghost exchange. CPU twin: impose_odd_current_ghosts (single block: every face is physical). Every ghost write
   !< reads a cell that is interior in the component's own direction, so all three components are imposed in one race-free
   !< kernel.
   class(prism_fnl_pic_conserving_object), intent(in)    :: self   !< Conserving current object.
   real(R8P),                              intent(inout) :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Real field.
   integer(I4P),                           intent(in)    :: var_jx !< Jx slot in q.
   integer(I4P), pointer                                 :: off(:,:)
   integer(I4P)                                          :: b, i, j, k, ni, nj, nk, ngc, nbl, vx, vy, vz, nx, ny, nz
   logical                                               :: xm, xp, ym, yp, zm, zp
   logical                                               :: px, py, pz

   ni = self%bni ; nj = self%bnj ; nk = self%bnk ; ngc = self%ngc ; nbl = self%bblocks
   nx = self%ni ; ny = self%nj ; nz = self%nk
   vx = var_jx ; vy = var_jx + 1_I4P ; vz = var_jx + 2_I4P
   px = self%is_periodic(1) ; py = self%is_periodic(2) ; pz = self%is_periodic(3)
   off => self%off_gpu
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q,off) firstprivate(vx,vy,vz,nx,ny,nz,px,py,pz) &
   !$acc& private(xm,xp,ym,yp,zm,zp)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q,off) firstprivate(vx,vy,vz,nx,ny,nz,px,py,pz) private(xm,xp,ym,yp,zm,zp)
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do b=1, nbl
      ! physical faces only: periodic faces keep the exchanged images
      xm = off(b,1) == 0_I4P .and. (.not.px) ; xp = off(b,1) + ni == nx .and. (.not.px)
      ym = off(b,2) == 0_I4P .and. (.not.py) ; yp = off(b,2) + nj == ny .and. (.not.py)
      zm = off(b,3) == 0_I4P .and. (.not.pz) ; zp = off(b,3) + nk == nz .and. (.not.pz)
      if (i < 1  .and. xm) q(b,i,j,k,vx) = -q(b,1-i,j,k,vx)
      if (i > ni .and. xp) q(b,i,j,k,vx) = -q(b,2*ni+1-i,j,k,vx)
      if (j < 1  .and. ym) q(b,i,j,k,vy) = -q(b,i,1-j,k,vy)
      if (j > nj .and. yp) q(b,i,j,k,vy) = -q(b,i,2*nj+1-j,k,vy)
      if (k < 1  .and. zm) q(b,i,j,k,vz) = -q(b,i,j,1-k,vz)
      if (k > nk .and. zp) q(b,i,j,k,vz) = -q(b,i,j,2*nk+1-k,vz)
   enddo
   enddo
   enddo
   enddo
   endsubroutine impose_current_ghosts

   subroutine reduce_virtual4(self, a)
   !< Sum a rank-4 virtual-block device array over the ranks, in place: every rank adds the sources of its own particles,
   !< the sum is the source of all the particles. Staged through the host: an MPI_Allreduce on device buffers is not
   !< portable (the UCC/UCX collectives of HPC-X read device memory as host memory on WSL, SIGSEGV). No-op on one rank.
   !< Collective.
   class(prism_fnl_pic_conserving_object), intent(in)    :: self       !< Conserving current object.
   real(R8P),                              intent(inout) :: a(:,:,:,:) !< Device array.
   real(R8P), allocatable                                :: buf(:,:,:,:)

   if (mpih%procs_number == 1_I4P) return
   allocate(buf(size(a,1),size(a,2),size(a,3),size(a,4)))
   call dev_memcpy_from_device(dst=buf, src=a)
   call MPI_ALLREDUCE(MPI_IN_PLACE, buf, size(buf), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   call dev_memcpy_to_device(dst=a, src=buf)
   endsubroutine reduce_virtual4

   subroutine reduce_virtual5(self, a)
   !< Rank-5 twin of reduce_virtual4.
   class(prism_fnl_pic_conserving_object), intent(in)    :: self         !< Conserving current object.
   real(R8P),                              intent(inout) :: a(:,:,:,:,:) !< Device array.
   real(R8P), allocatable                                :: buf(:,:,:,:,:)

   if (mpih%procs_number == 1_I4P) return
   allocate(buf(size(a,1),size(a,2),size(a,3),size(a,4),size(a,5)))
   call dev_memcpy_from_device(dst=buf, src=a)
   call MPI_ALLREDUCE(MPI_IN_PLACE, buf, size(buf), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   call dev_memcpy_to_device(dst=a, src=buf)
   endsubroutine reduce_virtual5

   ! private methods
   subroutine build_virtual_geometry(self, field, grid)
   !< Build the virtual block of the whole domain: sizes, origin, spacing, cell centers and the origin of every local real
   !< block. Requirements (error otherwise): one refinement level (all blocks of all ranks with the same spacing), real
   !< blocks of all ranks tiling the domain exactly once. With one block the virtual block is the block itself. The virtual
   !< cell centers are the real ones of the blocks that hold them (gathered over the ranks): the Esirkepov weights use the
   !< field geometry of the run, identical on every rank.
   class(prism_fnl_pic_conserving_object), intent(inout) :: self  !< Conserving current object.
   type(field_object),                     intent(in)    :: field !< Host field.
   type(grid_object),                      intent(in)    :: grid  !< Grid.
   integer(I4P), allocatable                             :: off(:,:), slot(:,:,:)
   real(R8P),    allocatable                             :: vx(:), vy(:), vz(:)
   integer(I4P)                                          :: b, d, i, nbl, nvirt(3), nblk(3), nbs(3)
   real(R8P)                                             :: dx_min(3), dx_max(3)
   real(R8P), parameter                                  :: UNSET = -huge(1._R8P)

   nbl = field%blocks_number
   nblk = [self%bni, self%bnj, self%bnk]
   self%vemin = grid%domain_emin
   dx_min = huge(1._R8P) ; dx_max = -huge(1._R8P)
   do b=1, nbl
      dx_min = min(dx_min, field%dxyz(:,b)) ; dx_max = max(dx_max, field%dxyz(:,b))
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, dx_min, 3, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, dx_max, 3, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, mpih%error)
   if (any(dx_min /= dx_max)) &
      call mpih%error_stop(msg=': the FNL conserving PIC current needs a single refinement level (uniform block spacing)')
   self%vdx = dx_min
   do d=1, 3
      nvirt(d) = nint((grid%domain_emax(d) - grid%domain_emin(d)) / self%vdx(d), kind=I4P)
      if (mod(nvirt(d), nblk(d)) /= 0_I4P) &
         call mpih%error_stop(msg=': the FNL conserving PIC current: blocks do not tile the domain')
      nbs(d) = nvirt(d) / nblk(d)
   enddo
   allocate(off(1:max(1_I4P,self%bnb),1:3), slot(0:nbs(1)-1,0:nbs(2)-1,0:nbs(3)-1))
   off = 0_I4P ; slot = 0_I4P
   do b=1, nbl
      do d=1, 3
         off(b,d) = nint((field%emin(d,b) - self%vemin(d)) / self%vdx(d), kind=I4P)
         if (off(b,d) < 0_I4P .or. off(b,d) + nblk(d) > nvirt(d) .or. mod(off(b,d), nblk(d)) /= 0_I4P) &
            call mpih%error_stop(msg=': the FNL conserving PIC current: block origin off the virtual block lattice')
      enddo
      slot(off(b,1)/nblk(1),off(b,2)/nblk(2),off(b,3)/nblk(3)) = slot(off(b,1)/nblk(1),off(b,2)/nblk(2),off(b,3)/nblk(3)) + 1
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, slot, size(slot), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   if (any(slot /= 1_I4P)) &
      call mpih%error_stop(msg=': the FNL conserving PIC current: blocks do not tile the domain exactly once')
   ! virtual cell centers: the real ones, so that the Esirkepov weights use the field geometry of this run (blocks sharing a
   ! line of cells hold the same centers: the max over the ranks picks them exactly)
   allocate(vx(1:nvirt(1)), vy(1:nvirt(2)), vz(1:nvirt(3)))
   vx = UNSET ; vy = UNSET ; vz = UNSET
   do b=1, nbl
      do i=1, nblk(1) ; vx(off(b,1)+i) = field%x_cell(i,b) ; enddo
      do i=1, nblk(2) ; vy(off(b,2)+i) = field%y_cell(i,b) ; enddo
      do i=1, nblk(3) ; vz(off(b,3)+i) = field%z_cell(i,b) ; enddo
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, vx, nvirt(1), MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, vy, nvirt(2), MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, vz, nvirt(3), MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, mpih%error)
   if (any(vx == UNSET) .or. any(vy == UNSET) .or. any(vz == UNSET)) &
      call mpih%error_stop(msg=': the FNL conserving PIC current: virtual cell centers not covered by the blocks')
   self%ni = nvirt(1) ; self%nj = nvirt(2) ; self%nk = nvirt(3)
   self%bblocks = nbl
   call dev_assign_to_device(src=off, dst=self%off_gpu)
   call dev_assign_to_device(src=vx,  dst=self%vx_gpu)
   call dev_assign_to_device(src=vy,  dst=self%vy_gpu)
   call dev_assign_to_device(src=vz,  dst=self%vz_gpu)
   endsubroutine build_virtual_geometry

   subroutine scatter_current(self, q, var_jx)
   !< Real q(J) = virtual J at the same global cell, for every real cell including ghosts: interface ghosts receive the
   !< neighbour values, physical ghosts the virtual (single-block) odd/zero ghosts. Other variables are untouched.
   class(prism_fnl_pic_conserving_object), intent(in)    :: self   !< Conserving current object.
   real(R8P),                              intent(inout) :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Real field.
   integer(I4P),                           intent(in)    :: var_jx !< Jx slot in q.
   real(R8P),    pointer                                 :: jv(:,:,:,:,:)
   integer(I4P), pointer                                 :: off(:,:)
   integer(I4P)                                          :: b, i, j, k, d, ni, nj, nk, ngc, nbl, vj

   ni = self%bni ; nj = self%bnj ; nk = self%bnk ; ngc = self%ngc ; nbl = self%bblocks ; vj = var_jx
   jv => self%jv_gpu ; off => self%off_gpu
   !$acc parallel loop independent gang vector collapse(5) DEVICEVAR(q,jv,off) firstprivate(vj)
   !$omp OMPLOOP collapse(5) DEVICEPTR(q,jv,off) firstprivate(vj)
   do d=1, 3
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do b=1, nbl
      q(b,i,j,k,vj+d-1) = jv(1,off(b,1)+i,off(b,2)+j,off(b,3)+k,d)
   enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine scatter_current

   subroutine impose_virtual_ghosts(self, q, var_jx)
   !< Odd mirror current ghosts of the virtual block (all its faces are physical). CPU twin: impose_odd_current_ghosts.
   class(prism_fnl_pic_conserving_object), intent(in)    :: self   !< Conserving current object.
   real(R8P),                              intent(inout) :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Virtual field.
   integer(I4P),                           intent(in)    :: var_jx !< Jx slot in q.
   integer(I4P)                                          :: b, i, j, k, ni, nj, nk, ngc, nbl, vx, vy, vz

   ni = self%ni ; nj = self%nj ; nk = self%nk ; ngc = self%ngc ; nbl = self%blocks_number
   vx = var_jx ; vy = var_jx + 1_I4P ; vz = var_jx + 2_I4P
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q) firstprivate(vx,vy,vz)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q) firstprivate(vx,vy,vz)
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do b=1, nbl
      if (i < 1)  q(b,i,j,k,vx) = -q(b,1-i,j,k,vx)
      if (i > ni) q(b,i,j,k,vx) = -q(b,2*ni+1-i,j,k,vx)
      if (j < 1)  q(b,i,j,k,vy) = -q(b,i,1-j,k,vy)
      if (j > nj) q(b,i,j,k,vy) = -q(b,i,2*nj+1-j,k,vy)
      if (k < 1)  q(b,i,j,k,vz) = -q(b,i,j,1-k,vz)
      if (k > nk) q(b,i,j,k,vz) = -q(b,i,j,2*nk+1-k,vz)
   enddo
   enddo
   enddo
   enddo
   endsubroutine impose_virtual_ghosts

   subroutine impose_virtual_ghosts_periodic(self, q, var_jx)
   !< Current ghosts of the virtual block with periodic axes: along a periodic axis every component takes the periodic
   !< image, along the others the normal component is the odd mirror (impose_virtual_ghosts). One kernel per axis, in
   !< order, so that the edge and corner ghosts take the images of already filled ghosts. CPU twin: impose_odd_current_ghosts.
   class(prism_fnl_pic_conserving_object), intent(in)    :: self   !< Conserving current object.
   real(R8P),                              intent(inout) :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Virtual field.
   integer(I4P),                           intent(in)    :: var_jx !< Jx slot in q.
   integer(I4P)                                          :: b, i, j, k, g, v, ni, nj, nk, ngc, nbl, vj
   logical                                               :: per

   ni = self%ni ; nj = self%nj ; nk = self%nk ; ngc = self%ngc ; nbl = self%blocks_number ; vj = var_jx
   per = self%is_periodic(1)
   !$acc parallel loop independent gang vector collapse(5) DEVICEVAR(q) firstprivate(vj,per)
   !$omp OMPLOOP collapse(5) DEVICEPTR(q) firstprivate(vj,per)
   do v=0, 2
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do g=1, ngc
   do b=1, nbl
      if (per) then
         q(b,1-g,j,k,vj+v)  = q(b,ni+1-g,j,k,vj+v)
         q(b,ni+g,j,k,vj+v) = q(b,g,j,k,vj+v)
      elseif (v == 0) then
         q(b,1-g,j,k,vj)  = -q(b,g,j,k,vj)
         q(b,ni+g,j,k,vj) = -q(b,ni+1-g,j,k,vj)
      endif
   enddo
   enddo
   enddo
   enddo
   enddo
   per = self%is_periodic(2)
   !$acc parallel loop independent gang vector collapse(5) DEVICEVAR(q) firstprivate(vj,per)
   !$omp OMPLOOP collapse(5) DEVICEPTR(q) firstprivate(vj,per)
   do v=0, 2
   do k=1-ngc, nk+ngc
   do g=1, ngc
   do i=1-ngc, ni+ngc
   do b=1, nbl
      if (per) then
         q(b,i,1-g,k,vj+v)  = q(b,i,nj+1-g,k,vj+v)
         q(b,i,nj+g,k,vj+v) = q(b,i,g,k,vj+v)
      elseif (v == 1) then
         q(b,i,1-g,k,vj+1)  = -q(b,i,g,k,vj+1)
         q(b,i,nj+g,k,vj+1) = -q(b,i,nj+1-g,k,vj+1)
      endif
   enddo
   enddo
   enddo
   enddo
   enddo
   per = self%is_periodic(3)
   !$acc parallel loop independent gang vector collapse(5) DEVICEVAR(q) firstprivate(vj,per)
   !$omp OMPLOOP collapse(5) DEVICEPTR(q) firstprivate(vj,per)
   do v=0, 2
   do g=1, ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do b=1, nbl
      if (per) then
         q(b,i,j,1-g,vj+v)  = q(b,i,j,nk+1-g,vj+v)
         q(b,i,j,nk+g,vj+v) = q(b,i,j,g,vj+v)
      elseif (v == 2) then
         q(b,i,j,1-g,vj+2)  = -q(b,i,j,g,vj+2)
         q(b,i,j,nk+g,vj+2) = -q(b,i,j,nk+1-g,vj+2)
      endif
   enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine impose_virtual_ghosts_periodic

   subroutine check_source_support(self, h, slo, shi, particle)
   !< Stop if the directional source of one particle has nonzero cells outside the support box built from its trajectory
   !< end points. Box and source must describe the same trajectory; a mismatch would otherwise surface downstream as a
   !< misleading "current outside theoretical support" tail error. CPU twin: the same check in
   !< solve_modified_current_with_cleanup.
   class(prism_fnl_pic_conserving_object), intent(in) :: self              !< Conserving current object.
   real(R8P),                              intent(in) :: h(1:,1:,1:,1:,1:) !< Particle source [1,NX,NY,NZ,3].
   integer(I4P),                           intent(in) :: slo(3)            !< Support box, lower cell.
   integer(I4P),                           intent(in) :: shi(3)            !< Support box, upper cell.
   integer(I4P),                           intent(in) :: particle          !< Particle index (report only).
   integer(I4P)                                       :: b, i, j, k, d, ni, nj, nk, nbl
   integer(I4P)                                       :: e1, e2, e3, f1, f2, f3

   ni = self%ni ; nj = self%nj ; nk = self%nk ; nbl = self%blocks_number
   e1 = huge(1_I4P) ; e2 = huge(1_I4P) ; e3 = huge(1_I4P)
   f1 = -huge(1_I4P) ; f2 = -huge(1_I4P) ; f3 = -huge(1_I4P)
   !$acc parallel loop independent gang vector collapse(5) DEVICEVAR(h) reduction(min:e1,e2,e3) reduction(max:f1,f2,f3)
   !$omp OMPLOOP collapse(5) DEVICEPTR(h) reduction(min:e1,e2,e3) reduction(max:f1,f2,f3)
   do d=1, 3
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, nbl
      if (h(b,i,j,k,d) /= 0._R8P) then
         e1 = min(e1,i) ; e2 = min(e2,j) ; e3 = min(e3,k)
         f1 = max(f1,i) ; f2 = max(f2,j) ; f3 = max(f3,k)
      endif
   enddo
   enddo
   enddo
   enddo
   enddo
   if (e1 > f1) return ! empty source
   if (e1 < slo(1) .or. e2 < slo(2) .or. e3 < slo(3) .or. f1 > shi(1) .or. f2 > shi(2) .or. f3 > shi(3)) then
      write(*,'(a)') 'ERROR: particle source outside the support box of its trajectory end points'
      write(*,'(a,i0)') 'particle = ',particle
      write(*,'(a,3(1x,i0))') 'source_lo = ',e1,e2,e3
      write(*,'(a,3(1x,i0))') 'source_hi = ',f1,f2,f3
      write(*,'(a,3(1x,i0))') 'box_lo    = ',slo
      write(*,'(a,3(1x,i0))') 'box_hi    = ',shi
      call mpih%error_stop(msg=': particle source outside its trajectory support box in Esirkepov current cleanup')
   endif
   endsubroutine check_source_support

   ! private kernels
   subroutine zero_h(self, h)
   !< h = 0 over the interior directional source.
   class(prism_fnl_pic_conserving_object), intent(in)    :: self              !< Conserving current object.
   real(R8P),                              intent(inout) :: h(1:,1:,1:,1:,1:) !< Directional source [nb,ni,nj,nk,3].
   integer(I4P)                                          :: b, i, j, k, d, ni, nj, nk, nbl

   ni = self%ni ; nj = self%nj ; nk = self%nk ; nbl = self%blocks_number
   !$acc parallel loop independent gang vector collapse(5) DEVICEVAR(h)
   !$omp OMPLOOP collapse(5) DEVICEPTR(h)
   do d=1, 3
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, nbl
      h(b,i,j,k,d) = 0._R8P
   enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine zero_h

   subroutine zero_current(self, q, var_jx)
   !< q(J) = 0 over the full extent (CPU: q(var_Jx:var_Jz,:,:,:,:) = 0).
   class(prism_fnl_pic_conserving_object), intent(in)    :: self   !< Conserving current object.
   real(R8P),                              intent(inout) :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field.
   integer(I4P),                           intent(in)    :: var_jx !< Jx slot in q.
   integer(I4P)                                          :: b, i, j, k, d, ni, nj, nk, ngc, nbl, vj

   ni = self%ni ; nj = self%nj ; nk = self%nk ; ngc = self%ngc ; nbl = self%blocks_number ; vj = var_jx
   !$acc parallel loop independent gang vector collapse(5) DEVICEVAR(q) firstprivate(vj)
   !$omp OMPLOOP collapse(5) DEVICEPTR(q) firstprivate(vj)
   do d=0, 2
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do b=1, nbl
      q(b,i,j,k,vj+d) = 0._R8P
   enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine zero_current

   subroutine filter_axis(self, input, output, axis)
   !< Separable P convolution along one axis, zero extension outside the block interior. CPU twin: filter_modified_axis
   !< (taps summed from -radius to +radius, as on the CPU).
   class(prism_fnl_pic_conserving_object), intent(in)    :: self              !< Conserving current object.
   real(R8P),                              intent(in)    :: input(1:,1:,1:,1:)  !< Input [nb,ni,nj,nk].
   real(R8P),                              intent(inout) :: output(1:,1:,1:,1:) !< Output [nb,ni,nj,nk].
   integer(I4P),                           intent(in)    :: axis              !< Filter axis.
   real(R8P), pointer                                    :: p(:)
   real(R8P)                                             :: acc
   integer(I4P)                                          :: b, i, j, k, s, ii, jj, kk, ni, nj, nk, nbl, radius

   ni = self%ni ; nj = self%nj ; nk = self%nk ; nbl = self%blocks_number ; radius = self%radius
   p => self%p_gpu
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(input,output,p) firstprivate(axis,radius) &
   !$acc& private(acc,ii,jj,kk)
   !$omp OMPLOOP collapse(4) DEVICEPTR(input,output,p) firstprivate(axis,radius) private(acc,ii,jj,kk)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, nbl
      acc = 0._R8P
      !$acc loop seq
      do s=-radius, radius
         ii = i ; jj = j ; kk = k
         if (axis == 1) ii = i + s
         if (axis == 2) jj = j + s
         if (axis == 3) kk = k + s
         if (ii < 1 .or. ii > ni .or. jj < 1 .or. jj > nj .or. kk < 1 .or. kk > nk) cycle
         acc = acc + p(s+8)*input(b,ii,jj,kk)
      enddo
      output(b,i,j,k) = acc
   enddo
   enddo
   enddo
   enddo
   endsubroutine filter_axis

   subroutine transposed_filter_axis(self, input, output, axis)
   !< Transpose of filter_axis with the coherent gather kernel: output_m = sum_s p(s) input_{m-s}, zero extension. CPU twin:
   !< transposed_filter_axis of prism_pic_object.
   class(prism_fnl_pic_conserving_object), intent(in)    :: self              !< Conserving current object.
   real(R8P),                              intent(in)    :: input(1:,1:,1:,1:)  !< Input [nb,ni,nj,nk].
   real(R8P),                              intent(inout) :: output(1:,1:,1:,1:) !< Output [nb,ni,nj,nk].
   integer(I4P),                           intent(in)    :: axis              !< Filter axis.
   real(R8P), pointer                                    :: p(:)
   real(R8P)                                             :: acc
   integer(I4P)                                          :: b, i, j, k, s, ii, jj, kk, ni, nj, nk, nbl, radius

   ni = self%ni ; nj = self%nj ; nk = self%nk ; nbl = self%blocks_number ; radius = self%gradius
   p => self%gp_gpu
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(input,output,p) firstprivate(axis,radius) &
   !$acc& private(acc,ii,jj,kk)
   !$omp OMPLOOP collapse(4) DEVICEPTR(input,output,p) firstprivate(axis,radius) private(acc,ii,jj,kk)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, nbl
      acc = 0._R8P
      !$acc loop seq
      do s=-radius, radius
         ii = i ; jj = j ; kk = k
         if (axis == 1) ii = i - s
         if (axis == 2) jj = j - s
         if (axis == 3) kk = k - s
         if (ii < 1 .or. ii > ni .or. jj < 1 .or. jj > nj .or. kk < 1 .or. kk > nk) cycle
         acc = acc + p(s+8)*input(b,ii,jj,kk)
      enddo
      output(b,i,j,k) = acc
   enddo
   enddo
   enddo
   enddo
   endsubroutine transposed_filter_axis

   ! device helpers
   pure subroutine line_cell(dir, m, t1, t2, i, j, k)
   !< Cell (i,j,k) of the m-th element of the line (t1,t2) along dir.
   !$acc routine seq
   !$omp declare target
   integer(I4P), intent(in)  :: dir, m, t1, t2
   integer(I4P), intent(out) :: i, j, k

   select case(dir)
   case(1) ; i = m  ; j = t1 ; k = t2
   case(2) ; i = t1 ; j = m  ; k = t2
   case default ; i = t1 ; j = t2 ; k = m
   endselect
   endsubroutine line_cell

   pure subroutine shape_weights_dev(shp, filter, sigma, cutoff_limit, radius, x, ic, span, dx, xc, w, perr)
   !< 1D Esirkepov weights of one endpoint on the trajectory box. CPU twins: compute_esirkepov_gaussian_shape_1d,
   !< closed_shape_1d (unfiltered B-splines, closed weights of bspline_closed_weights), esirkepov_shape_weight_1d.
   !$acc routine seq
   !$omp declare target
   integer(I4P), intent(in)    :: shp          !< Shape.
   logical,      intent(in)    :: filter       !< Binomial filter on the Esirkepov shape.
   real(R8P),    intent(in)    :: sigma        !< Gaussian width.
   real(R8P),    intent(in)    :: cutoff_limit !< Gaussian cutoff.
   integer(I4P), intent(in)    :: radius       !< Support radius.
   real(R8P),    intent(in)    :: x            !< Particle coordinate.
   integer(I4P), intent(in)    :: ic           !< Particle cell, local box index.
   integer(I4P), intent(in)    :: span         !< Box cells.
   real(R8P),    intent(in)    :: dx           !< Cell size.
   real(R8P),    intent(in)    :: xc(MAX_SPAN) !< Cell centers of the box.
   real(R8P),    intent(out)   :: w(MAX_SPAN)  !< Weights on the box.
   integer(I4P), intent(inout) :: perr         !< Error code.
   real(R8P)                   :: weight_sum, shift, pu_error, xs(7), ws(7)
   integer(I4P)                :: m, smin, smax, ns

   do m=1, MAX_SPAN
      w(m) = 0._R8P
   enddo
   if (shp == SHAPE_GAUSSIAN) then
      shift = dx / sigma
      do m=ic-radius, ic+radius
         w(m) = effective_gaussian_weight_dev(r=(x-xc(m))/sigma, shift=shift, cutoff_limit=cutoff_limit, filter=filter)
      enddo
      weight_sum = 0._R8P
      do m=1, span
         weight_sum = weight_sum + w(m)
      enddo
      if (weight_sum <= tiny(1._R8P)) then
         perr = max(perr, ERR_GAUSS)
      else
         do m=1, span
            w(m) = w(m) / weight_sum
         enddo
      endif
   elseif (.not. filter) then
      ! Unfiltered B-spline (SHAPE_* = order): closed weights on the deposition stencil, the same ones as the rho deposit.
      ! The stencil lies inside the box (box radius >= stencil half-width).
      call set_bspline_stencil_dev(order=shp, x_p=x, x_c=xc(ic), i_p=ic, i_min=smin, i_max=smax)
      ns = smax - smin + 1_I4P
      do m=1, ns
         xs(m) = xc(smin+m-1)
      enddo
      call bspline_closed_weights_dev(order=shp, x_p=x, x_cell=xs, n=ns, dx=dx, w=ws, pu_error=pu_error)
      do m=1, ns
         w(smin+m-1) = ws(m)
      enddo
      ! NGP is discontinuous: a raw weight of 0 on the owning cell (ceiling/x_cell tie) is closed to 1, not an error
      if (shp > SHAPE_NGP .and. pu_error > TAU_PU) perr = max(perr, ERR_PU)
   else
      do m=1, span
         w(m) = esirkepov_shape_weight_dev(shp=shp, filter=filter, x=x, xc=xc(m), dx=dx)
      enddo
   endif
   endsubroutine shape_weights_dev

   pure function esirkepov_shape_weight_dev(shp, filter, x, xc, dx) result(w)
   !< Optionally binomial-filtered 1D particle shape. CPU twin: esirkepov_shape_weight_1d.
   !$acc routine seq
   !$omp declare target
   integer(I4P), intent(in) :: shp
   logical,      intent(in) :: filter
   real(R8P),    intent(in) :: x, xc, dx
   real(R8P)                :: w
   real(R8P), parameter     :: tol = 64._R8P*epsilon(1._R8P)
   real(R8P)                :: r, wm, wc, wp

   if (shp == SHAPE_NGP) then
      r = (x-xc)/dx
      w = 0._R8P
      if (filter) then
         if (r > -0.5_R8P+tol .and. r <= 0.5_R8P+tol) then
            w = 0.5_R8P
         elseif ((r > 0.5_R8P-tol .and. r <= 1.5_R8P+tol) .or. (r > -1.5_R8P-tol .and. r <= -0.5_R8P+tol)) then
            w = 0.25_R8P
         endif
      else
         if (r > -0.5_R8P .and. r <= 0.5_R8P) w = 1._R8P
      endif
      return
   endif
   wc = particle_shape_weight_dev(shp=shp, x=x, xc=xc, dx=dx)
   w = wc
   if (filter) then
      wm = particle_shape_weight_dev(shp=shp, x=x, xc=xc-dx, dx=dx)
      wp = particle_shape_weight_dev(shp=shp, x=x, xc=xc+dx, dx=dx)
      w = 0.25_R8P*wm + 0.5_R8P*wc + 0.25_R8P*wp
   endif
   endfunction esirkepov_shape_weight_dev

   pure function particle_shape_weight_dev(shp, x, xc, dx) result(w)
   !< 1D particle shape (polynomial forms of the CPU). CPU twin: particle_shape_weight_1d (Gaussian not reached here).
   !$acc routine seq
   !$omp declare target
   integer(I4P), intent(in) :: shp
   real(R8P),    intent(in) :: x, xc, dx
   real(R8P)                :: w
   real(R8P)                :: ar

   w = 0._R8P
   ar = abs((x - xc) / dx)
   select case(shp)
   case(SHAPE_NGP)
      if (ar <= 0.5_R8P) w = 1._R8P
   case(SHAPE_CIC)
      if (ar <= 1._R8P) w = 1._R8P - ar
   case(SHAPE_TSC)
      if (ar <= 0.5_R8P) then
         w = 0.75_R8P - ar**2
      elseif (ar <= 1.5_R8P) then
         w = 0.5_R8P * (1.5_R8P - ar)**2
      endif
   case(SHAPE_CUBIC)
      if (ar <= 1._R8P) then
         w = 2._R8P/3._R8P - ar**2 + 0.5_R8P*ar**3
      elseif (ar <= 2._R8P) then
         w = (2._R8P - ar)**3 / 6._R8P
      endif
   case(SHAPE_QUARTIC)
      if (ar <= 0.5_R8P) then
         w = 0.25_R8P*ar**4 - 5._R8P/8._R8P*ar**2 + 115._R8P/192._R8P
      elseif (ar <= 1.5_R8P) then
         w = -ar**4/6._R8P + 5._R8P*ar**3/6._R8P - 5._R8P*ar**2/4._R8P + 5._R8P*ar/24._R8P + 55._R8P/96._R8P
      elseif (ar <= 2.5_R8P) then
         w = (2.5_R8P - ar)**4 / 24._R8P
      endif
   case(SHAPE_QUINTIC)
      if (ar <= 1._R8P) then
         w = -ar**5/12._R8P + ar**4/4._R8P - ar**2/2._R8P + 11._R8P/20._R8P
      elseif (ar <= 2._R8P) then
         w = ar**5/24._R8P - 3._R8P*ar**4/8._R8P + 5._R8P*ar**3/4._R8P - 7._R8P*ar**2/4._R8P + &
             5._R8P*ar/8._R8P + 17._R8P/40._R8P
      elseif (ar <= 3._R8P) then
         w = (3._R8P - ar)**5 / 120._R8P
      endif
   case(SHAPE_SEXTIC)
      w = bspline_weight_dev(order=6_I4P, r=(x - xc)/dx)
   endselect
   endfunction particle_shape_weight_dev


   pure function effective_gaussian_weight_dev(r, shift, cutoff_limit, filter) result(weight)
   !< Optionally binomial-filtered 1D Gaussian weight. Same expressions as prism_pic_object%effective_gaussian_weight.
   !$acc routine seq
   !$omp declare target
   real(R8P), intent(in) :: r, shift, cutoff_limit
   logical,   intent(in) :: filter
   real(R8P)             :: weight

   if (filter) then
      weight = 0.25_R8P*gaussian_weight_dev(r=r + shift, cutoff_limit=cutoff_limit) &
             + 0.50_R8P*gaussian_weight_dev(r=r,         cutoff_limit=cutoff_limit) &
             + 0.25_R8P*gaussian_weight_dev(r=r - shift, cutoff_limit=cutoff_limit)
   else
      weight = gaussian_weight_dev(r=r, cutoff_limit=cutoff_limit)
   endif
   endfunction effective_gaussian_weight_dev

   pure function gaussian_weight_dev(r, cutoff_limit) result(weight)
   !< Unnormalized 1D Gaussian with compact cutoff.
   !$acc routine seq
   !$omp declare target
   real(R8P), intent(in) :: r, cutoff_limit
   real(R8P)             :: weight

   weight = 0._R8P
   if (abs(r) <= cutoff_limit) weight = exp(-0.5_R8P*r*r)
   endfunction gaussian_weight_dev

   ! host helpers
   function support_radius(self) result(rx)
   !< Integer support radius of the particle shape. CPU twin: particle_weighting_radius.
   class(prism_fnl_pic_conserving_object), intent(in) :: self !< Conserving current object.
   integer(I4P)                                       :: rx   !< Support radius.

   select case(self%shape)
   case(SHAPE_NGP)                  ; rx = 0_I4P
   case(SHAPE_CIC, SHAPE_TSC)       ; rx = 1_I4P
   case(SHAPE_CUBIC, SHAPE_QUARTIC) ; rx = 2_I4P
   case(SHAPE_QUINTIC, SHAPE_SEXTIC); rx = 3_I4P
   case(SHAPE_GAUSSIAN)             ; rx = self%gaussian_support_cells
   case default                     ; rx = 0_I4P
   endselect
   endfunction support_radius

   subroutine check_particle_error(ierr)
   !< Translate a device particle error code into the CPU error message.
   integer(I4P), intent(in) :: ierr !< Error code.

   select case(ierr)
   case(ERR_OUTSIDE)  ; call mpih%error_stop(msg=': PIC particle outside the domain in Esirkepov current')
   case(ERR_CHARGE)   ; call mpih%error_stop(msg=': esirkepov particle charge changed during an RK stage')
   case(ERR_BOUNDARY) ; call mpih%error_stop(msg=': esirkepov particle shape touches a domain boundary; charge flux BC is '// &
                                                 'unsupported')
   case(ERR_FILTER)   ; call mpih%error_stop(msg=': esirkepov_Q<N> P filter support reaches a domain boundary')
   case(ERR_SPAN)     ; call mpih%error_stop(msg=': esirkepov trajectory support exceeds the FNL limit MAX_SPAN='// &
                                                 trim(str(MAX_SPAN,.true.))//' cells')
   case(ERR_GAUSS)    ; call mpih%error_stop(msg=': empty Gaussian deposition support in Esirkepov current')
   case(ERR_PU)       ; call mpih%error_stop(msg=': esirkepov B-spline weights do not sum to one on the deposition stencil')
   endselect
   endsubroutine check_particle_error

   subroutine modified_kernels(esirkepov_q, p, qw, radius, qfirst, nq)
   !< P = A Q (base-shape nodes to centers) and Q (face prefix sums to centers) of esirkepov_Q<N>. CPU twin: modified_kernels
   !< (offsets and the half-cell shift of the odd filters documented there).
   integer(I4P), intent(in)  :: esirkepov_q
   real(R8P),    intent(out) :: p(-7:7), qw(0:9)
   integer(I4P), intent(out) :: radius, qfirst, nq

   p = 0._R8P ; qw = 0._R8P
   select case(esirkepov_q)
   case(5_I4P)
      p(-5:4) = real([1,2,3,49,1985,1985,49,3,2,1],R8P)/4080._R8P
      qw(0:4) = real([1,10,46,10,1],R8P)/68._R8P
      radius = 6_I4P ; qfirst = -3_I4P ; nq = 5_I4P
   case(6_I4P)
      p(-5:5) = real([1,1,1,118,1378,2522,1378,118,1,1,1],R8P)/5520._R8P
      qw(0:5) = real([1,9,36,36,9,1],R8P)/92._R8P
      radius = 5_I4P ; qfirst = -3_I4P ; nq = 6_I4P
   case(8_I4P)
      p(-6:6) = real([1,1,1,152,1167,3327,4622,3327,1167,152,1,1,1],R8P)/13920._R8P
      qw(0:7) = real([1,9,36,70,70,36,9,1],R8P)/232._R8P
      radius = 6_I4P ; qfirst = -4_I4P ; nq = 8_I4P
   case(10_I4P)
      p(-7:7) = real([1,1,1,166,1111,3487,6567,8052,6567,3487,1111,166,1,1,1],R8P)/30720._R8P
      qw(0:9) = real([1,9,36,84,126,126,84,36,9,1],R8P)/512._R8P
      radius = 7_I4P ; qfirst = -5_I4P ; nq = 10_I4P
   case default
      radius = 0_I4P ; qfirst = 0_I4P ; nq = 0_I4P
      call mpih%error_stop(msg=': unsupported esirkepov_Q filter '//trim(str(esirkepov_q, .true.)))
   endselect
   endsubroutine modified_kernels

   subroutine odd_mirror_index(idx, n, col, sgn)
   !< Odd mirror of a line index. CPU twin: odd_mirror_index.
   integer(I4P), intent(in)  :: idx, n
   integer(I4P), intent(out) :: col
   real(R8P),    intent(out) :: sgn

   if (idx < 1_I4P) then
      col = 1_I4P - idx
      sgn = -1._R8P
   elseif (idx > n) then
      col = 2_I4P*n + 1_I4P - idx
      sgn = -1._R8P
   else
      col = idx
      sgn = 1._R8P
   endif
   if (col < 1_I4P .or. col > n) call mpih%error_stop(msg=': centered derivative line solve needs more interior cells')
   endsubroutine odd_mirror_index

   subroutine build_centered_reconstruction_matrix(hs, n, a)
   !< Face-reconstruction matrix R such that FD1_CC J = (R_{i+1/2}-R_{i-1/2})/dx. CPU twin of the same name.
   integer(I4P), intent(in)    :: hs, n
   real(R8P),    intent(inout) :: a(0:,1:)
   integer(I4P)                :: face, m, col
   real(R8P)                   :: coeff, sgn

   a(0:n,1:n) = 0._R8P
   do face=0, n
      do m=1, hs
         coeff = FV1_CC(m,hs)
         call odd_mirror_index(idx=face+m, n=n, col=col, sgn=sgn)
         a(face,col) = a(face,col) + coeff * sgn
         call odd_mirror_index(idx=face+1-m, n=n, col=col, sgn=sgn)
         a(face,col) = a(face,col) + coeff * sgn
      enddo
   enddo
   endsubroutine build_centered_reconstruction_matrix

   subroutine factorize_matrix_pivot(a, pivot, n)
   !< LU with partial pivoting, in place. CPU twin of the same name.
   real(R8P),    intent(inout) :: a(1:,1:)
   integer(I4P), intent(out)   :: pivot(1:)
   integer(I4P), intent(in)    :: n
   integer(I4P)                :: i, j, k, p
   real(R8P)                   :: pivot_abs, factor, tmp

   do k=1, n-1
      p = k
      pivot_abs = abs(a(k,k))
      do i=k+1, n
         if (abs(a(i,k)) > pivot_abs) then
            p = i
            pivot_abs = abs(a(i,k))
         endif
      enddo
      if (pivot_abs <= 100._R8P*tiny(1._R8P)) call mpih%error_stop(msg=': singular PIC current derivative line solve')
      pivot(k) = p
      if (p /= k) then
         do j=1, n
            tmp = a(k,j) ; a(k,j) = a(p,j) ; a(p,j) = tmp
         enddo
      endif
      do i=k+1, n
         factor = a(i,k) / a(k,k)
         a(i,k) = factor
         do j=k+1, n
            a(i,j) = a(i,j) - factor * a(k,j)
         enddo
      enddo
   enddo
   if (abs(a(n,n)) <= 100._R8P*tiny(1._R8P)) call mpih%error_stop(msg=': singular PIC current derivative line solve')
   pivot(n) = n
   endsubroutine factorize_matrix_pivot
endmodule adam_prism_fnl_pic_conserving_object
