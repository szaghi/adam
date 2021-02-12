!< ADAM, Euler equations system class definition, GPU backend.
module adam_equation_euler_gpu_object
!< ADAM, Euler equations system class definition, GPU backend.

use adam_adam_object
use adam_base_gpu_object
use adam_field_object
use adam_grid_object
use adam_parameters
use adam_weno_library_gpu
use FiNeR
use PENF
use MPI
use CUDAFOR
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none
private
public :: equation_euler_gpu_object
public :: BC_EXTRAPOLATION
public :: BC_INFLOW

integer(I4P), parameter :: BC_EXTRAPOLATION = 1_I4P
integer(I4P), parameter :: BC_INFLOW        = 2_I4P

type :: equation_euler_gpu_object
   !< Euler equations system class definition, GPU backend.
   !<
   !< Multifluids is modeled by the standard thermodynamic model.
   !<
   !< The conservative varibales are arranged as follows:
   !<```
   !< q(1):    rho(1)
   !< q(2):    rho(2),
   !< ...
   !< q(ns):   rho(ns), specific density of last specie
   !< q(ns+1): rho * u
   !< q(ns+2): rho * v
   !< q(ns+3): rho * w
   !< q(ns+4): rho * E
   !<```
   !< Where `rho(s)` is the specific density of s-th specie, `rho=sum(rho(s))`, `[u,v,w]` is the velocity
   !< vector and `E` it the total specific internal energy. The auxiliary variables array is arranged as follows:
   !<```
   !< q_aux(1): c(1)
   !< q_aux(2): c(2)
   !< ...
   !< q_aux(ns): c(ns), specie concentration of last specie
   !< q_aux(ns+1): rho=sum(rho(s))
   !< q_aux(ns+2): u
   !< q_aux(ns+3): v
   !< q_aux(ns+4): w
   !< q_aux(ns+5): g
   !< q_aux(ns+6): p
   !<```
   !< Where `c` is the species concentration, i.e. rho/rho(s), `p` is the pressure and
   !< `g` is the specific heat ratio of the mixture, i.e.
   !<```
   !< cp = sum(rho(s)/rho * cp(s))
   !< cv = sum(rho(s)/rho * cv(s))
   !< g = cp / cv
   !<```
   type(adam_object),  pointer :: adam=>null()          !< ADAM.
   type(field_object), pointer :: field=>null()         !< The field.
   type(grid_object),  pointer :: grid=>null()          !< The grid.
   integer(I4P),       pointer :: ngc=>null()           !< Number of ghost cells.
   integer(I4P),       pointer :: ni=>null()            !< Number of cells in i direction.
   integer(I4P),       pointer :: nj=>null()            !< Number of cells in j direction.
   integer(I4P),       pointer :: nk=>null()            !< Number of cells in k direction.
   integer(I4P),       pointer :: nb=>null()            !< Total blocks number for MPI.
   integer(I4P),       pointer :: blocks_number=>null() !< Actual blocks number.
   integer(I4P),       pointer :: nv=>null()            !< Number of variables.
   type(base_gpu_object)       :: base_gpu              !< The base GPU handler.
   integer(I4P)                :: myrank=0_I4P          !< MPI rank process.
   integer(I4P)                :: procs_number=1_I4P    !< Number of MPI processes.
   integer(I4P)                :: error=0_I4P           !< Error traping flag.
   ! equation data
   integer(I4P)           :: ns=1_I4P               !< Number of fluid species.
   real(R8P), allocatable :: cp0(:)                 !< Specific heat at constant pressure of initial species.
   real(R8P), allocatable :: cv0(:)                 !< Specific heat at constant pressure of initial species.
   real(R8P)              :: dt=0._R8P              !< Maximum time step accordingly to CFL criterion.
   real(R8P)              :: CFL=0.3_R8P            !< CFL limit.
   logical                :: null_xyz(3)=[.false.,&
                                          .false.,&
                                          .false.]  !< Flag triggering 1D/2D simulations.
   real(R8P), allocatable :: q_aux(:,:,:,:,:)       !< Auxiliary cell centered variables.
   ! Runge-Kutta data
   integer(I4P)           :: nrk=3_I4P !< Runge-Kutta stages number.
   real(R8P), allocatable :: alph(:,:) !< RK alpha coefficients.
   real(R8P), allocatable :: beta(:)   !< RK beta coefficients.
   real(R8P), allocatable :: gamm(:)   !< RK gamma coefficients.
   ! cuf data
   integer(I4P), allocatable, device :: weno_stencils_gpu    !< WENO stencils number/dimension.
   real(R8P),    allocatable, device :: cp0_gpu(:)           !< Specific heat at constant pressure of initial species.
   real(R8P),    allocatable, device :: cv0_gpu(:)           !< Specific heat at constant pressure of initial species.
   real(R8P),    allocatable, device :: fp_gpu(:,:,:,:,:)    !< Positive fluxes.
   real(R8P),    allocatable, device :: fm_gpu(:,:,:,:,:)    !< Negative fluxes.
   real(R8P),    allocatable, device :: f_i_gpu(:,:,:,:,:)   !< Fluxes for i direction.
   real(R8P),    allocatable, device :: f_j_gpu(:,:,:,:,:)   !< Fluxes for i direction.
   real(R8P),    allocatable, device :: f_k_gpu(:,:,:,:,:)   !< Fluxes for i direction.
   real(R8P),    allocatable, device :: f_gpu(:,:,:,:,:)     !< Convective fluxes.
   real(R8P),    allocatable, device :: dxyz_gpu(:,:)        !< Space steps.
   real(R8P),    allocatable, device :: alph_gpu(:,:)        !< RK alpha coefficients.
   real(R8P),    allocatable, device :: beta_gpu(:)          !< RK beta coefficients.
   real(R8P),    allocatable, device :: gamm_gpu(:)          !< RK gamma coefficients.
   real(R8P),    allocatable, device :: q_aux_gpu(:,:,:,:,:) !< Auxiliary cell centered variables.
   real(R8P),    allocatable, device :: q_gpu(:,:,:,:,:)     !< Field cell centered variables stages.
   real(R8P),    allocatable, device :: q_s_gpu(:,:,:,:,:,:) !< RK Field cell centered variables stages.
   contains
      ! public methods
      procedure, pass(self) :: amr_update              !< Do AMR update.
      procedure, pass(self) :: compute_aux             !< Compute auxiliary variables.
      procedure, pass(self) :: compute_dt              !< Compute time step.
      procedure, pass(self) :: copy_cpu_gpu            !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu            !< Copy data from GPU to CPU.
      procedure, pass(self) :: destroy                 !< Destroy the equation.
      procedure, pass(self) :: initialize              !< Initialize the equation.
      procedure, pass(self) :: mark_by_grad_rho        !< Mark blocks to be refined/derefined by a `grad(rho)` value.
      procedure, pass(self) :: integrate               !< Runge Kutta integration of equation.
      procedure, pass(self) :: print_progress          !< Print simulation progress.
      procedure, pass(self) :: refine_uniform          !< Refine all blocks uniformly.
      procedure, pass(self) :: runge_kutta_initialize  !< Initialize Runge-Kutta data.
      procedure, pass(self) :: save_hdf5               !< Save simulation data in HDF5 format.
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions of equation.
      procedure, pass(self) :: set_initial_conditions  !< Set initial conditions of equation.
      procedure, pass(self) :: update_ghost_gpu        !< Update ghost cells and set boundary conditions.
      ! operators
      generic :: assignment(=) => eq_assign_eq      !< Overload `=`.
      procedure, pass(lhs), private :: eq_assign_eq !< Operator `=`.
endtype equation_euler_gpu_object

contains
   ! public methods
   subroutine amr_update(self, iterations)
   !< Do AMR update.
   class(equation_euler_gpu_object), intent(inout)        :: self            !< The equation.
   integer(I4P),                     intent(in), optional :: iterations      !< Number of AMR iterations.
   integer(I4P)                                           :: iterations_     !< Number of AMR iterations, local var.
   logical                                                :: is_grid_changed !< Flag to check grid changes.
   integer(I4P)                                           :: i               !< Counter.

   iterations_ = 1 ; if (present(iterations)) iterations_ = iterations
   amr: do i=1, iterations_
      call self%mark_by_grad_rho(grad_tol=0.05_R8P, delta_fine=0.006_R8P, delta_coarse=0.015_R8P)
      call self%update_ghost_gpu(q_gpu=self%q_gpu)
      call self%copy_gpu_cpu()
      call self%adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false., is_grid_changed=is_grid_changed)
      call self%copy_cpu_gpu
      if (.not.is_grid_changed) exit amr
   enddo amr
   endsubroutine amr_update

   subroutine compute_aux(self, q_gpu, q_aux_gpu)
   !< Compute auxiliary variables.
   class(equation_euler_gpu_object), intent(in)          :: self          !< The equation.
   real(R8P),                        intent(in),  device :: q_gpu(1:,                    &
                                                                  1-self%ngc:,&
                                                                  1-self%ngc:,&
                                                                  1-self%ngc:,&
                                                                  1:)     !< Conservative variables.
   real(R8P),                        intent(out), device :: q_aux_gpu(1:,                    &
                                                                      1-self%ngc:,&
                                                                      1-self%ngc:,&
                                                                      1-self%ngc:,&
                                                                      1:) !< Auxiliary variables.

   associate(blocks_number=>self%field%blocks_number, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, ns=>self%ns)
      call compute_aux_cuf(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, ns=ns, &
                           cp0_gpu=self%cp0_gpu, cv0_gpu=self%cv0_gpu, q_gpu=q_gpu, q_aux_gpu=q_aux_gpu)
   endassociate
   endsubroutine compute_aux

   subroutine compute_dt(self)
   !< Compute maximum time step accordingly to CFL stabilty criterion.
   class(equation_euler_gpu_object), intent(inout) :: self !< The equation.
   real(R8P)                                       :: umax !< Maximum speed of waves propagation.
   integer(I4P)                                    :: b    !< Counter.

   associate(blocks_number=>self%field%blocks_number, dxyz=>self%field%dxyz, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             ngc=>self%ngc, ns=>self%ns, q=>self%field%q, dt=>self%dt, CFL=>self%CFL)
      call self%compute_aux(q_gpu=self%q_gpu, q_aux_gpu=self%q_aux_gpu)
      dt = huge(1._R8P)
      do b=1, self%field%blocks_number
         call compute_umax_cuf(b, ni=ni, nj=nj, nk=nk, ngc=ngc, ns=ns,   &
                               dx=dxyz(1,b), dy=dxyz(2,b), dz=dxyz(3,b), &
                               q_aux_gpu=self%q_aux_gpu, umax=umax)
         dt = min(dt, minval(dxyz(:,b)) / umax * CFL)
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, self%error)
   endassociate
   endsubroutine compute_dt

   subroutine copy_cpu_gpu(self)
   !< Copy data from CPU to GPU.
   class(equation_euler_gpu_object), intent(inout) :: self        !< The base backend.
   real(R8P), allocatable                          :: dxyz_t(:,:) !< Space steps transposed.
   integer(I4P)                                    :: i, b        !< Counter.

   allocate(dxyz_t(1:self%field%nb,3))
   do b=1, self%field%blocks_number
      do i=1, 3
         dxyz_t(b,i) = self%field%dxyz(i,b)
      enddo
   enddo
   self%dxyz_gpu = dxyz_t
   call self%base_gpu%copy_transpose_cpu_gpu(q_cpu=self%field%q, q_gpu=self%q_gpu)
   call self%base_gpu%copy_cpu_gpu
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self, compute_q_aux)
   !< Copy data from GPU to CPU.
   class(equation_euler_gpu_object), intent(inout)        :: self          !< The base backend.
   logical,                          intent(in), optional :: compute_q_aux !< Flag to compute auxiliary variables.

   call self%base_gpu%copy_transpose_gpu_cpu(nv=self%nv, q_gpu=self%q_gpu, q_cpu=self%field%q)
   if (present(compute_q_aux)) then
      if (compute_q_aux) then
         call self%compute_aux(q_gpu=self%q_gpu, q_aux_gpu=self%q_aux_gpu)
         call self%base_gpu%copy_transpose_gpu_cpu(nv=self%ns+6, q_gpu=self%q_aux_gpu, q_cpu=self%q_aux)
      endif
   endif
   endsubroutine copy_gpu_cpu

   subroutine destroy(self)
   !< Destroy the equation.
   class(equation_euler_gpu_object), intent(inout) :: self  !< The equation.
   type(equation_euler_gpu_object)                 :: fresh !< Fresh equation.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, adam, ns, nrk, cp0, cv0, CFL, null_xyz, weno_stencils, fields_gpu_number)
   !< Initialize the equation.
   class(equation_euler_gpu_object), intent(inout)        :: self              !< The equation.
   type(adam_object),                intent(in), target   :: adam              !< ADAM.
   integer(I4P),                     intent(in), optional :: ns                !< Species number.
   integer(I4P),                     intent(in), optional :: nrk               !< Runge-Kutta stages number.
   real(R8P),                        intent(in), optional :: cp0(:)            !< Initial specific heats at constant pressure.
   real(R8P),                        intent(in), optional :: cv0(:)            !< Initial specific heats at constant volume.
   real(R8P),                        intent(in), optional :: CFL               !< CFL value.
   logical,                          intent(in), optional :: null_xyz(3)       !< Flag triggering 1D/2D simulations.
   integer(I4P),                     intent(in), optional :: weno_stencils     !< Number of WENO stencils.
   integer(I4P),                     intent(in), optional :: fields_gpu_number !< Number of fields allocated on GPU.
   integer(I4P)                                           :: v                 !< Counter.

   ! CPU data
   call self%destroy
   self%adam          => adam
   self%field         => adam%field
   self%grid          => adam%field%grid
   self%ngc           => adam%field%grid%ngc
   self%ni            => adam%field%grid%ni
   self%nj            => adam%field%grid%nj
   self%nk            => adam%field%grid%nk
   self%nb            => adam%field%nb
   self%blocks_number => adam%field%blocks_number
   self%nv            => adam%field%nv
   if (present(ns)) self%ns = ns
   if (self%nv - self%ns /= 4) then
      write(stderr, '(A)') 'ADAM-ERROR: field%nv must be euler%ns+4'
      call MPI_FINALIZE(self%error)
      stop
   endif
   call self%base_gpu%initialize(field=self%field, nv_aux=self%ns+6, fields_gpu_number=fields_gpu_number)
   if (present(cp0)) then
      self%cp0 = cp0
   else
      allocate(self%cp0(self%ns))
      self%cp0 = 1040._R8P
   endif
   if (present(cv0)) then
      self%cv0 = cv0
   else
      allocate(self%cv0(self%ns))
      self%cv0 = 742.85_R8P
   endif
   if (present(nrk)) self%nrk = nrk
   if (present(CFL)) self%CFL = CFL
   if (present(null_xyz)) self%null_xyz = null_xyz

   if (present(weno_stencils)) then
      self%weno_stencils_gpu = weno_stencils
   else
      self%weno_stencils_gpu = 1
   endif
   call weno_initialize(weno_stencils=self%weno_stencils_gpu)
   call self%runge_kutta_initialize
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   ! allocate large array
   associate(nv=>self%nv, ns=>self%ns, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, nb=>self%nb, nrk=>self%nrk)
   ! CPU data
   allocate(self%q_aux(1:ns+6,       &
                       1-ngc:ni+ngc, &
                       1-ngc:nj+ngc, &
                       1-ngc:nk+ngc, 1:nb))
   ! GPU data
   allocate(self%fp_gpu(1:nb,         &
                        1-ngc:ni+ngc, &
                        1-ngc:nj+ngc, &
                        1-ngc:nk+ngc, 1:nv))
   allocate(self%fm_gpu(1:nb,         &
                        1-ngc:ni+ngc, &
                        1-ngc:nj+ngc, &
                        1-ngc:nk+ngc, 1:nv))
   allocate(self%f_i_gpu(1-ngc:ni+ngc, &
                         1:nb,         &
                         1-ngc:nj+ngc, &
                         1-ngc:nk+ngc, 1:nv))
   allocate(self%f_j_gpu(1-ngc:nj+ngc, &
                         1:nb,         &
                         1-ngc:ni+ngc, &
                         1-ngc:nk+ngc, 1:nv))
   allocate(self%f_k_gpu(1-ngc:nk+ngc, &
                         1:nb,         &
                         1-ngc:ni+ngc, &
                         1-ngc:nj+ngc, 1:nv))
   allocate(self%f_gpu(1:nb, &
                       0:ni, &
                       0:nj, &
                       0:nk, 1:nv))
   allocate(self%dxyz_gpu(1:nb, 1:3))
   allocate(self%q_aux_gpu(1:nb,         &
                           1-ngc:ni+ngc, &
                           1-ngc:nj+ngc, &
                           1-ngc:nk+ngc, 1:ns+6))
   allocate(self%q_gpu(1:nb,         &
                       1-ngc:ni+ngc, &
                       1-ngc:nj+ngc, &
                       1-ngc:nk+ngc, 1:nv))
   allocate(self%q_s_gpu(1:nb,         &
                         1-ngc:ni+ngc, &
                         1-ngc:nj+ngc, &
                         1-ngc:nk+ngc, 1:nv, 1:nrk))
   endassociate
   ! copy data that is not variable during the simulation
   self%cp0_gpu  = self%cp0
   self%cv0_gpu  = self%cv0
   self%alph_gpu = self%alph
   self%beta_gpu = self%beta
   self%gamm_gpu = self%gamm
   endsubroutine initialize

   subroutine mark_by_grad_rho(self, grad_tol, delta_fine, delta_coarse, threshold)
   !< Mark blocks to be refined/derefined by a `grad(rho)` value.
   class(equation_euler_gpu_object), intent(inout)        :: self           !< The equation.
   real(R8P),                        intent(in)           :: grad_tol       !< Gradiend tolerance value.
   real(R8P),                        intent(in)           :: delta_fine     !< Maximum cell delta in fine grids.
   real(R8P),                        intent(in)           :: delta_coarse   !< Minimum cell delta in coarse grids.
   real(R8P),                        intent(in), optional :: threshold      !< Threshold for sphere proximity.
   real(R8P)                                              :: threshold_     !< Threshold for sphere proximity, local var.
   real(R8P)                                              :: max_cell_delta !< Maximum cell delta.
   real(R8P)                                              :: grad_rho       !< Value (max) of gradient of rho.
   integer(I4P)                                           :: b              !< Counter.

   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   self%field%refinements_needed = [(TO_NOT_TOUCH,b=1,self%blocks_number)]
   call self%update_ghost_gpu(q_gpu=self%q_gpu)
   associate (ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
              blocks_number=>self%blocks_number, ns=>self%ns, dxyz=>self%field%dxyz)
      call self%compute_aux(q_gpu=self%q_gpu, q_aux_gpu=self%q_aux_gpu)
      do b=1, blocks_number
         grad_rho = gradient_cuf(b=b, ni=ni, nj=nj, nk=nk, ngc=ngc, ns=ns, &
                                 dx=dxyz(1,b), dy=dxyz(2,b), dz=dxyz(3,b), q_gpu=self%q_aux_gpu)

         max_cell_delta = max_cell_delta_grad(grad=grad_rho)

         if (maxval(dxyz(:,b)) > max_cell_delta) then
            self%field%refinements_needed(b) = TO_BE_REFINED
         elseif (maxval(dxyz(:,b)) * threshold_ < max_cell_delta) then
            self%field%refinements_needed(b) = TO_BE_DEREFINED
         else
            self%field%refinements_needed(b) = TO_NOT_TOUCH
         endif
      enddo
   endassociate
   contains
      function gradient_cuf(b, ni, nj, nk, ngc, ns, dx, dy, dz, q_gpu) result(gradient)
      !< Gradient done by CUF threads.
      integer(I4P), intent(in)         :: b                                 !< Block index.
      integer(I4P), intent(in)         :: ni                                !< Grid cells number in I direction.
      integer(I4P), intent(in)         :: nj                                !< Grid cells number in J direction.
      integer(I4P), intent(in)         :: nk                                !< Grid cells number in K direction.
      integer(I4P), intent(in)         :: ngc                               !< Ghost cells number.
      integer(I4P), intent(in)         :: ns                                !< Species number.
      real(R8P),    intent(in)         :: dx                                !< X space step.
      real(R8P),    intent(in)         :: dy                                !< Y space step.
      real(R8P),    intent(in)         :: dz                                !< Z space step.
      real(R8P),    intent(in), device :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field component to which apply gradient.
      real(R8P)                        :: gradient                          !< Maximum gradient of q.
      real(R8P)                        :: grad                              !< Current gradient of q.
      integer(I4P)                     :: i, j, k                           !< Counter.
      integer(I4P)                     :: iercuda                           !< Error trapping flag for CUDAFortran.

      gradient = 0._R8P
      !$cuf kernel do(3) <<<*,*>>>
      do k=1, nk
         do j=1, nj
            do i=1, ni
               grad = sqrt(((q_gpu(b,i+1,j,k,ns+1) - q_gpu(b,i-1,j,k,ns+1))/(2*dx))**2 + &
                           ((q_gpu(b,i,j+1,k,ns+1) - q_gpu(b,i,j-1,k,ns+1))/(2*dy))**2 + &
                           ((q_gpu(b,i,j,k+1,ns+1) - q_gpu(b,i,j,k-1,ns+1))/(2*dz))**2)
               gradient = max(gradient, grad)

            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
      endfunction gradient_cuf

      function max_cell_delta_grad(grad) result(delta)
      !< Return the maximum cell delta given a gradient tollerance.
      real(R8P), intent(in) :: grad  !< Gradient value.
      real(R8P)             :: delta !< Maximum cell delta admissible.

      if (grad > grad_tol) then
         delta = delta_fine
      else
         delta = delta_coarse
      endif
      endfunction max_cell_delta_grad
   endsubroutine mark_by_grad_rho

   subroutine integrate(self, t, do_ghost_syncro, residual)
   !< Runge Kutta integration of field.
   class(equation_euler_gpu_object), intent(inout)         :: self             !< The equation.
   real(R8P),                        intent(in)            :: t                !< Time.
   logical,                          intent(in),  optional :: do_ghost_syncro  !< Flag to do syncrous ghost update.
   real(R8P),                        intent(out), optional :: residual         !< Global residual.
   logical                                                 :: do_ghost_syncro_ !< Flag to do syncrous ghost update, local var.
   integer(I4P)                                            :: s                !< Counter.

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   associate(alph=>self%alph, beta=>self%beta, gamm=>self%gamm, dt=>self%dt, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             ngc=>self%ngc, nv=>self%nv, nrk=>self%nrk, ns=>self%ns, blocks_number=>self%blocks_number,             &
             inner_blocks_number=>self%field%inner_blocks_number,                                                   &
             alph_gpu=>self%alph_gpu, beta_gpu=>self%beta_gpu)
   do s=1, nrk
      call compute_rk_stage_gpu_cuf(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                    alph_gpu=alph_gpu, dt=dt, s=s, q_gpu=self%q_gpu, q_s_gpu=self%q_s_gpu)
      if (do_ghost_syncro_) then
         call self%update_ghost_gpu(q_gpu=self%q_s_gpu(:,:,:,:,:,s)) ! all ghosts
         call self%compute_aux(q_gpu=self%q_s_gpu(:,:,:,:,:,s), q_aux_gpu=self%q_aux_gpu)
         call compute_residuals_gpu_cuf(ni=ni, nj=nj, nk=nk, ngc=ngc, ns=ns, blocks_number=blocks_number,          &
                                        null_x=self%null_xyz(1), null_y=self%null_xyz(2), null_z=self%null_xyz(3), &
                                        dx_gpu            = self%dxyz_gpu(:,1),                                    &
                                        dy_gpu            = self%dxyz_gpu(:,2),                                    &
                                        dz_gpu            = self%dxyz_gpu(:,3),                                    &
                                        q_aux_gpu         = self%q_aux_gpu,                                        &
                                        weno_stencils_gpu = self%weno_stencils_gpu,                                &
                                        fp_gpu            = self%fp_gpu,                                           &
                                        fm_gpu            = self%fm_gpu,                                           &
                                        f_i_gpu           = self%f_i_gpu,                                          &
                                        f_j_gpu           = self%f_j_gpu,                                          &
                                        f_k_gpu           = self%f_k_gpu,                                          &
                                        f_gpu             = self%f_gpu,                                            &
                                        q_gpu             = self%q_s_gpu(:,:,:,:,:,s))
      else
         ! TODO
      endif
      if (present(residual).and.s==self%nrk) then
         ! TODO
      endif
   enddo
   call advance_q_gpu_cuf(ni=ni, nj=nj, nk=nk, ngc=ngc, nrk=nrk, nv=nv, blocks_number=blocks_number, &
                           beta_gpu=beta_gpu, dt=dt, q_s_gpu=self%q_s_gpu, q_gpu=self%q_gpu)
   endassociate
   endsubroutine integrate

   subroutine print_progress(self, t, time, time_max)
   !< Print simulation progress.
   class(equation_euler_gpu_object), intent(in) :: self     !< The equation.
   integer(I4P),                     intent(in) :: t        !< Time iteration.
   real(R8P),                        intent(in) :: time     !< Time.
   real(R8P)                                    :: time_max !< Maximum time of integration.

   print '(A)', 'blocks number: '//trim(str(self%adam%tree%nodes_number, .true.))
   print '(A)', 'time step:     '//trim(str(self%dt, .true.))
   print '(A)', 'time:          '//trim(str(time, .true.))
   print '(A)', 't:             '//trim(str(t,.true.))
   print '(A)', 'progress:      '//trim(str(int(time/time_max * 100), .true.))//'%'
   endsubroutine print_progress

   subroutine refine_uniform(self, refinement_levels)
   !< Refine all blocks uniformly.
   class(equation_euler_gpu_object), intent(inout) :: self              !< The equation.
   integer(I4P),                     intent(in)    :: refinement_levels !< Number of refinement to be performed.
   integer(I4P)                                    :: l                 !< Counter.

   do l=1, refinement_levels
      call self%adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
      call self%adam%amr_update(do_blocks_reorder=.false., do_mpi_redistribute=.true.)
   enddo
   endsubroutine

   subroutine runge_kutta_initialize(self)
   !< Initialize Runge-Kutta data.
   class(equation_euler_gpu_object), intent(inout) :: self !< The equation.

   allocate(self%alph(self%nrk,self%nrk), self%beta(self%nrk), self%gamm(self%nrk))
   select case(self%nrk)
   case(1_I4P)
      self%alph(:,:) = reshape([1._R8P], [1,1])
      self%beta(:) = [1._R8P]
      self%gamm(:) = [0._R8P]
   case(3_I4P)
      self%alph(:,:) = reshape([0._R8P, 1._R8P, 0.25_R8P, &
                                0._R8P, 0._R8P, 0.25_R8P, &
                                0._R8P, 0._R8P,0._R8P], [3,3])
      self%beta(:) = [1._R8P/6._R8P, &
                      1._R8P/6._R8P, &
                      2._R8P/3._R8P]
      self%gamm(:) = [0._R8P, &
                      1._R8P, &
                      0._R8P]
   endselect
   endsubroutine runge_kutta_initialize

   subroutine save_hdf5(self, output_basename, t, time)
   !< Save simulation data in HDF5 format.
   class(equation_euler_gpu_object), intent(inout) :: self            !< The equation.
   character(*),                     intent(in)    :: output_basename !< Output base name.
   integer(I4P),                     intent(in)    :: t               !< Time iteration.
   real(R8P),                        intent(in)    :: time            !< Time.

   call self%copy_gpu_cpu(compute_q_aux=.true.)
   call self%adam%save_hdf5(basename=trim(output_basename)//trim(strz(t,9)),  &
                            q=self%field%q,                                   &
                            q_aux=self%q_aux,                                 &
                            q_name=['rho  ','rho-u','rho-v','rho-w','rho-E'], &
                            q_aux_name=['c1','r ','u ','v ','w ','g ','p '],  &
                            with_cell_morton=.true.)
   endsubroutine save_hdf5

   subroutine set_boundary_conditions(self, q_gpu)
   !< Set boundary conditions of equation.
   class(equation_euler_gpu_object), intent(in)            :: self                  !< The equation.
   real(R8P),                        intent(inout), device :: q_gpu(1:,         &
                                                                    1-self%ngc:,&
                                                                    1-self%ngc:,&
                                                                    1-self%ngc:,1:) !< Conservative variables.

   if (allocated(self%base_gpu%local_map_bc_crown_gpu)) call set_bc(nv=self%nv, ngc=self%ngc, &
                                                                    local_map_bc=self%base_gpu%local_map_bc_crown_gpu)
   contains
      subroutine set_bc(nv, ngc, local_map_bc)
      integer(I4P), intent(in)         :: nv                  !< Number of variables.
      integer(I4P), intent(in)         :: ngc                 !< Ghost cells number.
      integer(I8P), intent(in), device :: local_map_bc(:,:,:) !< Local map for BC ghost cells.
      integer(I4P)                     :: b                   !< Counter.
      integer(I4P)                     :: c, i, j, k, v       !< Counter.
      integer(I4P)                     :: idelta              !< IJK delta step for extrapolation.
      integer(I4P)                     :: jdelta              !< IJK delta step for extrapolation.
      integer(I4P)                     :: kdelta              !< IJK delta step for extrapolation.
      integer(I4P)                     :: bc_type             !< Boundary condition type.
      integer(I4P)                     :: crown               !< Crown counter.
      integer(I4P)                     :: iercuda             !< Error trapping flag for CUDAFortran.

      do crown=1, ngc
         !$cuf kernel do(1) <<<*,*>>>
         do c=1, size(local_map_bc, dim=1)
            b = local_map_bc(c, 1 ,crown)
            if (b>0) then
               i       = local_map_bc(c, 2 ,crown)
               j       = local_map_bc(c, 3 ,crown)
               k       = local_map_bc(c, 4 ,crown)
               idelta  = local_map_bc(c, 5 ,crown)
               jdelta  = local_map_bc(c, 6 ,crown)
               kdelta  = local_map_bc(c, 7 ,crown)
               bc_type = local_map_bc(c, 8 ,crown)
               if (bc_type == BC_EXTRAPOLATION) then
                  do v=1, nv
                     q_gpu(b,i,j,k,v) = q_gpu(b,i-idelta,j-jdelta,k-kdelta,v)
                  enddo
               else
               endif
            endif
         enddo
         !@cuf iercuda=cudaDeviceSynchronize()
      enddo
      endsubroutine set_bc
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self, nic_regions, emin_icr, emax_icr, rho_icr, velocity_icr, pressure_icr)
   !< Set initial conditions of field.
   class(equation_euler_gpu_object), intent(inout) :: self                         !< The equation.
   integer(I4P),                     intent(in)    :: nic_regions                  !< Number of initial conditions regions.
   real(R8P),                        intent(in)    :: emin_icr(:,:), emax_icr(:,:) !< Initial conditions regions extents.
   real(R8P),                        intent(in)    :: rho_icr(:)                   !< Initial conditions regions density.
   real(R8P),                        intent(in)    :: velocity_icr(:,:)            !< Initial conditions regions velocity.
   real(R8P),                        intent(in)    :: pressure_icr(:)              !< Initial conditions regions pressure.
   integer(I4P)                                    :: b, i, j, k, icr              !< Counter.

   associate(blocks_number=>self%blocks_number, q=>self%field%q, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             ngc=>self%ngc, x_cell=>self%field%x_cell, y_cell=>self%field%y_cell, z_cell=>self%field%z_cell)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               do icr=1, nic_regions
                  if (x_cell(i,b)>=emin_icr(1,icr).and.x_cell(i,b)<emax_icr(1,icr).and. &
                      y_cell(j,b)>=emin_icr(2,icr).and.y_cell(j,b)<emax_icr(2,icr).and. &
                      z_cell(k,b)>=emin_icr(3,icr).and.z_cell(k,b)<emax_icr(3,icr)) then
                     q(1,i,j,k,b) = rho_icr(icr)
                     q(2,i,j,k,b) = q(1,i,j,k,b) * velocity_icr(1,icr)
                     q(3,i,j,k,b) = q(1,i,j,k,b) * velocity_icr(2,icr)
                     q(4,i,j,k,b) = q(1,i,j,k,b) * velocity_icr(3,icr)
                     q(5,i,j,k,b) = q(1,i,j,k,b) * E(p=pressure_icr(icr),                                          &
                                                     r=q(1,i,j,k,b),                                               &
                                                     u=sqrt(dot_product(velocity_icr(:,icr),velocity_icr(:,icr))), &
                                                     g=self%cp0(1)/self%cv0(1))

                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   endassociate
   call self%copy_cpu_gpu
      !subroutine set_ic_riemann_2D(x, y, q_ic)
      !!< Set initial conditions for 1D Riemann problem.
      !real(R8P),    intent(in)  :: x        !< X (or Y or Z) coordinate.
      !real(R8P),    intent(in)  :: y        !< Y (or X or Z) coordinate.
      !real(R8P),    intent(out) :: q_ic(1:) !< Initial conditions to be set.
      !real(R8P)                 :: mod_u    !< Module of velocity vector.

      !select case(trim(ic_type_))
      !case('kt-02')
      !   if (x<=0.5_R8P.and.y<=0.5_R8P) then
      !      q_ic(1) = 0.138_R8P
      !      q_ic(2) = q_ic(1) * 1.206_R8P
      !      q_ic(3) = 0._R8P
      !      q_ic(4) = q_ic(1) * 1.206_R8P
      !      mod_u   = sqrt((q_ic(2)/q_ic(1))**2 + (q_ic(3)/q_ic(1))**2 + (q_ic(4)/q_ic(1))**2)
      !      q_ic(5) = q_ic(1) * E(p=0.029_R8P, r=q_ic(1), u=mod_u, g=self%cp0(1)/self%cv0(1))
      !   elseif (x>0.5_R8P.and.y<=0.5_R8P) then
      !      q_ic(1) = 0.5323_R8P
      !      q_ic(2) = 0._R8P
      !      q_ic(3) = 0._R8P
      !      q_ic(4) = q_ic(1) * 1.206_R8P
      !      mod_u   = sqrt((q_ic(2)/q_ic(1))**2 + (q_ic(3)/q_ic(1))**2 + (q_ic(4)/q_ic(1))**2)
      !      q_ic(5) = q_ic(1) * E(p=0.3_R8P, r=q_ic(1), u=mod_u, g=self%cp0(1)/self%cv0(1))
      !   elseif (x<=0.5_R8P.and.y>0.5_R8P) then
      !      q_ic(1) = 0.5323_R8P
      !      q_ic(2) = q_ic(1) * 1.206_R8P
      !      q_ic(3) = 0._R8P
      !      q_ic(4) = 0._R8P
      !      mod_u   = sqrt((q_ic(2)/q_ic(1))**2 + (q_ic(3)/q_ic(1))**2 + (q_ic(4)/q_ic(1))**2)
      !      q_ic(5) = q_ic(1) * E(p=0.3_R8P, r=q_ic(1), u=mod_u, g=self%cp0(1)/self%cv0(1))
      !   else
      !      q_ic(1) = 1.5_R8P
      !      q_ic(2) = 0._R8P
      !      q_ic(3) = 0._R8P
      !      q_ic(4) = 0._R8P
      !      mod_u   = sqrt((q_ic(2)/q_ic(1))**2 + (q_ic(3)/q_ic(1))**2 + (q_ic(4)/q_ic(1))**2)
      !      q_ic(5) = q_ic(1) * E(p=1.5_R8P, r=q_ic(1), u=mod_u, g=self%cp0(1)/self%cv0(1))
      !   endif
      !endselect
      !endsubroutine set_ic_riemann_2D
   endsubroutine set_initial_conditions

   subroutine update_ghost_gpu(self, q_gpu, step)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(equation_euler_gpu_object), intent(inout)         :: self            !< The equation.
   real(R8P),                        intent(inout), device :: q_gpu(1:,         &
                                                                    1-self%ngc:,&
                                                                    1-self%ngc:,&
                                                                    1-self%ngc:,&
                                                                    1:)       !< Conservative variables.
   integer(I4P),                     intent(in), optional  :: step            !< Step to be perfordmed in asyncronous comp.
   logical                                                 :: do_local_update !< Flag for triggering local update.
   logical                                                 :: do_set_bc       !< Flag for triggering setting bc.

   ! perform local update if step is not speficied or if first step is selected
   do_local_update = .false.
   do_set_bc       = .false.
   if (.not.present(step)) then
      do_local_update = .true.
      do_set_bc       = .true.
   else
      if (step==1) do_local_update = .true.
      if (step==3) do_set_bc       = .true.
   endif

   if (do_local_update) call self%base_gpu%update_ghost_local_gpu(q_gpu=q_gpu)
                        call self%base_gpu%update_ghost_mpi_gpu(q_gpu=q_gpu, step=step)
   if (do_set_bc)       call self%set_boundary_conditions(q_gpu=q_gpu)
   endsubroutine update_ghost_gpu

   ! operators
   ! =
   subroutine eq_assign_eq(lhs, rhs)
   !< Operator `=`.
   class(equation_euler_gpu_object), intent(inout) :: lhs !< Left hand side.
   type(equation_euler_gpu_object),  intent(in)    :: rhs !< Right hand side.

   lhs%adam          => rhs%adam
   lhs%field         => rhs%field
   lhs%grid          => rhs%grid
   lhs%ni            => rhs%ni
   lhs%nj            => rhs%nj
   lhs%nk            => rhs%nk
   lhs%ngc           => rhs%ngc
   lhs%nb            => rhs%nb
   lhs%blocks_number => rhs%blocks_number
   lhs%nv            => rhs%nv
   lhs%base_gpu = rhs%base_gpu
   lhs%myrank = rhs%myrank
   lhs%procs_number = rhs%procs_number
   lhs%error = rhs%error
   lhs%ns = rhs%ns
   lhs%dt = rhs%dt
   lhs%CFL = rhs%CFL
   lhs%null_xyz = rhs%null_xyz
   lhs%nrk = rhs%nrk
   call assign_allocatable(lhs=lhs%q_aux, rhs=rhs%q_aux )
   call assign_allocatable(lhs=lhs%cp0 ,  rhs=rhs%cp0   )
   call assign_allocatable(lhs=lhs%cv0 ,  rhs=rhs%cv0   )
   call assign_allocatable(lhs=lhs%alph,  rhs=rhs%alph  )
   call assign_allocatable(lhs=lhs%beta,  rhs=rhs%beta  )
   call assign_allocatable(lhs=lhs%gamm,  rhs=rhs%gamm  )
   call assign_allocatable_gpu(lhs=lhs%cp0_gpu,    rhs=rhs%cp0_gpu    )
   call assign_allocatable_gpu(lhs=lhs%cv0_gpu,    rhs=rhs%cv0_gpu    )
   call assign_allocatable_gpu(lhs=lhs%f_gpu,      rhs=rhs%f_gpu      )
   call assign_allocatable_gpu(lhs=lhs%fp_gpu,     rhs=rhs%fp_gpu     )
   call assign_allocatable_gpu(lhs=lhs%fm_gpu,     rhs=rhs%fm_gpu     )
   call assign_allocatable_gpu(lhs=lhs%dxyz_gpu,   rhs=rhs%dxyz_gpu   )
   call assign_allocatable_gpu(lhs=lhs%alph_gpu,   rhs=rhs%alph_gpu   )
   call assign_allocatable_gpu(lhs=lhs%beta_gpu,   rhs=rhs%beta_gpu   )
   call assign_allocatable_gpu(lhs=lhs%gamm_gpu,   rhs=rhs%gamm_gpu   )
   call assign_allocatable_gpu(lhs=lhs%q_aux_gpu,  rhs=rhs%q_aux_gpu  )
   call assign_allocatable_gpu(lhs=lhs%q_gpu,      rhs=rhs%q_gpu      )
   call assign_allocatable_gpu(lhs=lhs%q_s_gpu,    rhs=rhs%q_s_gpu    )
   endsubroutine eq_assign_eq

   ! non TBP cuf procedures
   subroutine advance_q_gpu_cuf(ni, nj, nk, ngc, nv, nrk, blocks_number, beta_gpu, dt, q_s_gpu, q_gpu)
   !< Advance q_gpu by means of RK stages.
   integer(I4P), intent(in)            :: ni                                     !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                     !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                     !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                    !< Ghost grid number.
   integer(I4P), intent(in)            :: nv                                     !< Number of conservative varibales.
   integer(I4P), intent(in)            :: nrk                                    !< Number of RK stages.
   integer(I4P), intent(in)            :: blocks_number                          !< Number of blocks.
   real(R8P),    intent(in),    device :: beta_gpu(:)                            !< RK betaa coefficients.
   real(R8P),    intent(in)            :: Dt                                     !< Time step.
   real(R8P),    intent(in),    device :: q_s_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:) !< RK stage.
   real(R8P),    intent(inout), device ::   q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative variables.
   integer(I4P)                        :: i, j, k, b, s, v                       !< Counter.
   integer(I4P)                        :: iercuda                                !< Error trapping flag for CUDAFortran.

   do s=1, nrk
      !$cuf kernel do(5) <<<*,*>>>
      do v=1, nv
         do k=1-ngc, nk+ngc
            do j=1-ngc, nj+ngc
               do i=1-ngc, ni+ngc
                  do b=1, blocks_number
                     q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + q_s_gpu(b,i,j,k,v,s) * dt * beta_gpu(s)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   enddo
   endsubroutine advance_q_gpu_cuf

   subroutine compute_aux_cuf(ni, nj, nk, ngc, ns, blocks_number, cp0_gpu, cv0_gpu, q_gpu, q_aux_gpu)
   !< Compute auxiliary variables by means of CUF threads.
   integer(I4P), intent(in)          :: ni                                    !< Grid cells number in I direction.
   integer(I4P), intent(in)          :: nj                                    !< Grid cells number in J direction.
   integer(I4P), intent(in)          :: nk                                    !< Grid cells number in K direction.
   integer(I4P), intent(in)          :: ngc                                   !< Ghost cells number.
   integer(I4P), intent(in)          :: ns                                    !< Number of fluid species.
   integer(I4P), intent(in)          :: blocks_number                         !< Number of blocks.
   real(R8P),    intent(in),  device :: cp0_gpu(:)                            !< Specific heat at constant pressure of species.
   real(R8P),    intent(in),  device :: cv0_gpu(:)                            !< Specific heat at constant pressure of species.
   real(R8P),    intent(in),  device ::     q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(out), device :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   integer(I4P)                      :: b, i, j, k, s                         !< Counter.
   real(R8P)                         :: cp, cv                                !< Specific heats.
   integer(I4P)                      :: iercuda                               !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(4) <<<*,*>>>
   do k=1-ngc, nk+ngc
      do j=1-ngc, nj+ngc
         do i=1-ngc, ni+ngc
            do b=1, blocks_number
               q_aux_gpu(b,i,j,k,ns+1) = 0._R8P
               do s=1, ns
                  q_aux_gpu(b,i,j,k,ns+1) = q_aux_gpu(b,i,j,k,ns+1) + q_gpu(b,i,j,k,s)
               enddo
               do s=1, ns
                  q_aux_gpu(b,i,j,k,s) = q_gpu(b,i,j,k,s) / q_aux_gpu(b,i,j,k,ns+1)
               enddo
               q_aux_gpu(b,i,j,k,ns+2) = q_gpu(b,i,j,k,ns+1) / q_aux_gpu(b,i,j,k,ns+1)
               q_aux_gpu(b,i,j,k,ns+3) = q_gpu(b,i,j,k,ns+2) / q_aux_gpu(b,i,j,k,ns+1)
               q_aux_gpu(b,i,j,k,ns+4) = q_gpu(b,i,j,k,ns+3) / q_aux_gpu(b,i,j,k,ns+1)
               cp = 0._R8P
               cv = 0._R8P
               do s=1, ns
                  cp = cp + q_aux_gpu(b,i,j,k,s) * cp0_gpu(s)
                  cv = cv + q_aux_gpu(b,i,j,k,s) * cv0_gpu(s)
               enddo
               q_aux_gpu(b,i,j,k,ns+5) = cp / cv
               q_aux_gpu(b,i,j,k,ns+6) = (q_gpu(b,i,j,k,ns+4) - 0.5_R8P * q_aux_gpu(b,i,j,k,ns+1) *      &
                                                                         (q_aux_gpu(b,i,j,k,ns+2)**2 +   &
                                                                          q_aux_gpu(b,i,j,k,ns+3)**2 +   &
                                                                          q_aux_gpu(b,i,j,k,ns+4)**2)) * &
                                         (q_aux_gpu(b,i,j,k,ns+5) - 1._R8P)
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_aux_cuf

   subroutine compute_residuals_gpu_cuf(ni, nj, nk, ngc, ns, blocks_number, null_x, null_y, null_z, &
                                        dx_gpu, dy_gpu, dz_gpu, q_aux_gpu, weno_stencils_gpu,       &
                                        fp_gpu, fm_gpu, f_i_gpu, f_j_gpu, f_k_gpu, f_gpu, q_gpu)
   !< Compute residuals of equation.
   integer(I4P), intent(in)            :: ni                                    !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                    !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                    !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                   !< Ghost cells number.
   integer(I4P), intent(in)            :: ns                                    !< Number of species.
   integer(I4P), intent(in)            :: blocks_number                         !< Number of blocks.
   logical,      intent(in)            :: null_x                                !< Nullify x direction.
   logical,      intent(in)            :: null_y                                !< Nullify y direction.
   logical,      intent(in)            :: null_z                                !< Nullify z direction.
   real(R8P),    intent(in),    device :: dx_gpu(1:)                            !< X space steps.
   real(R8P),    intent(in),    device :: dy_gpu(1:)                            !< Y space steps.
   real(R8P),    intent(in),    device :: dz_gpu(1:)                            !< Z space steps.
   real(R8P),    intent(in),    device :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   integer(I4P), intent(in),    device :: weno_stencils_gpu                     !< WENO stencils number/dimension.
   real(R8P),    intent(inout), device ::    fp_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Positive fluxes.
   real(R8P),    intent(inout), device ::    fm_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Negative fluxes.
   real(R8P),    intent(inout), device ::  f_i_gpu(1-ngc:,1:,1-ngc:,1-ngc:,1:)  !< Fluxes for i direction.
   real(R8P),    intent(inout), device ::  f_j_gpu(1-ngc:,1:,1-ngc:,1-ngc:,1:)  !< Fluxes for i direction.
   real(R8P),    intent(inout), device ::  f_k_gpu(1-ngc:,1:,1-ngc:,1-ngc:,1:)  !< Fluxes for i direction.
   real(R8P),    intent(inout), device ::  f_gpu(1:,0:,0:,0:,1:)                !< Convective fluxes.
   real(R8P),    intent(inout), device :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Conservative variables.
   integer(I4P)                        :: b, i, j, k, v                         !< Counter.
   integer(I4P)                        :: iercuda                               !< Error trapping flag for CUDAFortran.

   associate(s=>weno_stencils_gpu)
   ! initialize residuals
   !$cuf kernel do(5) <<<*,*>>>
   do v=1, ns+4
      do k=1, nk
         do j=1, nj
            do i=1,ni
               do b=1, blocks_number
                  q_gpu(b,i,j,k,v) = 0._R8P
               enddo
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()

   ! accumulate difference of fluxes in i direction
   if (.not.null_x) then
      !$cuf kernel do(4) <<<*,*>>>
      do k=1, nk
         do j=1, nj
            do i=1-ngc, ni+ngc
               do b=1, blocks_number
                  call fluxes_pm(r=q_aux_gpu(    b,i,j,k,ns+1), &
                                 u=q_aux_gpu(    b,i,j,k,ns+2), &
                                 v=q_aux_gpu(    b,i,j,k,ns+3), &
                                 w=q_aux_gpu(    b,i,j,k,ns+4), &
                                 g=q_aux_gpu(    b,i,j,k,ns+5), &
                                 p=q_aux_gpu(    b,i,j,k,ns+6), &
                                 fp_rho  =fp_gpu(b,i,j,k,ns  ), &
                                 fp_rho_u=fp_gpu(b,i,j,k,ns+1), &
                                 fp_rho_v=fp_gpu(b,i,j,k,ns+2), &
                                 fp_rho_w=fp_gpu(b,i,j,k,ns+3), &
                                 fp_rho_E=fp_gpu(b,i,j,k,ns+4), &
                                 fm_rho  =fm_gpu(b,i,j,k,ns  ), &
                                 fm_rho_u=fm_gpu(b,i,j,k,ns+1), &
                                 fm_rho_v=fm_gpu(b,i,j,k,ns+2), &
                                 fm_rho_w=fm_gpu(b,i,j,k,ns+3), &
                                 fm_rho_E=fm_gpu(b,i,j,k,ns+4))
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      !$cuf kernel do(5) <<<*,*>>>
      do v=1, ns+4
         do k=1, nk
            do j=1, nj
               do i=1-ngc, ni+ngc
                  do b=1, blocks_number
                     f_i_gpu(i,b,j,k,v) = fp_gpu(b,i,j,k,v)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      !$cuf kernel do(4) <<<*,*>>>
      do v=1, ns+4
         do k=1, nk
            do j=1, nj
               do b=1, blocks_number
                  do i=0,ni
                     call reconstruct_weno(side=weno_l_side,               &
                                           s=s,                            &
                                           q=f_i_gpu(i+1-s:i-1+s,b,j,k,v), &
                                           qr=f_gpu(b,i,j,k,v))
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      !$cuf kernel do(5) <<<*,*>>>
      do v=1, ns+4
         do k=1, nk
            do j=1, nj
               do i=1,ni
                  do b=1, blocks_number
                     q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + (f_gpu(b,i-1,j,k,v) - f_gpu(b,i,j,k,v)) / dx_gpu(b)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      !$cuf kernel do(5) <<<*,*>>>
      do v=1, ns+4
         do k=1, nk
            do j=1, nj
               do i=1-ngc, ni+ngc
                  do b=1, blocks_number
                     f_i_gpu(i,b,j,k,v) = fm_gpu(b,i,j,k,v)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      !$cuf kernel do(4) <<<*,*>>>
      do v=1, ns+4
         do k=1, nk
            do j=1, nj
               do b=1, blocks_number
                  do i=0,ni
                     call reconstruct_weno(side=weno_r_side,                   &
                                           s=s,                                &
                                           q=f_i_gpu(i+1-s+1:i-1+s+1,b,j,k,v), &
                                           qr=f_gpu(b,i,j,k,v))
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      !$cuf kernel do(5) <<<*,*>>>
      do v=1, ns+4
         do k=1, nk
            do j=1, nj
               do i=1,ni
                  do b=1, blocks_number
                     q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + (f_gpu(b,i-1,j,k,v) - f_gpu(b,i,j,k,v)) / dx_gpu(b)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   endif

   ! accumulate difference of fluxes in j direction
   if (.not.null_y) then
      !$cuf kernel do(4) <<<*,*>>>
      do k=1, nk
         do j=1-ngc, nj+ngc
            do i=1,ni
               do b=1, blocks_number
                  call fluxes_pm(r=q_aux_gpu(    b,i,j,k,ns+1), &
                                 u=q_aux_gpu(    b,i,j,k,ns+3), &
                                 v=q_aux_gpu(    b,i,j,k,ns+2), &
                                 w=q_aux_gpu(    b,i,j,k,ns+4), &
                                 g=q_aux_gpu(    b,i,j,k,ns+5), &
                                 p=q_aux_gpu(    b,i,j,k,ns+6), &
                                 fp_rho  =fp_gpu(b,i,j,k,ns  ), &
                                 fp_rho_u=fp_gpu(b,i,j,k,ns+2), &
                                 fp_rho_v=fp_gpu(b,i,j,k,ns+1), &
                                 fp_rho_w=fp_gpu(b,i,j,k,ns+3), &
                                 fp_rho_E=fp_gpu(b,i,j,k,ns+4), &
                                 fm_rho  =fm_gpu(b,i,j,k,ns  ), &
                                 fm_rho_u=fm_gpu(b,i,j,k,ns+2), &
                                 fm_rho_v=fm_gpu(b,i,j,k,ns+1), &
                                 fm_rho_w=fm_gpu(b,i,j,k,ns+3), &
                                 fm_rho_E=fm_gpu(b,i,j,k,ns+4))
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      !$cuf kernel do(5) <<<*,*>>>
      do v=1, ns+4
         do k=1, nk
            do j=1-ngc, nj+ngc
               do i=1, ni
                  do b=1, blocks_number
                     f_j_gpu(j,b,i,k,v) = fp_gpu(b,i,j,k,v)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      !$cuf kernel do(4) <<<*,*>>>
      do v=1, ns+4
         do k=1, nk
            do i=1,ni
               do b=1, blocks_number
                  do j=0, nj
                     call reconstruct_weno(side=weno_l_side,               &
                                           s=s,                            &
                                           q=f_j_gpu(j+1-s:j-1+s,b,i,k,v), &
                                           qr=f_gpu(b,i,j,k,v))
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      !$cuf kernel do(5) <<<*,*>>>
      do v=1, ns+4
         do k=1, nk
            do j=1, nj
               do i=1,ni
                  do b=1, blocks_number
                     q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + (f_gpu(b,i,j-1,k,v) - f_gpu(b,i,j,k,v)) / dy_gpu(b)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      !$cuf kernel do(5) <<<*,*>>>
      do v=1, ns+4
         do k=1, nk
            do j=1-ngc, nj+ngc
               do i=1, ni
                  do b=1, blocks_number
                     f_j_gpu(j,b,i,k,v) = fm_gpu(b,i,j,k,v)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      !$cuf kernel do(5) <<<*,*>>>
      do v=1, ns+4
         do k=1, nk
            do i=1, ni
               do b=1, blocks_number
                  do j=0,nj
                     call reconstruct_weno(side=weno_r_side,                   &
                                           s=s,                                &
                                           q=f_j_gpu(j+1-s+1:j-1+s+1,b,i,k,v), &
                                           qr=f_gpu(b,i,j,k,v))
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      !$cuf kernel do(5) <<<*,*>>>
      do v=1, ns+4
         do k=1, nk
            do j=1, nj
               do i=1,ni
                  do b=1, blocks_number
                     q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + (f_gpu(b,i,j-1,k,v) - f_gpu(b,i,j,k,v)) / dy_gpu(b)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   endif

   ! accumulate difference of fluxes in k direction
   if (.not.null_z) then
      !$cuf kernel do(4) <<<*,*>>>
      do k=1-ngc, nk+ngc
         do j=1, nj
            do i=1,ni
               do b=1, blocks_number
                  call fluxes_pm(r=q_aux_gpu(    b,i,j,k,ns+1), &
                                 u=q_aux_gpu(    b,i,j,k,ns+4), &
                                 v=q_aux_gpu(    b,i,j,k,ns+2), &
                                 w=q_aux_gpu(    b,i,j,k,ns+3), &
                                 g=q_aux_gpu(    b,i,j,k,ns+5), &
                                 p=q_aux_gpu(    b,i,j,k,ns+6), &
                                 fp_rho  =fp_gpu(b,i,j,k,ns  ), &
                                 fp_rho_u=fp_gpu(b,i,j,k,ns+3), &
                                 fp_rho_v=fp_gpu(b,i,j,k,ns+1), &
                                 fp_rho_w=fp_gpu(b,i,j,k,ns+2), &
                                 fp_rho_E=fp_gpu(b,i,j,k,ns+4), &
                                 fm_rho  =fm_gpu(b,i,j,k,ns  ), &
                                 fm_rho_u=fm_gpu(b,i,j,k,ns+3), &
                                 fm_rho_v=fm_gpu(b,i,j,k,ns+1), &
                                 fm_rho_w=fm_gpu(b,i,j,k,ns+2), &
                                 fm_rho_E=fm_gpu(b,i,j,k,ns+4))
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      !$cuf kernel do(5) <<<*,*>>>
      do v=1, ns+4
         do k=1-ngc, nk+ngc
            do j=1, nj
               do i=1, ni
                  do b=1, blocks_number
                     f_k_gpu(k,b,i,j,v) = fp_gpu(b,i,j,k,v)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      !$cuf kernel do(4) <<<*,*>>>
      do v=1, ns+4
         do j=1, nj
            do i=1,ni
               do b=1, blocks_number
                  do k=0, nk
                     call reconstruct_weno(side=weno_l_side,               &
                                           s=s,                            &
                                           q=f_k_gpu(k+1-s:k-1+s,b,i,j,v), &
                                           qr=f_gpu(b,i,j,k,v))
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      !$cuf kernel do(5) <<<*,*>>>
      do v=1, ns+4
         do k=1, nk
            do j=1, nj
               do i=1,ni
                  do b=1, blocks_number
                     q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + (f_gpu(b,i,j,k-1,v) - f_gpu(b,i,j,k,v)) / dz_gpu(b)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      !$cuf kernel do(5) <<<*,*>>>
      do v=1, ns+4
         do k=1-ngc, nk+ngc
            do j=1, nj
               do i=1, ni
                  do b=1, blocks_number
                     f_k_gpu(k,b,i,j,v) = fm_gpu(b,i,j,k,v)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      !$cuf kernel do(5) <<<*,*>>>
      do v=1, ns+4
         do j=1, nj
            do i=1,ni
               do b=1, blocks_number
                  do k=0, nk
                     call reconstruct_weno(side=weno_r_side,                   &
                                           s=s,                                &
                                           q=f_k_gpu(k+1-s+1:k-1+s+1,b,i,j,v), &
                                           qr=f_gpu(b,i,j,k,v))
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      !$cuf kernel do(5) <<<*,*>>>
      do v=1, ns+4
         do k=1, nk
            do j=1, nj
               do i=1,ni
                  do b=1, blocks_number
                     q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + (f_gpu(b,i,j,k-1,v) - f_gpu(b,i,j,k,v)) / dz_gpu(b)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   endif
   endassociate
   endsubroutine compute_residuals_gpu_cuf

   subroutine compute_rk_stage_gpu_cuf(ni, nj, nk, ngc, nv, blocks_number, alph_gpu, dt, s, q_gpu, q_s_gpu)
   !< Initialize RK stage with q_gpu.
   integer(I4P), intent(in)            :: ni                                     !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                     !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                     !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                    !< Ghost cells number.
   integer(I4P), intent(in)            :: nv                                     !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number                          !< Number of blocks.
   real(R8P),    intent(in),    device :: alph_gpu(:,:)                          !< RK alpha coefficients.
   real(R8P),    intent(in)            :: dt                                     !< Time step.
   integer(I4P), intent(in)            :: s                                      !< Stage to initialize.
   real(R8P),    intent(in),    device ::   q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative field.
   real(R8P),    intent(inout), device :: q_s_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:) !< RK stage.
   integer(I4P)                        :: i, j, k, b, v, ss                      !< Counter.
   integer(I4P)                        :: iercuda                                !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(5) <<<*,*>>>
   do v=1, nv
      do k=1, nk
         do j=1, nj
            do i=1, ni
               do b=1, blocks_number
                  q_s_gpu(b,i,j,k,v,s) = q_gpu(b,i,j,k,v)
               enddo
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   do ss=1, s - 1
      !$cuf kernel do(5) <<<*,*>>>
      do v=1, nv
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do b=1, blocks_number
                     q_s_gpu(b,i,j,k,v,s) = q_s_gpu(b,i,j,k,v,s) + (q_s_gpu(b,i,j,k,v,ss) * (dt * alph_gpu(s, ss)))
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   enddo
   endsubroutine compute_rk_stage_gpu_cuf

   subroutine compute_umax_cuf(b, ni, nj, nk, ngc, ns, dx, dy, dz, q_aux_gpu, umax)
   !< Compute maximum speed by means of CUF threads.
   integer(I4P), intent(in)         :: b                                     !< Block index.
   integer(I4P), intent(in)         :: ni                                    !< Grid cells number in I direction.
   integer(I4P), intent(in)         :: nj                                    !< Grid cells number in J direction.
   integer(I4P), intent(in)         :: nk                                    !< Grid cells number in K direction.
   integer(I4P), intent(in)         :: ngc                                   !< Ghost cells number.
   integer(I4P), intent(in)         :: ns                                    !< Number of species.
   real(R8P),    intent(in)         :: dx                                    !< X space step.
   real(R8P),    intent(in)         :: dy                                    !< Y space step.
   real(R8P),    intent(in)         :: dz                                    !< Z space step.
   real(R8P),    intent(in), device :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary varibales.
   real(R8P),    intent(out)        :: umax                                  !< Maximum speed.
   real(R8P)                        :: ss                                    !< Speed of sound.
   integer(I4P)                     :: i, j, k                               !< Counter.
   integer(I4P)                     :: iercuda                               !< Error trapping flag for CUDAFortran.

   umax = 0._R8P
   !$cuf kernel do(3) <<<*,*>>>
   do k=1, nk
      do j=1, nj
         do i=1, ni
            ss = a(p=q_aux_gpu(b,i,j,k,ns+6), r=q_aux_gpu(b,i,j,k,ns+1), g=q_aux_gpu(b,i,j,k,ns+5))
            umax = max(umax, abs(q_aux_gpu(b,i,j,k,ns+2)) + ss, &
                             abs(q_aux_gpu(b,i,j,k,ns+3)) + ss, &
                             abs(q_aux_gpu(b,i,j,k,ns+4)) + ss)

         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_umax_cuf

   ! non type-bound kernel procedures
   attributes(device) subroutine fluxes_pm(r, u, v, w, g, p,                              &
                                           fp_rho, fp_rho_u, fp_rho_v, fp_rho_w, fp_rho_E,&
                                           fm_rho, fm_rho_u, fm_rho_v, fm_rho_w, fm_rho_E)
   !< Compute positive and negative fluxes in cell center.
   real(R8P), intent(in)  :: r             !< Density.
   real(R8P), intent(in)  :: u             !< Normal velocity.
   real(R8P), intent(in)  :: v             !< First tangential component of velocity.
   real(R8P), intent(in)  :: w             !< Second tangential component of velocity.
   real(R8P), intent(in)  :: p             !< Pressure.
   real(R8P), intent(in)  :: g             !< Specific heats ratio.
   real(R8P), intent(out) :: fp_rho        !< Positive flux of rho.
   real(R8P), intent(out) :: fp_rho_u      !< Positive Flux of rho*u.
   real(R8P), intent(out) :: fp_rho_v      !< Positive Flux of rho*v.
   real(R8P), intent(out) :: fp_rho_w      !< Positive Flux of rho*w.
   real(R8P), intent(out) :: fp_rho_E      !< Positive Flux of rho*E.
   real(R8P), intent(out) :: fm_rho        !< Negative flux of rho.
   real(R8P), intent(out) :: fm_rho_u      !< Negative Flux of rho*u.
   real(R8P), intent(out) :: fm_rho_v      !< Negative Flux of rho*v.
   real(R8P), intent(out) :: fm_rho_w      !< Negative Flux of rho*w.
   real(R8P), intent(out) :: fm_rho_E      !< Negative Flux of rho*E.
   real(R8P)              :: c2, c, Hh     !< Dummy.
   real(R8P)              :: alfa_lam      !< Dummy.

   c2 = a2(p=p, r=r, g=g)
   c  = sqrt(c2)
   Hh = Hv2(p=p, r=r, u2=(u**2+v**2+w**2), g=g)

   if (u+c < 0._R8P) then
      fm_rho   = r*u
      fm_rho_u = fm_rho*u + p
      fm_rho_v = fm_rho*v
      fm_rho_w = fm_rho*w
      fm_rho_E = fm_rho*Hh
      fp_rho   = 0._R8P
      fp_rho_u = 0._R8P
      fp_rho_v = 0._R8P
      fp_rho_w = 0._R8P
      fp_rho_E = 0._R8P
   else if (u-c > 0d0) then
      fp_rho   = r*u
      fp_rho_u = fp_rho*u + p
      fp_rho_v = fp_rho*v
      fp_rho_w = fp_rho*w
      fp_rho_E = fp_rho*Hh
      fm_rho   = 0._R8P
      fm_rho_u = 0._R8P
      fm_rho_v = 0._R8P
      fm_rho_w = 0._R8P
      fm_rho_E = 0._R8P
   else
      alfa_lam = 0.5_R8P*r/g
      if (u > 0._R8P) then
         alfa_lam = alfa_lam*(u-c)

         fm_rho   = alfa_lam
         fm_rho_u = alfa_lam*(u-c)
         fm_rho_v = alfa_lam*v
         fm_rho_w = alfa_lam*w
         fm_rho_E = alfa_lam*(Hh + c*u)

         fp_rho    = r*u
         fp_rho_u  = fp_rho*u + p
         fp_rho_v  = fp_rho*v
         fp_rho_w  = fp_rho*w
         fp_rho_E  = fp_rho*Hh

         fp_rho    = fp_rho   - fm_rho
         fp_rho_u  = fp_rho_u - fm_rho_u
         fp_rho_v  = fp_rho_v - fm_rho_v
         fp_rho_w  = fp_rho_w - fm_rho_w
         fp_rho_E  = fp_rho_E - fm_rho_E
      else
         alfa_lam = alfa_lam*(u+c)

         fp_rho   = alfa_lam
         fp_rho_u = alfa_lam*(u+c)
         fp_rho_v = alfa_lam*v
         fp_rho_w = alfa_lam*w
         fp_rho_E = alfa_lam*(Hh - c*u)

         fm_rho    = r*u
         fm_rho_u  = fm_rho*u + p
         fm_rho_v  = fm_rho*v
         fm_rho_w  = fm_rho*w
         fm_rho_E  = fm_rho*Hh

         fm_rho    = fm_rho   - fp_rho
         fm_rho_u  = fm_rho_u - fp_rho_u
         fm_rho_v  = fm_rho_v - fp_rho_v
         fm_rho_w  = fm_rho_w - fp_rho_w
         fm_rho_E  = fm_rho_E - fp_rho_E
      endif
   endif
   endsubroutine fluxes_pm

   attributes(device) function a(p, r, g) result(ss)
   !< Return speed of sound for an ideal calorically perfect gas.
   real(R8P), intent(in) :: p  !< Pressure.
   real(R8P), intent(in) :: r  !< Density.
   real(R8P), intent(in) :: g  !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
   real(R8P)             :: ss !< Speed of sound.

   ss = sqrt(g*p/r)
   endfunction a

   attributes(device) function a2(p, r, g) result(ss)
   !< Return square of speed of sound for an ideal calorically perfect gas.
   real(R8P), intent(in) :: p  !< Pressure.
   real(R8P), intent(in) :: r  !< Density.
   real(R8P), intent(in) :: g  !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
   real(R8P)             :: ss !< Speed of sound.

   ss = g*p/r
   endfunction a2

   attributes(host,device) function E(p, r, u, g) result(energy)
   !< Return total specific energy (per unit of mass).
   !<$$
   !<  E = \frac{p}{{\left( {\g  - 1} \right)\r }} + \frac{{u^2 }}{2}
   !<$$
   real(R8P), intent(in) :: p      !< Pressure.
   real(R8P), intent(in) :: r      !< Density.
   real(R8P), intent(in) :: u      !< Module of velocity vector.
   real(R8P), intent(in) :: g      !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
   real(R8P)             :: energy !< Total specific energy (per unit of mass).

   energy = p/((g - 1._R8P) * r) + 0.5_R8P * u*u
   endfunction E

   attributes(device) function H(p, r, u, g) result(entalpy)
   !< Return total specific entalpy (per unit of mass).
   !<$$
   !<  H = \frac{{\g p}}{{\left( {\g  - 1} \right)\r }} + \frac{{u^2 }}{2}
   !<$$
   real(R8P), intent(in) :: g       !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
   real(R8P), intent(in) :: p       !< Pressure.
   real(R8P), intent(in) :: r       !< Density.
   real(R8P), intent(in) :: u       !< Module of velocity vector.
   real(R8P)             :: entalpy !< Total specific entalpy (per unit of mass).

   entalpy = g * p / ((g - 1._R_P) * r) + 0.5_R_P * u*u
   endfunction H

   attributes(device) function Hv2(p, r, u2, g) result(entalpy)
   !< Return total specific entalpy (per unit of mass).
   !<$$
   !<  H = \frac{{\g p}}{{\left( {\g  - 1} \right)\r }} + \frac{{u^2 }}{2}
   !<$$
   real(R8P), intent(in) :: g       !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
   real(R8P), intent(in) :: p       !< Pressure.
   real(R8P), intent(in) :: r       !< Density.
   real(R8P), intent(in) :: u2      !< Square of module of velocity vector.
   real(R8P)             :: entalpy !< Total specific entalpy (per unit of mass).

   entalpy = g * p / ((g - 1._R_P) * r) + 0.5_R_P * u2
   endfunction Hv2

   attributes(device) function p(r, a, g) result(pressure)
   !< Return pressure for an ideal calorically perfect gas.
   real(R8P), intent(in) :: r        !< Density.
   real(R8P), intent(in) :: a        !< Speed of sound.
   real(R8P), intent(in) :: g        !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
   real(R8P)             :: pressure !< Pressure.

   pressure = r*a*a/g
   endfunction p

   attributes(device) function r(p, a, g) result(density)
   !< Return density for an ideal calorically perfect gas.
   real(R8P), intent(in) :: p       !< Pressure.
   real(R8P), intent(in) :: a       !< Speed of sound.
   real(R8P), intent(in) :: g       !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
   real(R8P)             :: density !< Density.

   density = g*p/(a*a)
   endfunction r
endmodule adam_equation_euler_gpu_object
