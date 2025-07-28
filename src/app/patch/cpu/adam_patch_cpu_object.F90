!< ADAM, PATCH Poisson solver with Adpative mesh Refinement for HPC computing class definition, CPU backend.
module adam_patch_cpu_object
!< ADAM, PATCH Poisson solver with Adpative mesh Refinement for HPC computing class definition, CPU backend.

! ADAM modules
use adam_common_library
! patch modules
use adam_patch_common_library
! third party modules
use penf
use mpi

implicit none
private
public :: patch_cpu_object

type, extends(patch_common_object) :: patch_cpu_object
   !< Maxwell equations system class definition, CPU backend.
   real(R8P), allocatable :: dq(:,:,:,:,:) !< Residuals right hand side.
   contains
      ! auxiliary methods
      procedure, pass(self) :: allocate_cpu !< Allocate CPU data.
      procedure, pass(self) :: finalize     !< Finalize simulation.
      procedure, pass(self) :: initialize   !< Initialize the equation.
      ! AMR methods
      procedure, pass(self) :: amr_update       !< Do AMR update.
      procedure, pass(self) :: compute_phi      !< Compute phi, distance from IB solid.
      procedure, pass(self) :: mark_by_geo      !< Mark blocks to be refined/derefined by a geometric constrain.
      procedure, pass(self) :: mark_by_grad_var !< Mark blocks to be refined/derefined by a `grad(var)` value.
      procedure, pass(self) :: move_phi         !< Move phi.
      procedure, pass(self) :: refine_uniform   !< Refine all blocks uniformly.
      ! IB methods
      procedure, pass(self) :: integrate_eikonal !< Integrate eikonal equation.
      ! IO methods
      procedure, pass(self) :: load_restart_files   !< Load restart files.
      procedure, pass(self) :: save_xh5f            !< Save simulation data in XH5F format.
      procedure, pass(self) :: save_residuals       !< Save residuals history.
      procedure, pass(self) :: save_restart_files   !< Save restart files.
      procedure, pass(self) :: save_simulation_data !< Save all simulation data.
      ! IC/BC
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions of equation.
      procedure, pass(self) :: set_initial_conditions  !< Set initial conditions (and coils) of equation.
      procedure, pass(self) :: update_ghost            !< Update ghost cells and set boundary conditions.
      ! numerical methods
      procedure, pass(self) :: integrate !< Integrate the equation.
      procedure, pass(self) :: simulate  !< Perform the simulation.
endtype patch_cpu_object

interface
   subroutine compute_smoothing_interface(ni, nj, nk, ngc, dx2, dy2, dz2, r, dq, q, dq_max, iterations)
   !< Compute smoothing, abstract interface.
   import :: I4P, R8P
   integer(I4P), intent(in)              :: ni,nj,nk,ngc !< Grid dimensions.
   real(R8P),    intent(in)              :: dx2,dy2,dz2  !< Square space steps.
   real(R8P),    intent(in)              :: r(1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:)    !< Rho distribution.
   real(R8P),    intent(inout), optional :: dq(1-ngc:,&
                                               1-ngc:,&
                                               1-ngc:)   !< Residuals.
   real(R8P),    intent(inout)           :: q(1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:)    !< Poisson field variable.
   real(R8P),    intent(inout), optional :: dq_max       !< Maximum residual.
   integer(I4P), intent(in),    optional :: iterations   !< Smoothing iterations.
   endsubroutine compute_smoothing_interface
endinterface
contains
   ! auxiliary methods
   subroutine allocate_cpu(self)
   !< Allocate CPU data.
   class(patch_cpu_object), intent(inout) :: self !< The equation.
   character(:), allocatable              :: msg_ !< Allocating message base.
   character(:), allocatable              :: msg  !< Allocating message.

   call self%mpih%print_message('patch_cpu_object%allocate_cpu start')
   msg_ = self%mpih%myrankstr//'patch_cpu_object%allocate_cpu '
   associate(nv=>self%nv, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, nb=>self%nb)
   msg = msg_//' dq '
   call allocate_variable(var=self%dq, ulb=reshape([1,nv,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,5]),msg=msg)
   self%dq = 0._R8P
   endassociate
   call self%mpih%print_message('patch_cpu_object%allocate_cpu finish')
   endsubroutine allocate_cpu

   subroutine finalize(self)
   !< Finalize simulation.
   class(patch_cpu_object), intent(inout) :: self !< The equation.

   call self%save_simulation_data
   if (self%mpih%myrank==0) call self%io%close_file_residuals
   call self%mpih%finalize
   endsubroutine finalize

   subroutine initialize(self, filename)
   !< Initialize the equation.
   class(patch_cpu_object), intent(inout) :: self     !< The equation.
   character(*),            intent(in)    :: filename !< Input file name.

   call self%mpih%initialize(do_mpi_init=.true., verbose=.true.)
   call self%mpih%print_message('patch_cpu_object%initialize start')
   call self%initialize_common(field = self%adam%field, filename=filename, memory_avail=self%mpih%memory_avail)
   call self%allocate_cpu
   if (self%io%restart) then
      call self%mpih%print_message('restart simulation from "'//trim(self%io%restart_basename)//'" files')
      call self%load_restart_files(t=self%time%it, time=self%time%time)
      call self%mpih%print_message('restart [t, time]: '//trim(str(self%time%it))//', '//trim(str(self%time%time)))
   else
      call self%mpih%print_message('impose initial conditions start')
      call self%set_initial_conditions
      self%time%time = 0._R8P
      self%time%it = 0
      call self%mpih%print_message('impose initial conditions finish')
   endif
   call self%save_simulation_data
   if (self%mpih%myrank==0) call self%io%open_file_residuals(nv=self%nv)
   call self%mpih%print_message('patch_cpu_object%initialize finish')
   endsubroutine initialize

   ! AMR methods
   subroutine amr_update(self)
   !< Do AMR update.
   class(patch_cpu_object), intent(inout) :: self                !< The equation.
   logical                                :: is_grid_changed     !< Flag to check grid changes for each marker.
   logical                                :: is_grid_changed_all !< Flag to check grid changes for each iter.
   integer(I4P)                           :: i, i_marker         !< Counter.
   type(amr_marker_object)                :: amr_marker          !< Current amr marker.

   amr: do i=1, self%amr%iters
      is_grid_changed_all = .false.
      do i_marker=1, self%amr%markers_number
         amr_marker = self%amr%markers(i_marker)
         call self%update_ghost(q=self%q)
         select case(amr_marker%mode)
         case(AMR_GEO)
            call self%mark_by_geo(delta_fine=amr_marker%delta_fine, delta_coarse=amr_marker%delta_coarse)
         case(AMR_GRAD)
            select case(amr_marker%field)
            case(1)
               call self%mark_by_grad_var(grad_tol=amr_marker%tol, delta_type=amr_marker%delta_type, &
                                          delta_fine=amr_marker%delta_fine,                          &
                                          delta_coarse=amr_marker%delta_coarse, ivar=amr_marker%ivar)
            case(2)
               call self%mark_by_grad_var(grad_tol=amr_marker%tol, delta_type=amr_marker%delta_type, &
                                          delta_fine=amr_marker%delta_fine,                          &
                                          delta_coarse=amr_marker%delta_coarse, ivar=amr_marker%ivar)
            endselect
         endselect
         call self%adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false., is_grid_changed=is_grid_changed, q=self%q)
         if (self%ib%solids_number > 0) call self%compute_phi()
         is_grid_changed_all = is_grid_changed_all.or.is_grid_changed
      enddo
      if (.not.is_grid_changed_all) then
          call self%mpih%print_message('AMR Grid stabilized after : '//trim(str(i))//' AMR iterations')
          exit amr
       elseif (i==self%amr%iters) then
          call self%mpih%print_message('AMR Grid is NOT stabilized after : '//trim(str(i))//' AMR iterations')
      endif
   enddo amr
   endsubroutine amr_update

   subroutine compute_phi(self)
   !< Compute phi, distance from IB solid.
   class(patch_cpu_object), intent(inout) :: self !< The equation.

   if (self%ib%solids_number>0) then
      call self%mpih%print_message('compute IB distance start')
      call self%ib%compute_phi
      call self%mpih%print_message('compute IB distance finish')
   endif
   endsubroutine compute_phi

   subroutine mark_by_geo(self, delta_fine, delta_coarse, threshold, do_init)
   !< Mark blocks to be refined/derefined by a geometric constrain.
   class(patch_cpu_object), intent(inout)        :: self           !< The equation.
   real(R8P),               intent(in)           :: delta_fine     !< Maximum cell delta in fine grids.
   real(R8P),               intent(in)           :: delta_coarse   !< Minimum cell delta in coarse grids.
   real(R8P),               intent(in), optional :: threshold      !< Threshold for sphere proximity.
   real(R8P)                                     :: threshold_     !< Threshold for sphere proximity, local var.
   real(R8P)                                     :: max_cell_delta !< Maximum cell delta.
   real(R8P)                                     :: distance       !< Value (max) of gradient of rho.
   integer(I4P)                                  :: b              !< Counter.
   logical, optional,       intent(in)           :: do_init
   logical                                       :: do_init_

   do_init_ = .true.    ; if (present(do_init)) do_init_ = do_init
   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   if (do_init_) self%field%refinements_needed = [(TO_BE_DEREFINED,b=1,self%blocks_number)]
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
             blocks_number=>self%blocks_number, dxyz=>self%field%dxyz, phi=>self%ib%phi)
      do b=1, blocks_number
         distance = 1._R8P
         if (maxval(phi(1,:,:,:,b))*minval(phi(1,:,:,:,b)) < 0._R8P) then
            distance = 0._R8P
         endif
         max_cell_delta = max_cell_delta_dist(distance=distance)

         if (maxval(dxyz(:,b)) > max_cell_delta) then
            self%field%refinements_needed(b) = TO_BE_REFINED
         elseif (maxval(dxyz(:,b)) * threshold_ < max_cell_delta) then
            self%field%refinements_needed(b) = max(self%field%refinements_needed(b), TO_BE_DEREFINED)
         else
            self%field%refinements_needed(b) = max(self%field%refinements_needed(b), TO_NOT_TOUCH)
         endif
      enddo
   endassociate
   contains
      function max_cell_delta_dist(distance) result(delta)
      !< Return the maximum cell delta given a comparison distance.
      real(R8P),          intent(in) :: distance !< Comparison distance.
      real(R8P)                      :: delta    !< Maximum cell delta admissible.

      if (abs(distance) < epsilon(0._R8P)) then
         delta = delta_fine
      else
         delta = delta_coarse
      endif
      endfunction max_cell_delta_dist
   endsubroutine mark_by_geo

   subroutine mark_by_grad_var(self, grad_tol, delta_type, delta_fine, delta_coarse, ivar, threshold, do_init)
   !< Mark blocks to be refined/derefined by a `grad(var)` value.
   class(patch_cpu_object), intent(inout)        :: self                     !< The equation.
   real(R8P),               intent(in)           :: grad_tol                 !< Gradiend tolerance value.
   character(*),            intent(in)           :: delta_type               !< Delta criterion type.
   real(R8P),               intent(in)           :: delta_fine               !< Maximum cell delta in fine grids.
   real(R8P),               intent(in)           :: delta_coarse             !< Minimum cell delta in coarse grids.
   integer(I4P),            intent(in), optional :: ivar                     !< Variable for marking.
   real(R8P),               intent(in), optional :: threshold                !< Threshold for sphere proximity.
   logical,                 intent(in), optional :: do_init                  !< Re-initialize refinements queries.
   integer(I4P)                                  :: ivar_                    !< Variable for marking (local var).
   logical                                       :: do_init_                 !< Re-initialize refinements queries, local var.
   real(R8P)                                     :: threshold_               !< Threshold for sphere proximity, local var.
   real(R8P)                                     :: max_cell_delta           !< Maximum cell delta.
   real(R8P)                                     :: grad_var                 !< Value (max) of gradient of var.
   integer(I4P)                                  :: b                        !< Counter.
   real(R8P)                                     :: dc(1:self%blocks_number) !< Delta criterion.

   ! ivar_     = 1_R4P    ; if (present(ivar)) ivar_ = ivar
   ! do_init_ = .true.    ; if (present(do_init)) do_init_ = do_init
   ! threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   ! if (do_init_) self%field%refinements_needed = [(TO_BE_DEREFINED,b=1,self%blocks_number)]
   ! associate (ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
   !            blocks_number=>self%blocks_number, dxyz=>self%field%dxyz)
   !    select case(delta_type)
   !    case(AMR_DELTA_T_X)
   !       dc(1:blocks_number) = dxyz(1,1:blocks_number)
   !    case(AMR_DELTA_T_Y)
   !       dc(1:blocks_number) = dxyz(2,1:blocks_number)
   !    case(AMR_DELTA_T_Z)
   !       dc(1:blocks_number) = dxyz(3,1:blocks_number)
   !    case(AMR_DELTA_T_MAX)
   !       do b=1, blocks_number
   !          dc(b) = maxval(dxyz(:,b))
   !       enddo
   !    endselect
   !    call self%update_ghost(q=self%q)
   !    ! call self%compute_q_auxiliary(q=self%q, q_aux=self%q_aux)
   !    do b=1, blocks_number
   !       call compute_q_gradient(b=b, ni=ni, nj=nj, nk=nk, ngc=ngc, &
   !                               dx=dxyz(1,b), dy=dxyz(2,b), dz=dxyz(3,b), q=self%q_aux, ivar=ivar_, gradient=grad_var)
   !       max_cell_delta = max_cell_delta_grad(grad=grad_var)
   !       if ((dc(b)) > max_cell_delta) then
   !          self%field%refinements_needed(b) = TO_BE_REFINED
   !       elseif ((dc(b)) * threshold_ < max_cell_delta) then
   !          self%field%refinements_needed(b) = max(self%field%refinements_needed(b), TO_BE_DEREFINED)
   !       else
   !          self%field%refinements_needed(b) = max(self%field%refinements_needed(b), TO_NOT_TOUCH)
   !       endif
   !    enddo
   ! endassociate
   contains
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
   endsubroutine mark_by_grad_var

   subroutine move_phi(self, velocity, s)
   !< Move phi and the actual tree representation.
   class(patch_cpu_object), intent(inout) :: self        !< The equation.
   real(R8P),               intent(in)    :: velocity(3) !< Velocity of the movement.
   integer(I4P),            intent(in)    :: s           !< Solid index.

   if (self%ib%solids_number>0) then
      call self%mpih%print_message('move IB distance start')
      call self%ib%move_phi(velocity=velocity, s=s)
      call self%mpih%print_message('move IB distance finish')
   endif
   endsubroutine move_phi

   subroutine refine_uniform(self, refinement_levels)
   !< Refine all blocks uniformly.
   class(patch_cpu_object), intent(inout) :: self              !< The equation.
   integer(I4P),            intent(in)    :: refinement_levels !< Number of refinement to be performed.
   integer(I4P)                           :: l                 !< Counter.

   do l=1, refinement_levels
      call self%adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
      call self%adam%amr_update(do_blocks_reorder=.false., do_mpi_redistribute=.true., q=self%q)
   enddo
   endsubroutine

   ! IB methods
   subroutine integrate_eikonal(self, q)
   !< Integrate eikonal equation.
   class(patch_cpu_object), intent(inout) :: self      !< The equation.
   real(R8P),               intent(inout) :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1:)     !< Conservative variables.
   integer(I4P)                           :: i_eikonal !< Counter.

   associate(blocks_number=>self%blocks_number, solids_number=>self%ib%solids_number)
      if (blocks_number > 0) then
         if (solids_number > 0) then
            call self%update_ghost(q=q)
            do i_eikonal=1, self%ib%n_eikonal
               call self%mpih%barrier
               call self%ib%evolve_eikonal(q=q)
               call self%update_ghost(q=q)
            enddo
            call self%ib%invert_eikonal(q=q)
            call self%mpih%barrier
         endif
      endif
   endassociate
   endsubroutine integrate_eikonal

   ! IO methods
   subroutine load_restart_files(self, t, time)
   !< Save restart files.
   class(patch_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P),            intent(out)   :: t    !< Time iteration.
   real(R8P),               intent(out)   :: time !< Time.

   call self%adam%load_restart_files(basename=self%io%restart_basename, t=t, time=time, q=self%q)
   call self%adam%make_comm_local_maps_ghost_bc
   endsubroutine load_restart_files

   subroutine save_xh5f(self, output_basename, with_ghost)
   !< Save simulation data in HDF5 format.
   class(patch_cpu_object), intent(inout)        :: self             !< The equation.
   character(*),            intent(in), optional :: output_basename  !< Output basename.
   logical,                 intent(in), optional :: with_ghost       !< Flag to save ghost cells.
   character(:), allocatable                     :: output_basename_ !< Output basename, local var.

   call self%mpih%barrier(tictoc=.true.)
   call self%mpih%print_message('save HDF5 files t: '//trim(str(self%time%it,.true.))//', time: '//&
                                trim(str(self%time%time,.true.)))
   output_basename_ = trim(self%io%output_basename)//'-'//trim(strz(self%time%it,9))
   if (present(output_basename)) output_basename_ = trim(output_basename)
   call self%adam%io%save_xh5f(basename=trim(output_basename_), &
                               q=self%q, q_name=self%q_name,    &
                               with_ghost=with_ghost,           &
                               with_cell_morton=.true.,         &
                               t=self%time%it, time=self%time%time)
   call self%mpih%barrier(tictoc=.true.)
   endsubroutine save_xh5f

   subroutine save_residuals(self)
   !< Save residuals history.
   class(patch_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: v    !< Counter.

   if (self%time%is_to_save(it_save=self%io%residuals_save)) then
      call self%field%compute_normL2_residuals(dq=self%dq, norm=self%field%residuals)
      do v=1, self%nv
         call MPI_ALLREDUCE(MPI_IN_PLACE, self%field%residuals(v), 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, self%mpih%error)
         self%field%residuals(v) = sqrt(self%field%residuals(v))/sqrt(real(self%ni*self%nj*self%nk, R8P))
      enddo
      if (self%mpih%myrank==0) call self%io%save_residuals(it=self%time%it, time=self%time%time, &
                                                           blocks_number=self%blocks_number, residuals=self%field%residuals)
   endif
   endsubroutine save_residuals

   subroutine save_restart_files(self)
   !< Save restart files.
   class(patch_cpu_object), intent(inout) :: self !< The equation.

   call self%mpih%barrier(tictoc=.true.)
   call self%mpih%print_message('save restart files t: '//trim(str(self%time%it,.true.))//', time: '//&
                                trim(str(self%time%time,.true.)))
   call self%adam%save_restart_files(basename=self%io%restart_basename, t=self%time%it, time=self%time%time, q=self%q)
   call self%save_xh5f(output_basename=self%io%restart_basename)
   call self%mpih%barrier(tictoc=.true.)
   endsubroutine save_restart_files

   subroutine save_simulation_data(self)
   !< Save all simulation data.
   class(patch_cpu_object), intent(inout) :: self !< The equation.

   if ((self%time%is_to_save(it_save=self%io%it_save)).or.      &
       (self%time%is_to_save(it_save=self%io%restart_save))) then
      call self%update_ghost(q=self%q)
      if (self%time%is_to_save(it_save=self%io%it_save)) call self%save_xh5f
      if (mod(self%time%it,self%io%restart_save)==0) call self%save_restart_files
   endif
   endsubroutine save_simulation_data

   ! IC/BC
   subroutine set_boundary_conditions(self, q)
   !< Set boundary conditions of equation.
   class(patch_cpu_object), intent(in)    :: self                 !< The equation.
   real(R8P),               intent(inout) :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,1:)    !< Conservative variables.
   integer(I4P)                           :: b, c, i, j, k, v     !< Counter.
   integer(I4P)                           :: idelta,jdelta,kdelta !< IJK delta step for extrapolation.
   integer(I4P)                           :: bc_type              !< Boundary condition type.
   integer(I4P)                           :: crown                !< Crown counter.
   integer(I4P)                           :: fec                  !< Boundary fec (1 to 26).
   integer(I4P)                           :: fec_1_6              !< Boundary fec (1 to 6).

   associate(local_map_bc_crown=>self%field%maps%local_map_bc_crown, ngc=>self%ngc, q_bc_vars=>self%bc%q)
   if (allocated(self%field%maps%local_map_bc_crown)) then
      do crown=1, ngc
         do c=1, size(local_map_bc_crown, dim=1)
            b = local_map_bc_crown(c, 1 ,crown)
            if (b>0) then
               i       = local_map_bc_crown(c, 2 ,crown)
               j       = local_map_bc_crown(c, 3 ,crown)
               k       = local_map_bc_crown(c, 4 ,crown)
               bc_type = local_map_bc_crown(c, 8 ,crown)
               fec     = local_map_bc_crown(c, 9 ,crown)
               fec_1_6 = fec_1_6_array(fec)
               if (bc_type == BC_DIRICHLET) then
                  q(1,i,j,k,b) = q_bc_vars(1, fec_1_6)
               endif
            endif
         enddo
      enddo
   endif
   endassociate
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self)
   !< Set initial conditions and coils on field.
   class(patch_cpu_object), intent(inout) :: self !< The equation.

   call self%ic%set_initial_conditions(field=self%field, r=self%r)
   endsubroutine set_initial_conditions

   subroutine update_ghost(self, q, step)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(patch_cpu_object), intent(inout)        :: self            !< The equation.
   real(R8P),               intent(inout)        :: q(1:,         &
                                                      1-self%ngc:,&
                                                      1-self%ngc:,&
                                                      1-self%ngc:,&
                                                      1:)           !< Conservative variables.
   integer(I4P),            intent(in), optional :: step            !< Step to be perfordmed in asyncronous comp.
   logical                                       :: do_local_update !< Flag for triggering local update.
   logical                                       :: do_set_bc       !< Flag for triggering setting bc.

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

   if (do_local_update) call self%field%update_ghost_local(q=q)
                        call self%field%update_ghost_mpi(q=q, step=step)
   if (do_set_bc)       call self%set_boundary_conditions(q=q)
   endsubroutine update_ghost

   ! numerical methods
   subroutine integrate(self, compute_smoothing)
   !< Integrate the equation.
   class(patch_cpu_object), intent(inout) :: self              !< The equation.
   procedure(compute_smoothing_interface) :: compute_smoothing !< Smoothing procedure.
   real(R8P)                              :: timing(1:2)       !< Tic toc timing.
   real(R8P)                              :: timing_step(1:2)  !< Tic toc timing.
   real(R8P)                              :: tolerance         !< Maximum residual tolerance.
   real(R8P)                              :: dq_max            !< Maximum residual.
   integer(I4P)                           :: i                 !< Counter.
   integer(I4P)                           :: b                 !< Counter.
   real(R8P)                              :: dx2,dy2,dz2       !< Square space steps.

   tolerance = 0.000001_R8P
   call self%mpih%barrier(tictoc=.true., timing=timing(1), single=.true.)
   integration: do
      call self%mpih%barrier(tictoc=.true., timing=timing_step(1), single=.true.)
      self%time%it = self%time%it + 1

      if (self%io%save_memory_status) then
         call save_memory_status(file_name='memory_cpu-'//self%mpih%myrankstr//'.dat', tag=str(self%time%it,.true.))
      endif

      if (self%blocks_number>0) then
         do b=1, self%blocks_number
            dx2 = self%field%dxyz(1,b)*self%field%dxyz(1,b)
            dy2 = self%field%dxyz(2,b)*self%field%dxyz(2,b)
            dz2 = self%field%dxyz(3,b)*self%field%dxyz(3,b)
            call compute_smoothing(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, dx2=dx2, dy2=dy2, dz2=dz2, &
                                   r=self%r(1,:,:,:,b), dq=self%dq(1,:,:,:,b), q=self%q(1,:,:,:,b), dq_max=dq_max)
         enddo
      endif

      if (dq_max < tolerance) then
         call self%mpih%print_message('convergence reached at iteration '//trim(str(self%time%it,.true.)))
         exit
      endif
      self%dq = self%q

      call self%time%print_progress(nodes_number=self%adam%tree%nodes_number)

      call self%save_simulation_data

      if (((self%time%it_max <= 0).and.(self%time%time >= self%time%time_max)).or.&
         ((self%time%it>=self%time%it_max).and.(self%time%it_max > 0))) exit integration

      call self%mpih%barrier(tictoc=.true., timing=timing_step(2), single=.true.)
   enddo integration
   call self%mpih%barrier(tictoc=.true., timing=timing(2), single=.true.)
   endsubroutine integrate

   subroutine simulate(self, filename)
   !< Perform the simulation.
   class(patch_cpu_object), intent(inout) :: self     !< The equation.
   character(*),            intent(in)    :: filename !< Input file name.

   call self%initialize(filename=filename)

   select case(self%time%smoothing)
   case(SMOOTHING_MULTIGRID)
      call self%integrate(compute_smoothing=compute_smoothing_multigrid)
   case(SMOOTHING_GAUSS_SEIDEL)
      call self%integrate(compute_smoothing=compute_smoothing_gauss_seidel)
   case(SMOOTHING_SOR)
      call self%integrate(compute_smoothing=compute_smoothing_sor)
   case(SMOOTHING_SOR_OMP)
      call self%integrate(compute_smoothing=compute_smoothing_sor_omp)
   endselect

   call self%finalize
   endsubroutine simulate

   ! non TBP
   subroutine apply_bc_dirichlet(ni, nj, nk, ngc, q)
   integer(I4P), intent(in)    :: ni,nj,nk,ngc !< Grid dimensions.
   real(R8P),    intent(inout) :: q(1-ngc:,&
                                    1-ngc:,&
                                    1-ngc:)    !< Poisson field variable.
   q(0,:,:) = 0._R8P ; q(ni+1,:   ,:   ) = 0._R8P
   q(:,0,:) = 0._R8P ; q(:   ,nj+1,:   ) = 0._R8P
   q(:,:,0) = 0._R8P ; q(:   ,:   ,nk+1) = 0._R8P
   endsubroutine apply_bc_dirichlet

   subroutine compute_prolongation(ngc, nic, njc, nkc, coarse, fine)
   !< Compute trilinear prolongation.
   integer(I4P), intent(in)    :: ngc            !< Number of ghost cells.
   integer(I4P), intent(in)    :: nic,njc,nkc  !< Coarse grid dimensions.
   real(R8P),    intent(in)    :: coarse(1-ngc:,&
                                         1-ngc:,&
                                         1-ngc:) !< Coarse grid field variable.
   real(R8P),    intent(inout) :: fine(1-ngc:,&
                                       1-ngc:,&
                                       1-ngc:)   !< Fine grid field variable.
   integer(I4P)                :: i,j,k,b        !< Counter.
   integer(I4P)                :: ii,jj,kk       !< Counter.

   fine = 0._R8P
   do k = 1, nkc
   do j = 1, njc
   do i = 1, nic
      ii = 2*i - 1
      jj = 2*j - 1
      kk = 2*k - 1
      fine(ii,  jj,  kk  ) =              coarse(i,  j,  k  )
      fine(ii+1,jj,  kk  ) = 0.5_R8P   * (coarse(i,  j,  k  ) + coarse(i+1,j,  k  ))
      fine(ii,  jj+1,kk  ) = 0.5_R8P   * (coarse(i,  j,  k  ) + coarse(i,  j+1,k  ))
      fine(ii,  jj,  kk+1) = 0.5_R8P   * (coarse(i,  j,  k  ) + coarse(i,  j,  k+1))
      fine(ii+1,jj+1,kk  ) = 0.25_R8P  * (coarse(i,  j,  k  ) + coarse(i+1,j,  k  ) + &
                                          coarse(i,  j+1,k  ) + coarse(i+1,j+1,k  ))
      fine(ii+1,jj,  kk+1) = 0.25_R8P  * (coarse(i,  j,  k  ) + coarse(i+1,j,  k  ) + &
                                          coarse(i,  j,  k+1) + coarse(i+1,j,  k+1))
      fine(ii,  jj+1,kk+1) = 0.25_R8P  * (coarse(i,  j,  k  ) + coarse(i,  j+1,k  ) + &
                                          coarse(i,  j,  k+1) + coarse(i,  j+1,k+1))
      fine(ii+1,jj+1,kk+1) = 0.125_R8P * (coarse(i,  j,  k  ) + coarse(i+1,j,  k  ) + &
                                          coarse(i,  j+1,k  ) + coarse(i,  j,  k+1) + &
                                          coarse(i+1,j+1,k  ) + coarse(i+1,j,  k+1) + &
                                          coarse(i,  j+1,k+1) + coarse(i+1,j+1,k+1))
   enddo
   enddo
   enddo
   endsubroutine compute_prolongation

   subroutine compute_residuals(ni, nj, nk, ngc, dx2, dy2, dz2, r, q, dq)
   !< Compute residuals of `r-laplacian(q)`.
   integer(I4P), intent(in)    :: ni,nj,nk,ngc !< Grid dimensions.
   real(R8P),    intent(in)    :: dx2,dy2,dz2  !< Square space steps.
   real(R8P),    intent(in)    :: r(1-ngc:,&
                                    1-ngc:,&
                                    1-ngc:)    !< Rho distribution.
   real(R8P),    intent(in)    :: q(1-ngc:,&
                                    1-ngc:,&
                                    1-ngc:)    !< Poisson field variable.
   real(R8P),    intent(inout) :: dq(1-ngc:,&
                                     1-ngc:,&
                                     1-ngc:)   !< Residuals.
   real(R8P)                   :: laplacian    !< Laplacian value.
   integer(I4P)                :: i,j,k,b      !< Counter.

   do k = 1, nk
   do j = 1, nj
   do i = 1, ni
      laplacian = (q(i+1,j,  k  ) - 2._R8P * q(i,j,k) + q(i-1,j,  k  )) / dx2 + &
                  (q(i,  j+1,k  ) - 2._R8P * q(i,j,k) + q(i,  j-1,k  )) / dy2 + &
                  (q(i,  j,  k+1) - 2._R8P * q(i,j,k) + q(i,  j,  k-1)) / dz2
      dq(i,j,k) = r(i,j,k) - laplacian
   enddo
   enddo
   enddo
   call apply_bc_dirichlet(ni=ni, nj=nj, nk=nk, ngc=ngc, q=dq)
   endsubroutine compute_residuals

   subroutine compute_restriction(ngc, nic, njc, nkc, fine, coarse)
   !< Compute full weighting restriction.
   integer(I4P), intent(in)    :: ngc            !< Number of ghost cells.
   integer(I4P), intent(in)    :: nic,njc,nkc    !< Coarse grid dimensions.
   real(R8P),    intent(in)    :: fine(1-ngc:,&
                                       1-ngc:,&
                                       1-ngc:)   !< Fine grid field variable.
   real(R8P),    intent(inout) :: coarse(1-ngc:,&
                                         1-ngc:,&
                                         1-ngc:) !< Coarse grid field variable.
   integer(I4P)                :: i,j,k,b        !< Counter.
   integer(I4P)                :: ii,jj,kk       !< Counter.

   do k = 1, nkc
   do j = 1, njc
   do i = 1, nic
      ii = 2*i - 1
      jj = 2*j - 1
      kk = 2*k - 1
      coarse(i,j,k) = 0.125_R8P * (fine(ii,  jj,  kk  ) + fine(ii+1,jj,  kk  ) + &
                                   fine(ii,  jj+1,kk  ) + fine(ii,  jj,  kk+1) + &
                                   fine(ii+1,jj+1,kk  ) + fine(ii+1,jj,  kk+1) + &
                                   fine(ii,  jj+1,kk+1) + fine(ii+1,jj+1,kk+1))
   enddo
   enddo
   enddo
   call apply_bc_dirichlet(ni=nic, nj=njc, nk=nkc, ngc=ngc, q=coarse)
   endsubroutine compute_restriction

   subroutine compute_smoothing_gauss_seidel(ni, nj, nk, ngc, dx2, dy2, dz2, r, dq, q, dq_max, iterations)
   !< Compute smoothing by Gauss-Seidel method.
   integer(I4P), intent(in)              :: ni,nj,nk,ngc !< Grid dimensions.
   real(R8P),    intent(in)              :: dx2,dy2,dz2  !< Square space steps.
   real(R8P),    intent(in)              :: r(1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:)    !< Rho distribution.
   real(R8P),    intent(inout), optional :: dq(1-ngc:,&
                                               1-ngc:,&
                                               1-ngc:)   !< Residuals.
   real(R8P),    intent(inout)           :: q(1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:)    !< Poisson field variable.
   real(R8P),    intent(inout), optional :: dq_max       !< Maximum residual.
   integer(I4P), intent(in),    optional :: iterations   !< Smoothing iterations.
   integer(I4P)                          :: iterations_  !< Smoothing iterations, local var.
   real(R8P)                             :: dq_max_      !< Maximum residual, local var.
   real(R8P)                             :: q_old        !< Previous q.
   real(R8P)                             :: factor       !< Jacobi relaxation factor.
   integer(I4P)                          :: i,j,k,b,iter !< Counter.

   iterations_ = 1 ; if (present(iterations)) iterations_ = iterations
   factor = 1._R8P / (2._R8P * (1._R8P/dx2 + 1._R8P/dy2 + 1._R8P/dz2))
   dq_max_ = 0._R8P
   do iter=1, iterations_
      call apply_bc_dirichlet(ni=ni, nj=nj, nk=nk, ngc=ngc, q=q)
      do k = 1, nk
      do j = 1, nj
      do i = 1, ni
         q_old = q(i,j,k)
         q(i,j,k) = factor * ((q(i+1,j,  k  ) + q(i-1,j,  k  )) / dx2 + &
                              (q(i,  j+1,k  ) + q(i,  j-1,k  )) / dy2 + &
                              (q(i,  j,  k+1) + q(i,  j,  k-1)) / dz2 - &
                               r(i,  j,  k  ))
         dq_max_ = max(dq_max_, abs(q(i,j,k) - q_old))
      enddo
      enddo
      enddo
   enddo
   if (present(dq_max)) dq_max = dq_max_
   endsubroutine compute_smoothing_gauss_seidel

   subroutine compute_smoothing_multigrid(ni, nj, nk, ngc, dx2, dy2, dz2, r, dq, q, dq_max, iterations)
   integer(I4P), intent(in)              :: ni,nj,nk,ngc         !< Grid dimensions.
   real(R8P),    intent(in)              :: dx2,dy2,dz2          !< Square space steps.
   real(R8P),    intent(in)              :: r(1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:)            !< Rho distribution.
   real(R8P),    intent(inout), optional :: dq(1-ngc:,&
                                               1-ngc:,&
                                               1-ngc:)           !< Residuals.
   real(R8P),    intent(inout)           :: q(1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:)            !< Poisson field variable.
   real(R8P),    intent(inout), optional :: dq_max               !< Maximum residual.
   integer(I4P), intent(in),    optional :: iterations           !< Smoothing iterations.
   integer(I4P)                          :: iterations_          !< Smoothing iterations, local var.
   real(R8P)                             :: r_c(1-ngc:ni/2+ngc,&
                                                1-ngc:nj/2+ngc,&
                                                1-ngc:nk/2+ngc)  !< Rho distribution, coarse grid.
   real(R8P)                             :: dq_c(1-ngc:ni/2+ngc,&
                                                 1-ngc:nj/2+ngc,&
                                                 1-ngc:nk/2+ngc) !< Residuals, coarse grid.
   real(R8P)                             :: q_c(1-ngc:ni/2+ngc,&
                                                1-ngc:nj/2+ngc,&
                                                1-ngc:nk/2+ngc)  !< Poisson field variable, coarse grid.
   integer(I4P)                          :: iter                 !< Counter.

   iterations_ = 1 ; if (present(iterations)) iterations_ = iterations
   ! V-cicle
   do iter=1, iterations_
      call compute_smoothing_gauss_seidel(ni=ni,nj=nj,nk=nk,ngc=ngc,dx2=dx2,dy2=dy2,dz2=dz2,r=r,q=q,iterations=3)
      call compute_residuals(ni=ni,nj=nj,nk=nk,ngc=ngc,dx2=dx2,dy2=dy2,dz2=dz2,r=r,q=q,dq=dq)
      call compute_restriction(ngc=ngc, nic=ni/2, njc=nj/2, nkc=nk/2, fine=dq, coarse=r_c)
      q_c = 0._R8P
      call compute_smoothing_gauss_seidel(ni=ni/2,nj=nj/2,nk=nk/2,ngc=ngc,dx2=dx2*4,dy2=dy2*4,dz2=dz2*4,r=r_c,q=q_c,iterations=10)
      call compute_prolongation(ngc=ngc, nic=ni/2, njc=nj/2, nkc=nk/2, coarse=q_c, fine=dq)
      q = q + dq
      call compute_smoothing_gauss_seidel(ni=ni,nj=nj,nk=nk,ngc=ngc,dx2=dx2,dy2=dy2,dz2=dz2,r=r,q=q,dq_max=dq_max,iterations=3)
   enddo
   endsubroutine compute_smoothing_multigrid

   subroutine compute_smoothing_sor(ni, nj, nk, ngc, dx2, dy2, dz2, r, dq, q, dq_max, iterations)
   !< Compute smoothing by SOR, Successive Overrelaxation.
   integer(I4P), intent(in)              :: ni,nj,nk,ngc !< Grid dimensions.
   real(R8P),    intent(in)              :: dx2,dy2,dz2  !< Square space steps.
   real(R8P),    intent(in)              :: r(1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:)    !< Rho distribution.
   real(R8P),    intent(inout), optional :: dq(1-ngc:,&
                                               1-ngc:,&
                                               1-ngc:)   !< Residuals.
   real(R8P),    intent(inout)           :: q(1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:)    !< Poisson field variable.
   real(R8P),    intent(inout), optional :: dq_max       !< Maximum residual.
   integer(I4P), intent(in),    optional :: iterations   !< Smoothing iterations.
   integer(I4P)                          :: iterations_  !< Smoothing iterations, local var.
   real(R8P)                             :: dq_max_      !< Maximum residual, local var.
   real(R8P)                             :: factor       !< Jacobi relaxation factor.
   real(R8P)                             :: omega        !< Overrelaxation parameter.
   real(R8P)                             :: q_old        !< Previous q.
   integer(I4P)                          :: i,j,k,b,iter !< Counter.

   iterations_ = 1 ; if (present(iterations)) iterations_ = iterations
   omega = 1.8_R8P
   factor = 1._R8P / (2._R8P * (1._R8P/dx2 + 1._R8P/dy2 + 1._R8P/dz2))
   dq_max_ = 0._R8P
   do iter=1, iterations_
      call apply_bc_dirichlet(ni=ni, nj=nj, nk=nk, ngc=ngc, q=q)
      do k = 1, nk
      do j = 1, nj
      do i = 1, ni
         q_old = q(i,j,k)
         q(i,j,k) = factor * ((q(i+1,j,  k  ) +  q(i-1,j,  k  )) / dx2 + &
                              (q(i,  j+1,k  ) +  q(i,  j-1,k  )) / dy2 + &
                              (q(i,  j,  k+1) +  q(i,  j,  k-1)) / dz2 - &
                               r(i,  j,  k  ))
         q(i,j,k) = q_old + omega * (q(i,j,k) - q_old)
         dq_max_ = max(dq_max, abs(q(i,j,k) - q_old))
      enddo
      enddo
      enddo
   enddo
   if (present(dq_max)) dq_max = dq_max_
   endsubroutine compute_smoothing_sor

   subroutine compute_smoothing_sor_omp(ni, nj, nk, ngc, dx2, dy2, dz2, r, dq, q, dq_max, iterations)
   !< Compute smoothing by SOR, Successive Overrelaxation (parallel OpenMP) method: red-black 3D scheme.
   integer(I4P), intent(in)              :: ni,nj,nk,ngc !< Grid dimensions.
   real(R8P),    intent(in)              :: dx2,dy2,dz2  !< Square space steps.
   real(R8P),    intent(in)              :: r(1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:)    !< Rho distribution.
   real(R8P),    intent(inout), optional :: dq(1-ngc:,&
                                               1-ngc:,&
                                               1-ngc:)   !< Residuals.
   real(R8P),    intent(inout)           :: q(1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:)    !< Poisson field variable.
   real(R8P),    intent(inout), optional :: dq_max       !< Maximum residual.
   integer(I4P), intent(in),    optional :: iterations   !< Smoothing iterations.
   integer(I4P)                          :: iterations_  !< Smoothing iterations, local var.
   real(R8P)                             :: dq_max_      !< Maximum residual, local var.
   real(R8P)                             :: factor       !< Jacobi relaxation factor.
   real(R8P)                             :: omega        !< Overrelaxation parameter.
   integer(I4P)                          :: i,j,k,b,iter !< Counter.
   integer(I4P)                          :: color        !< Color counter.
   real(R8P)                             :: residual     !< Residual.

   iterations_ = 1 ; if (present(iterations)) iterations_ = iterations
   omega = 1.8_R8P
   factor = 1._R8P / (2._R8P * (1._R8P/dx2 + 1._R8P/dy2 + 1._R8P/dz2))
   dq_max_ = 0._R8P
   do iter=1, iterations_
      call apply_bc_dirichlet(ni=ni, nj=nj, nk=nk, ngc=ngc, q=q)
      do color=0, 7
         !$omp parallel do firstprivate(b,factor,omega,color) private(residual) shared(q,dq) reduction(max:dq_max) collapse(3)
         do k = 1, nk
         do j = 1, nj
         do i = 1, ni
            if (mod(i+j+k,8)==color) then
               residual = factor * ((q(i+1,j,  k  ) +  q(i-1,j,  k  )) / dx2 + &
                                    (q(i,  j+1,k  ) +  q(i,  j-1,k  )) / dy2 + &
                                    (q(i,  j,  k+1) +  q(i,  j,  k-1)) / dz2 - &
                                     r(i,  j,  k  ))
               dq(i,j,k) = q(i,j,k) + omega * (residual - q(i,j,k))
               dq_max_ = max(dq_max_, abs(residual - q(i,j,k)))
            endif
         enddo
         enddo
         enddo
         !$omp parallel do firstprivate(b,color) shared(q,dq) collapse(3)
         do k = 1, nk
         do j = 1, nj
         do i = 1, ni
            if (mod(i+j+k,8)==color) q(i,j,k) = dq(i,j,k)
         enddo
         enddo
         enddo
      enddo
   enddo
   if (present(dq_max)) dq_max = dq_max_
   endsubroutine compute_smoothing_sor_omp
endmodule adam_patch_cpu_object
