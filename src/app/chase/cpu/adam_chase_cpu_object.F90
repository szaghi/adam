!< ADAM, CHASE (CFD-HPC enabled-Adaptive mesh-Simulation code for-Euler equations) class definition, CPU backend.
module adam_chase_cpu_object
!< ADAM, CHASE (CFD-HPC enabled-Adaptive mesh-Simulation code for-Euler equations) class definition, CPU backend.

! ADAM modules
use adam_common_library
! chase modules
use adam_chase_common_library
! third party modules
use penf
use mpi

implicit none
private
public :: chase_cpu_object

type, extends(chase_common_object) :: chase_cpu_object !commentate procedure AMR e IB
   !< Maxwell equations system class definition, CPU backend.
   real(R8P), allocatable :: dq( :,:,:,:,:) !< Residuals right hand side.
   real(R8P), allocatable :: flx(:,:,:,:,:) !< Fluxes along x.
   real(R8P), allocatable :: fly(:,:,:,:,:) !< Fluxes along y.
   real(R8P), allocatable :: flz(:,:,:,:,:) !< Fluxes along z.
   contains
      ! auxiliary methods
      procedure, pass(self) :: allocate_cpu !< Allocate CPU data.
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
      procedure, pass(self) :: compute_dt          !< Compute time step.
      procedure, pass(self) :: compute_q_auxiliary !< Compute auxiliary variables.
      procedure, pass(self) :: compute_residuals   !< Compute residuals.
      procedure, pass(self) :: integrate           !< Perform one step integration.
      procedure, pass(self) :: simulate            !< Perform the simulation.
endtype chase_cpu_object

interface assign_omp
!< Assign array to scalar value with OpenMP threads.
module procedure assign_omp_R8P_5D
endinterface assign_omp

contains
   ! auxiliary methods
   subroutine allocate_cpu(self)
   !< Allocate CPU data.
   class(chase_cpu_object), intent(inout) :: self !< The equation.
   character(:), allocatable              :: msg_ !< Allocating message base.
   character(:), allocatable              :: msg  !< Allocating message.

   call self%mpih%print_message('chase_cpu_object%allocate_cpu start')
   msg_ = self%mpih%myrankstr//'chase_cpu_object%allocate_cpu '
   associate(nv=>self%nv, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
              nb=>self%nb, weno_s=>self%weno%S, solids_number=>self%ib%solids_number)
   msg = msg_//' dq '
   call allocate_variable(var=self%dq, ulb=reshape([1,nv,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,5]),msg=msg)
   self%dq = 0._R8P
   msg = msg_//' flx '
   call allocate_variable(var=self%flx,ulb=reshape([1,nv,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,5]),msg=msg)
   self%flx = 0._R8P
   msg = msg_//' fly '
   call allocate_variable(var=self%fly,ulb=reshape([1,nv,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,5]),msg=msg)
   self%fly = 0._R8P
   msg = msg_//' flz '
   call allocate_variable(var=self%flz,ulb=reshape([1,nv,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,5]),msg=msg)
   self%flz = 0._R8P
   endassociate
   call self%mpih%print_message('chase_cpu_object%allocate_cpu finish')
   endsubroutine allocate_cpu

   subroutine initialize(self, filename)
   !< Initialize the equation.
   class(chase_cpu_object), intent(inout) :: self     !< The equation.
   character(*),            intent(in)    :: filename !< Input file name.

   call self%mpih%initialize(do_mpi_init=.true., verbose=.true.)
   call self%mpih%print_message('chase_cpu_object%initialize start')
   call self%initialize_common(field = self%adam%field, filename=filename, memory_avail=self%mpih%memory_avail)
   call self%allocate_cpu
   print '(A)', self%mpih%description()
   call self%mpih%print_message('chase_cpu_object%initialize finish')
   endsubroutine initialize

   ! AMR methods
   subroutine amr_update(self)
   !< Do AMR update.
   class(chase_cpu_object), intent(inout) :: self                !< The equation.
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
   class(chase_cpu_object), intent(inout) :: self !< The equation.

   if (self%ib%solids_number>0) then
      call self%mpih%print_message('compute IB distance start')
      call self%ib%compute_phi
      call self%mpih%print_message('compute IB distance finish')
   endif
   endsubroutine compute_phi

   subroutine mark_by_geo(self, delta_fine, delta_coarse, threshold, do_init)
   !< Mark blocks to be refined/derefined by a geometric constrain.
   class(chase_cpu_object), intent(inout)        :: self           !< The equation.
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
   class(chase_cpu_object), intent(inout)        :: self                     !< The equation.
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

   ivar_     = 1_R4P    ; if (present(ivar)) ivar_ = ivar
   do_init_ = .true.    ; if (present(do_init)) do_init_ = do_init
   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   if (do_init_) self%field%refinements_needed = [(TO_BE_DEREFINED,b=1,self%blocks_number)]
   associate (ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
              blocks_number=>self%blocks_number, dxyz=>self%field%dxyz)
      select case(delta_type)
      case(AMR_DELTA_T_X)
         dc(1:blocks_number) = dxyz(1,1:blocks_number)
      case(AMR_DELTA_T_Y)
         dc(1:blocks_number) = dxyz(2,1:blocks_number)
      case(AMR_DELTA_T_Z)
         dc(1:blocks_number) = dxyz(3,1:blocks_number)
      case(AMR_DELTA_T_MAX)
         do b=1, blocks_number
            dc(b) = maxval(dxyz(:,b))
         enddo
      endselect
      call self%update_ghost(q=self%q)
      call self%compute_q_auxiliary(q=self%q, q_aux=self%q_aux)
      do b=1, blocks_number
         call compute_q_gradient(b=b, ni=ni, nj=nj, nk=nk, ngc=ngc, &
                                 dx=dxyz(1,b), dy=dxyz(2,b), dz=dxyz(3,b), q=self%q_aux, ivar=ivar_, gradient=grad_var)
         max_cell_delta = max_cell_delta_grad(grad=grad_var)
         if ((dc(b)) > max_cell_delta) then
            self%field%refinements_needed(b) = TO_BE_REFINED
         elseif ((dc(b)) * threshold_ < max_cell_delta) then
            self%field%refinements_needed(b) = max(self%field%refinements_needed(b), TO_BE_DEREFINED)
         else
            self%field%refinements_needed(b) = max(self%field%refinements_needed(b), TO_NOT_TOUCH)
         endif
      enddo
   endassociate
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
   class(chase_cpu_object), intent(inout) :: self        !< The equation.
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
   class(chase_cpu_object), intent(inout) :: self              !< The equation.
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
   class(chase_cpu_object), intent(inout) :: self      !< The equation.
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
   class(chase_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P),            intent(out)   :: t    !< Time iteration.
   real(R8P),               intent(out)   :: time !< Time.

   call self%adam%load_restart_files(basename=self%io%restart_basename, t=t, time=time, q=self%q)
   call self%adam%make_comm_local_maps_ghost_bc
   endsubroutine load_restart_files

   subroutine save_xh5f(self, output_basename, with_ghost)
   !< Save simulation data in HDF5 format.
   class(chase_cpu_object), intent(inout)        :: self             !< The equation.
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
   class(chase_cpu_object), intent(inout) :: self !< The equation.
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
   class(chase_cpu_object), intent(inout) :: self !< The equation.

   call self%mpih%barrier(tictoc=.true.)
   call self%mpih%print_message('save restart files t: '//trim(str(self%time%it,.true.))//', time: '//&
                                trim(str(self%time%time,.true.)))
   call self%adam%save_restart_files(basename=self%io%restart_basename, t=self%time%it, time=self%time%time, q=self%q)
   call self%save_xh5f(output_basename=self%io%restart_basename)
   call self%mpih%barrier(tictoc=.true.)
   endsubroutine save_restart_files

   subroutine save_simulation_data(self)
   !< Save all simulation data.
   class(chase_cpu_object), intent(inout) :: self !< The equation.

   if ((self%time%is_to_save(it_save=self%io%it_save)).or.      &
       (self%time%is_to_save(it_save=self%io%restart_save)).or. &
       (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max))) then
      call self%update_ghost(q=self%q)
      call self%compute_q_auxiliary(q=self%q, q_aux=self%q_aux)
      if (self%time%is_to_save(it_save=self%io%it_save)) call self%save_xh5f
      if (mod(self%time%it,self%io%restart_save)==0) call self%save_restart_files
      if (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max))&
         call self%slices%save_mat(basename=self%io%output_basename, &
                                   it=self%time%it,                  &
                                   it_max=self%time%it_max,          &
                                   time=self%time%time,              &
                                   time_max=self%time%time_max,      &
                                   adam=self%adam,                   &
                                   q=self%q,                         &
                                   q_name=['rho','rhu','rhv','rhw','rhe'])
   endif
   endsubroutine save_simulation_data

   ! IC/BC
   subroutine set_boundary_conditions(self, q)
   !< Set boundary conditions of equation.
   class(chase_cpu_object), intent(in)    :: self                 !< The equation.
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

   associate(local_map_bc_crown=>maps%local_map_bc_crown, &
             nv=>self%nv, ngc=>self%ngc, cv=>self%physics%eos%cv, R=>self%physics%eos%R, q_bc_vars=>self%bc%q)
   if (allocated(maps%local_map_bc_crown)) then
      do crown=1, ngc
         do c=1, size(local_map_bc_crown, dim=1)
            b = local_map_bc_crown(c, 1 ,crown)
            if (b>0) then
               i       = local_map_bc_crown(c, 2 ,crown)
               j       = local_map_bc_crown(c, 3 ,crown)
               k       = local_map_bc_crown(c, 4 ,crown)
               idelta  = local_map_bc_crown(c, 5 ,crown)
               jdelta  = local_map_bc_crown(c, 6 ,crown)
               kdelta  = local_map_bc_crown(c, 7 ,crown)
               bc_type = local_map_bc_crown(c, 8 ,crown)
               fec     = local_map_bc_crown(c, 9 ,crown)
               fec_1_6 = fec_1_6_array(fec)
               if (bc_type == BC_EXTRAPOLATION) then
                  do v=1, nv
                     q(v,i,j,k,b) = q(v,i-idelta,j-jdelta,k-kdelta,b)
                  enddo
               elseif (bc_type == BC_INFLOW) then
                   q(1,i,j,k,b) = q_bc_vars(1, fec_1_6)
                   q(2,i,j,k,b) = q_bc_vars(1, fec_1_6)* q_bc_vars(2, fec_1_6)
                   q(3,i,j,k,b) = q_bc_vars(1, fec_1_6)* q_bc_vars(3, fec_1_6)
                   q(4,i,j,k,b) = q_bc_vars(1, fec_1_6)* q_bc_vars(4, fec_1_6)
                   q(5,i,j,k,b) = q_bc_vars(1, fec_1_6)*                                &
                                  (cv*q_bc_vars(5, fec_1_6)/(q_bc_vars(1, fec_1_6)*R) + &
                                  0.5_R8P*(q_bc_vars(2, fec_1_6)**2+q_bc_vars(3, fec_1_6)**2+q_bc_vars(4, fec_1_6)**2))
               endif
            endif
         enddo
      enddo
   endif
   endassociate
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self)
   !< Set initial conditions and coils on field.
   class(chase_cpu_object), intent(inout) :: self !< The equation.

   call self%ic%set_initial_conditions(physics=self%physics, field=self%field, q=self%q)
   endsubroutine set_initial_conditions

   subroutine update_ghost(self, q, step)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(chase_cpu_object), intent(inout)        :: self            !< The equation.
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
   subroutine compute_dt(self)
   !< Compute maximum time step accordingly to CFL stabilty criterion.
   class(chase_cpu_object), intent(inout) :: self                            !< The equation.
   real(R8P)                              :: umax                            !< Maximum speed of waves propagation.
   real(R8P)                              :: ss                              !< Speed of sound.
   real(R8P)                              :: dx_locale, dy_locale, dz_locale !< Local space steps.
   integer(I4P)                           :: b, i, j, k                      !< Counter.

   call self%compute_q_auxiliary(q=self%q, q_aux=self%q_aux)
   self%time%dt = huge(1._R8P)
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, blocks_number=>self%blocks_number, mu=>self%physics%eos%mu, &
             dx=>self%field%dxyz(1,:), dy=>self%field%dxyz(2,:), dz=>self%field%dxyz(3,:), q_aux=>self%q_aux)
   umax = 0._R8P
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               dx_locale = dx(b)*0.5_R8P
               dy_locale = dy(b)*0.5_R8P
               dz_locale = dz(b)*0.5_R8P
               ss = q_aux(VAR_A,i,j,k,b)
               umax = max(umax, (abs(q_aux(VAR_U,i,j,k,b)) + ss)/dx_locale + &
                                (abs(q_aux(VAR_V,i,j,k,b)) + ss)/dy_locale + &
                                (abs(q_aux(VAR_W,i,j,k,b)) + ss)/dz_locale)
            enddo
         enddo
      enddo
   enddo
   endassociate
   self%time%dt = min(self%time%dt, self%time%CFL / umax)
   call MPI_ALLREDUCE(MPI_IN_PLACE, self%time%dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, self%mpih%error)
   endsubroutine compute_dt

   subroutine compute_q_auxiliary(self, q, q_aux)
   !< Compute auxiliary variables.
   class(chase_cpu_object), intent(in)  :: self       !< The equation.
   real(R8P),               intent(in)  :: q(1:,         &
                                             1-self%ngc:,&
                                             1-self%ngc:,&
                                             1-self%ngc:,&
                                             1:)      !< Conservative variables.
   real(R8P),               intent(out) :: q_aux(1:,         &
                                                 1-self%ngc:,&
                                                 1-self%ngc:,&
                                                 1-self%ngc:,&
                                                 1:)  !< Auxiliary variables.
   integer(I4P)                         :: b, i, j, k !< Counter.

   associate(eos=>self%physics%eos,ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, blocks_number=>self%blocks_number)
   do b=1, blocks_number
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
      call compute_auxiliary(cv=eos%cv, Rf=eos%R, g=eos%g, conservative=q(:,i,j,k,b), auxiliary=q_aux(:,i,j,k,b))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_q_auxiliary

   subroutine compute_residuals(self, q, q_aux, dq)
   !< Compute residuals of equation.
   class(chase_cpu_object), intent(inout) :: self      !< The equation.
   real(R8P),               intent(inout) :: q(1:,       &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1:)     !< Conservative variables.
   real(R8P),               intent(inout) :: q_aux(1:,       &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1:) !< Auxiliary variables.
   real(R8P),               intent(inout) :: dq(1:,         &
                                                1-self%ngc:,&
                                                1-self%ngc:,&
                                                1-self%ngc:,&
                                                1:)    !< Residuals.

   call self%update_ghost(q=q)
   call self%integrate_eikonal(q=q)
   call self%compute_q_auxiliary(q=q, q_aux=q_aux)
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv, blocks_number=>self%blocks_number, &
             dx=>self%field%dxyz(1,:), dy=>self%field%dxyz(2,:), dz=>self%field%dxyz(3,:),                         &
             phi=>self%ib%phi, flx=>self%flx, fly=>self%fly, flz=>self%flz,                                        &
             weno_s=>self%weno%S,                                                                                  &
             weno_a=>self%weno%a, weno_p=>self%weno%p, weno_d=>self%weno%d, ror_number=>self%weno%ror_number,      &
             ror_schemes=>self%weno%ror_schemes, ror_ivar=>self%weno%ror_ivar,                                     &
             ror_threshold=>self%weno%ror_threshold, enable_ror_stats=>self%weno%enable_ror_stats,                 &
             cell_scheme=>self%weno%cell_scheme, ror_stats=>self%weno%ror_stats, weno_zeps=>self%weno%zeps,        &
             solids_number=>self%ib%solids_number,                                                                 &
             g=>self%physics%eos%g, null_xyz=>self%grid%null_xyz)
   if (blocks_number > 0) then
      if (.not.null_xyz(1)) then
         call compute_fluxes_convective(dir=1,blocks_number=blocks_number,ni=ni,nj=nj,nk=nk,ngc=ngc,nv=nv,g=g,      &
                                        weno_s=weno_S,weno_zeps=weno_zeps,weno_a=weno_a,weno_p=weno_p,weno_d=weno_d,&
                                        q_aux=q_aux,fluxes=flx)
      else
         call assign_omp(blocks_number=blocks_number, ngc=ngc, lhs=flx, rhs=0._R8P)
      endif
      if (.not.null_xyz(2)) then
         call compute_fluxes_convective(dir=2,blocks_number=blocks_number,ni=ni,nj=nj,nk=nk,ngc=ngc,nv=nv,g=g,      &
                                        weno_s=weno_S,weno_zeps=weno_zeps,weno_a=weno_a,weno_p=weno_p,weno_d=weno_d,&
                                        q_aux=q_aux,fluxes=fly)
      else
         call assign_omp(blocks_number=blocks_number, ngc=ngc, lhs=fly, rhs=0._R8P)
      endif
      if (.not.null_xyz(3)) then
         call compute_fluxes_convective(dir=3,blocks_number=blocks_number,ni=ni,nj=nj,nk=nk,ngc=ngc,nv=nv,g=g,      &
                                        weno_s=weno_S,weno_zeps=weno_zeps,weno_a=weno_a,weno_p=weno_p,weno_d=weno_d,&
                                        q_aux=q_aux,fluxes=flz)
      else
         call assign_omp(blocks_number=blocks_number, ngc=ngc, lhs=flz, rhs=0._R8P)
      endif
      if (solids_number>0) then
         call compute_fluxes_difference(null_xyz=null_xyz,                                                                   &
                                        blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, ib_eps=1.e-12_R8P, &
                                        dx=dx, dy=dy, dz=dz, flx=flx, fly=fly, flz=flz, phi=phi, dq=dq)
      else
         call compute_fluxes_difference(null_xyz=null_xyz,                                                                   &
                                        blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, ib_eps=1.e-12_R8P, &
                                        dx=dx, dy=dy, dz=dz, flx=flx, fly=fly, flz=flz, dq=dq)
      endif
   endif
   endassociate
   endsubroutine compute_residuals

   subroutine integrate(self, do_ghost_syncro)
   !< Perform one step integration.
   class(chase_cpu_object), intent(inout)         :: self             !< The equation.
   logical,                 intent(in),  optional :: do_ghost_syncro  !< Flag to do syncrous ghost update.
   logical                                        :: do_ghost_syncro_ !< Flag to do syncrous ghost update, local var.
   integer(I4P)                                   :: s                !< Counter.

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   call self%rk%initialize_stages(q=self%q)
   select case(self%rk%scheme)
   case(RK_1, RK_2, RK_3)
      ! low storage RK working on q_rk(:,:,:,:,:,1)/q as stages, update q in place
      do s=1, self%rk%nrk
         call self%compute_residuals(q=self%q, q_aux=self%q_aux, dq=self%dq)
         if (s==1) call self%save_residuals
         if (self%ib%solids_number>0) then
            call self%rk%compute_stage_ls(s=s,dt=self%time%dt,phi=self%ib%phi,dq=self%dq,q=self%q)
         else
            call self%rk%compute_stage_ls(s=s,dt=self%time%dt,dq=self%dq,q=self%q)
         endif
      enddo
   case(RK_SSP_22, RK_SSP_33, RK_SSP_54)
      ! RK working on q_rk as stages
      do s=1, self%rk%nrk
         if (self%ib%solids_number>0) then
            call self%rk%compute_stage(s=s, dt=self%time%dt, phi=self%ib%phi)
         else
            call self%rk%compute_stage(s=s, dt=self%time%dt)
         endif
         call self%compute_residuals(q=self%rk%q_rk(:,:,:,:,:,s), q_aux=self%q_aux, dq=self%dq)
         if (s==1) call self%save_residuals
         if (self%ib%solids_number>0) then
            call self%rk%assign_stage(s=s, q=self%dq, phi=self%ib%phi)
         else
            call self%rk%assign_stage(s=s, q=self%dq)
         endif
      enddo
      if (self%ib%solids_number>0) then
         call self%rk%update_q(dt=self%time%dt, phi=self%ib%phi, q=self%q)
      else
         call self%rk%update_q(dt=self%time%dt, q=self%q)
      endif
   endselect
   endsubroutine integrate

   subroutine simulate(self, filename)
   !< Perform the simulation.
   class(chase_cpu_object), intent(inout) :: self             !< The equation.
   character(*),            intent(in)    :: filename         !< Input file name.
   real(R8P)                              :: timing(1:2)      !< Tic toc timing.
   real(R8P)                              :: timing_step(1:2) !< Tic toc timing.
   integer(I4P)                           :: i                !< Counter.

   ! initialization
   call self%initialize(filename=filename)
   if (self%io%restart) then
      call self%mpih%print_message('restart simulation from "'//trim(self%io%restart_basename)//'" files')
      call self%load_restart_files(t=self%time%it, time=self%time%time)
      call self%mpih%print_message('restart [t, time]: '//trim(str(self%time%it))//', '//trim(str(self%time%time)))
   else
      call self%mpih%print_message('impose initial conditions start')
      do i=1, self%ic%amr_iterations
         call self%mpih%print_message('  AMR/set IC iteration:'//trim(str(i,.true.)))
         call self%set_initial_conditions
         if (self%ib%solids_number > 0) call self%compute_phi()
         call self%amr_update
      enddo
      call self%set_initial_conditions
      self%time%time = 0._R8P
      self%time%it = 0
      call self%mpih%print_message('impose initial conditions finish')
   endif
   if (self%ib%solids_number > 0) call self%compute_phi()
   ! call self%amr_update
   call self%save_simulation_data
   if (self%mpih%myrank==0) call self%io%open_file_residuals(nv=self%nv)

   ! integration
   call self%mpih%barrier(tictoc=.true., timing=timing(1), single=.true.)
   integration: do
      call self%mpih%barrier(tictoc=.true., timing=timing_step(1), single=.true.)
      self%time%it = self%time%it + 1

      if (self%io%save_memory_status) then
         call save_memory_status(file_name='memory_cpu-'//self%mpih%myrankstr//'.dat', tag=str(self%time%it,.true.))
      endif

      if (mod(self%time%it,self%amr%frequency)==0) then
         call self%mpih%barrier(tictoc=.true.)
         call self%amr_update
         call self%mpih%barrier(tictoc=.true.)
      endif

      call self%compute_dt
      if ((self%time%it_max <= 0).and.(self%time%time+self%time%dt > self%time%time_max)) &
         self%time%dt=self%time%time_max-self%time%time

      call self%integrate

      self%time%time = self%time%time + self%time%dt
      call self%time%print_progress(nodes_number=self%adam%tree%nodes_number)

      call self%save_simulation_data

      if (((self%time%it_max <= 0).and.(self%time%time >= self%time%time_max)).or.&
         ((self%time%it>=self%time%it_max).and.(self%time%it_max > 0))) exit integration

      call self%mpih%barrier(tictoc=.true., timing=timing_step(2), single=.true.)
   enddo integration
   call self%mpih%barrier(tictoc=.true., timing=timing(2), single=.true.)
   call self%save_simulation_data
   if (self%mpih%myrank==0) call self%io%close_file_residuals
   call self%mpih%finalize
   endsubroutine simulate

   ! non TBP
   subroutine assign_omp_R8P_5D(blocks_number, ngc, lhs, rhs)
   !< Assign array to scalar value with OpenMP threads (kind R8P, rank 5)
   integer(I4P), intent(in)    :: blocks_number                   !< Number of blocks.
   integer(I4P), intent(in)    :: ngc                             !< Ghost cells number.
   real(R8P),    intent(inout) :: lhs(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Lest hand side.
   real(R8P),    intent(in)    :: rhs                             !< Right hand side.
   integer(I4P)                :: ni, nj, nk, nv, b, i, j, k, v   !< Counter.

   nv = ubound(lhs,dim=1)
   ni = ubound(lhs,dim=2) - ngc
   nj = ubound(lhs,dim=3) - ngc
   nk = ubound(lhs,dim=4) - ngc
   !$omp parallel do collapse(5) default(firstprivate) shared(lhs)
   do b=1, blocks_number
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do v=1    , nv
      lhs(v,i,j,k,b) = rhs
   enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine assign_omp_R8P_5D

   subroutine compute_fluxes_convective(dir,blocks_number,ni,nj,nk,ngc,nv,g,weno_s,weno_zeps,weno_a,weno_p,weno_d,q_aux,fluxes)
   !< Compute convective fluxes along direction `dir`.
   integer(I4P), intent(in)    :: dir                                !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)    :: blocks_number                      !< Number of blocks.
   integer(I4P), intent(in)    :: ni                                 !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                 !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                 !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                                 !< Number of conservative varibales.
   real(R8P),    intent(in)    :: g                                  !< Specific heats ratio.
   integer(I4P), intent(in)    :: weno_s                             !< Weno stencils number/dimension.
   real(R8P),    intent(in)    :: weno_zeps                          !< Parameter for avoiding division by zero in computing IS.
   real(R8P),    intent(in)    :: weno_a(1:,0:,1:)                   !< Optimal weights.
   real(R8P),    intent(in)    :: weno_p(1:,0:,0:,1:)                !< Polinomials coefficients.
   real(R8P),    intent(in)    :: weno_d(0:,0:,0:,1:)                !< Smoothness indicators coefficients.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Auxiliary variables.
   real(R8P),    intent(inout) :: fluxes(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Fluxes.
   integer(I4P)                :: si(3), si_i, si_j, si_k            !< Directional (1=x,2=y,3=z) increment.
   real(R8P)                   :: sir(3)                             !< Directional (1=x,2=y,3=z) increment, real.
   integer(I4P)                :: b, i, j, k                         !< Counter.

   select case(dir)
   case(1)
      si = [1,0,0]
   case(2)
      si = [0,1,0]
   case(3)
      si = [0,0,1]
   endselect
   sir = real(si,R8P)
   si_i = 1-si(1)
   si_j = 1-si(2)
   si_k = 1-si(3)

   do b=1, blocks_number
   do k=si_k, nk
   do j=si_j, nj
   do i=si_i, ni
      call compute_fluxes_convective_ri(dir=dir,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv,                                       &
                                        weno_s=weno_s,weno_zeps=weno_zeps,weno_a=weno_a,weno_p=weno_p, weno_d=weno_d,&
                                        si=si,sir=sir,g=g,q_aux=q_aux,fluxes=fluxes)
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_fluxes_convective

   subroutine compute_fluxes_convective_ri(dir,b,i,j,k,ngc,nv,                   &
                                           weno_s,weno_zeps,weno_a,weno_p,weno_d,&
                                           si,sir,g,q_aux,fluxes)
   !< Compute convective fluxes at right interface of b,i,j,k.
   integer(I4P), intent(in)    :: dir                                 !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)    :: b, i, j, k                          !< Counter.
   integer(I4P), intent(in)    :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                                  !< Number of conservative varibales.
   integer(I4P), intent(in)    :: weno_s                              !< Weno stencils number/dimension.
   real(R8P),    intent(in)    :: weno_zeps                           !< Parameter to avoid division by zero.
   real(R8P),    intent(in)    :: weno_a(1:,0:,1:)                    !< Optimal weights.
   real(R8P),    intent(in)    :: weno_p(1:,0:,0:,1:)                 !< Polinomials coefficients.
   real(R8P),    intent(in)    :: weno_d(0:,0:,0:,1:)                 !< Smoothness indicators coefficients.
   integer(I4P), intent(in)    :: si(3)                               !< Stencil increment.
   real(R8P),    intent(in)    :: sir(3)                              !< Stencil increment, real cast.
   real(R8P),    intent(in)    :: g                                   !< Specific heats ratio.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Auxiliary variables.
   real(R8P),    intent(inout) :: fluxes(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Fluxes.
   real(R8P)                   :: fmpc(1:2,1-weno_s:-1+weno_s,1:nv)   !< Fluxes -+ decomposition in c. space.
   real(R8P)                   :: fpmr(1:2,1:nv)                      !< Fluxes +- reconstructed.
   real(R8P)                   :: el(nv,nv), er(nv,nv)                !< Left and right eigenvalues.
   integer(I4P)                :: v, vv                               !< Counter.

   call compute_eigenvectors(si=si,sir=sir,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv,g=g,q_aux=q_aux,el=el,er=er)
   call decompose_fluxes_convective(weno_s=weno_s,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv,si=si,sir=sir,el=el,q_aux=q_aux,fmpc=fmpc)
   do v=1, nv
      call weno_reconstruct_upwind(S=weno_s, weno_a=weno_a, weno_p=weno_p, weno_d=weno_d,&
                                   weno_zeps=weno_zeps, V=fmpc(:,:,v), VR=fpmr(:,v))
   enddo
   ! back projection in conservative variables space
   do v=1, nv
      fluxes(v,i,j,k,b) = 0._R8P
      do vv=1,nv
         fluxes(v,i,j,k,b) = fluxes(v,i,j,k,b) + er(vv,v) * (fpmr(1,vv) + fpmr(2,vv))
      enddo
   enddo
   endsubroutine compute_fluxes_convective_ri

   subroutine compute_fluxes_difference(null_xyz, blocks_number, ni, nj, nk, ngc, nv, ib_eps, dx, dy, dz, flx, fly, flz, phi, dq)
   !< Compute fluxes difference.
   logical,      intent(in)           :: null_xyz(3)                     !< Nullified directions tags.
   integer(I4P), intent(in)           :: blocks_number                   !< Number of blocks.
   integer(I4P), intent(in)           :: ni                              !< Grid cells number in I direction.
   integer(I4P), intent(in)           :: nj                              !< Grid cells number in J direction.
   integer(I4P), intent(in)           :: nk                              !< Grid cells number in K direction.
   integer(I4P), intent(in)           :: ngc                             !< Ghost cells number.
   integer(I4P), intent(in)           :: nv                              !< Number of conservative varibales.
   real(R8P),    intent(in)           :: ib_eps                          !< Tolerance IB delta ratio.
   real(R8P),    intent(in)           :: dx(1:), dy(1:), dz(1:)          !< Space steps.
   real(R8P),    intent(in)           :: flx(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< X direction fluxes.
   real(R8P),    intent(in)           :: fly(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Y direction fluxes.
   real(R8P),    intent(in)           :: flz(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Z direction fluxes.
   real(R8P),    intent(in), optional :: phi(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< IB distance function.
   real(R8P),    intent(inout)        ::  dq(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Fluxes differences.
   real(R8P)                          :: delta_x, delta_y, delta_z       !< Space steps.
   real(R8P)                          :: dx_locale, dy_locale, dz_locale !< Local space steps.
   integer(I4P)                       :: b, i, j, k, v                   !< Counter.
   integer(I4P)                       :: all_solids                      !< Last phi index, all solids summary.
   real(R8P)                          :: qmx, qmy, qmz                   !< Momentum nullification scalar.

   qmx = 1._R8P ; if (null_xyz(1)) qmx = 0._R8P
   qmy = 1._R8P ; if (null_xyz(2)) qmy = 0._R8P
   qmz = 1._R8P ; if (null_xyz(3)) qmz = 0._R8P
   if (present(phi)) then
      all_solids = ubound(phi, dim=1)
      !$omp parallel do collapse(4) default(firstprivate) shared(dx,dy,dz,flx,fly,flz,phi,dq)
      do b=1,blocks_number
      do k=1,nk
      do j=1,nj
      do i=1,ni
         dx_locale = dx(b)
         if (phi(all_solids,i,j,k,b)<0.) then
            if (phi(all_solids,i+1,j,k,b)*phi(all_solids,i-1,j,k,b)<0) then
               if (phi(all_solids,i+1,j,k,b)>0.) then
                  delta_x = -phi(all_solids,i,j,k,b)/(phi(all_solids,i+1,j,k,b)-phi(all_solids,i,j,k,b)+ib_eps)*dx(b)
                  dx_locale = dx(b)/2 + delta_x
               else
                  delta_x = -phi(all_solids,i,j,k,b)/(phi(all_solids,i-1,j,k,b)-phi(all_solids,i,j,k,b)+ib_eps)*dx(b)
                  dx_locale = dx(b)/2 + delta_x
               endif
            endif
         endif
         dy_locale = dy(b)
         if (phi(all_solids,i,j,k,b)<0.) then
            if (phi(all_solids,i,j+1,k,b)*phi(all_solids,i,j-1,k,b)<0) then
               if (phi(all_solids,i,j+1,k,b)>0.) then
                  delta_y = -phi(all_solids,i,j,k,b)/(phi(all_solids,i,j+1,k,b)-phi(all_solids,i,j,k,b)+ib_eps)*dy(b)
                  dy_locale = dy(b)/2 + delta_y
               else
                  delta_y = -phi(all_solids,i,j,k,b)/(phi(all_solids,i,j-1,k,b)-phi(all_solids,i,j,k,b)+ib_eps)*dy(b)
                  dy_locale = dy(b)/2 + delta_y
               endif
            endif
         endif
         dz_locale = dz(b)
         if (phi(all_solids,i,j,k,b)<0.) then
            if (phi(all_solids,i,j,k+1,b)*phi(all_solids,i,j,k-1,b)<0) then
               if (phi(all_solids,i,j,k+1,b)>0.) then
                  delta_z = -phi(all_solids,i,j,k,b)/(phi(all_solids,i,j,k+1,b)-phi(all_solids,i,j,k,b)+ib_eps)*dz(b)
                  dz_locale = dz(b)/2 + delta_z
               else
                  delta_z = -phi(all_solids,i,j,k,b)/(phi(all_solids,i,j,k-1,b)-phi(all_solids,i,j,k,b)+ib_eps)*dz(b)
                  dz_locale = dz(b)/2 + delta_z
               endif
            endif
         endif
         do v=1, nv
            dq(v,i,j,k,b) = - (flx(v,i,j,k,b)-flx(v,i-1,j,k,b))/dx_locale &
                            - (fly(v,i,j,k,b)-fly(v,i,j-1,k,b))/dy_locale &
                            - (flz(v,i,j,k,b)-flz(v,i,j,k-1,b))/dz_locale
         enddo
         dq(2,i,j,k,b) = dq(2,i,j,k,b) * qmx
         dq(3,i,j,k,b) = dq(3,i,j,k,b) * qmy
         dq(4,i,j,k,b) = dq(4,i,j,k,b) * qmz
      enddo
      enddo
      enddo
      enddo
      !$omp end parallel do
   else
      !$omp parallel do collapse(4) default(firstprivate) shared(dx,dy,dz,flx,fly,flz,phi,dq)
      do b=1,blocks_number
      do k=1,nk
      do j=1,nj
      do i=1,ni
         do v=1, nv
            dq(v,i,j,k,b) = - (flx(v,i,j,k,b)-flx(v,i-1,j,k,b))/dx(b) &
                            - (fly(v,i,j,k,b)-fly(v,i,j-1,k,b))/dy(b) &
                            - (flz(v,i,j,k,b)-flz(v,i,j,k-1,b))/dz(b)
         enddo
         dq(2,i,j,k,b) = dq(2,i,j,k,b) * qmx
         dq(3,i,j,k,b) = dq(3,i,j,k,b) * qmy
         dq(4,i,j,k,b) = dq(4,i,j,k,b) * qmz
      enddo
      enddo
      enddo
      enddo
      !$omp end parallel do
   endif
   endsubroutine compute_fluxes_difference

   pure subroutine compute_eigenvectors(si,sir,b,i,j,k,ngc,nv,g,q_aux,el,er)
   ! Compute eigenvectors centered in inteface i,j,k/ip,jp,kp.
   integer(I4P), intent(in)    :: si(3)                             !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)    :: sir(3)                            !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)    :: b, i, j, k                        !< Cell indexes.
   integer(I4P), intent(in)    :: ngc                               !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                                !< Number of conservative varibales.
   real(R8P),    intent(in)    :: g                                 !< Specific heats ratio.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: el(1:nv,1:nv), er(1:nv,1:nv)      !< Left and right eigenvectors.
   real(R8P)                   :: uu, vv, ww, h, qq, c, ci, b1, b2  !< Roe average states.
   real(R8P)                   :: uvw, uvw_r1, uvw_r2               !< Velocity rotation accordingly dir.
   real(R8P), parameter        :: o_2=0.5_R8P                       !< Alias 1/2.
   real(R8P)                   :: ci1,ci2,ci3                       !< Alias sir(1)/c, sir(2)/c, sir(3)/c.

   call compute_roe_average(q_aux=q_aux, g=g, ngc=ngc, b=b, i=i, j=j, k=k, ip=i+si(1), jp=j+si(2), kp=k+si(3), &
                            uu=uu, vv=vv, ww=ww, h=h, qq=qq, c=c, ci=ci, b1=b1, b2=b2)

   uvw    = uu*sir(1)+vv*sir(2)+ww*sir(3)
   uvw_r1 = (vv-uvw*sir(2))*sir(1)+(uvw*sir(1)-uu)*sir(2)+(uu-uvw*sir(1))*sir(3)
   uvw_r2 = (uvw*sir(3)-ww)*sir(1)+(ww-uvw*sir(3))*sir(2)+(uvw*sir(2)-vv)*sir(3)
   ci1    = ci*sir(1)
   ci2    = ci*sir(2)
   ci3    = ci*sir(3)

   er(1,1)= 1._R8P      ; er(1,2)= 1._R8P ; er(1,3)= 1._R8P      ; er(1,4)= 0._R8P              ; er(1,5)= 0._R8P
   er(2,1)= uu-c*sir(1) ; er(2,2)= uu     ; er(2,3)= uu+c*sir(1) ; er(2,4)= sir(2)              ; er(2,5)=-sir(3)
   er(3,1)= vv-c*sir(2) ; er(3,2)= vv     ; er(3,3)= vv+c*sir(2) ; er(3,4)=-sir(1)              ; er(3,5)= 0._R8P
   er(4,1)= ww-c*sir(3) ; er(4,2)= ww     ; er(4,3)= ww+c*sir(3) ; er(4,4)= 0._R8P              ; er(4,5)= sir(1)
   er(5,1)= h-uvw*c     ; er(5,2)= qq     ; er(5,3)= h+uvw*c     ; er(5,4)= uu*sir(2)-vv*sir(1) ; er(5,5)= ww*sir(1)-uu*sir(3)

   el(1,1)= o_2*(b1+uvw*ci);el(1,2)=-o_2*(b2*uu+ci1);el(1,3)=-o_2*(b2*vv+ci2);el(1,4)=-o_2*(b2*ww+ci3);el(1,5)= o_2*b2
   el(2,1)= 1._R8P-b1      ;el(2,2)= b2*uu          ;el(2,3)= b2*vv          ;el(2,4)= b2*ww          ;el(2,5)=-b2
   el(3,1)= o_2*(b1-uvw*ci);el(3,2)=-o_2*(b2*uu-ci1);el(3,3)=-o_2*(b2*vv-ci2);el(3,4)=-o_2*(b2*ww-ci3);el(3,5)= o_2*b2
   el(4,1)= uvw_r1         ;el(4,2)= sir(2)-sir(3)  ;el(4,3)=-sir(1)         ;el(4,4)= 0._R8P         ;el(4,5)= 0._R8P
   el(5,1)= uvw_r2         ;el(5,2)= 0._R8P         ;el(5,3)= sir(3)         ;el(5,4)= sir(1)-sir(2)  ;el(5,5)= 0._R8P
   endsubroutine compute_eigenvectors

   pure subroutine compute_max_eigenvalues(si,weno_s,b,i,j,k,ngc,nv,q_aux,evmax)
   ! Compute maximum eigenvalues in the big stencil.
   integer(I4P), intent(in)  :: si(3)                             !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)  :: weno_s                            !< Weno stencils number/dimension.
   integer(I4P), intent(in)  :: b, i, j, k                        !< Cell indexes.
   integer(I4P), intent(in)  :: ngc                               !< Ghost cells number.
   integer(I4P), intent(in)  :: nv                                !< Number of conservative varibales.
   real(R8P),    intent(in)  :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(out) :: evmax                             !< Maximum eigenvalue in the big stencil.
   integer(I4P)              :: s, is, js, ks                     !< Counter.

   evmax = -1._R8P
   do s=1, 2*weno_s
      is = i + (s-weno_s) * si(1) ; js = j + (s-weno_s) * si(2) ; ks = k + (s-weno_s) * si(3)
      evmax = max(abs(q_aux(VAR_U*si(1)+VAR_V*si(2)+VAR_W*si(3),is,js,ks,b)) + &
                      q_aux(VAR_A                              ,is,js,ks,b), evmax)
   enddo
   endsubroutine compute_max_eigenvalues

   subroutine compute_q_gradient(b, ni, nj, nk, ngc, dx, dy, dz, q, ivar, gradient)
   !< Compute gradient of q(ivar).
   integer(I4P), intent(in)  :: b                             !< Block index.
   integer(I4P), intent(in)  :: ni                            !< Grid cells number in I direction.
   integer(I4P), intent(in)  :: nj                            !< Grid cells number in J direction.
   integer(I4P), intent(in)  :: nk                            !< Grid cells number in K direction.
   integer(I4P), intent(in)  :: ngc                           !< Ghost cells number.
   real(R8P),    intent(in)  :: dx                            !< X space step.
   real(R8P),    intent(in)  :: dy                            !< Y space step.
   real(R8P),    intent(in)  :: dz                            !< Z space step.
   real(R8P),    intent(in)  :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field component to which apply gradient.
   integer(I4P), intent(in)  :: ivar                          !< Index of variable for computing the gradient.
   real(R8P),    intent(out) :: gradient                      !< Maximum gradient of q.
   real(R8P)                 :: grad                          !< Current gradient of q.
   integer(I4P)              :: i, j, k                       !< Counter.
   real(R8P), parameter      :: tol=1.e-12                    !< Gradient denominator tolerance.

   gradient = 0._R8P
   !$omp parallel do collapse(3) default(firstprivate) shared(q) reduction(max:gradient)
   do k=1, nk
      do j=1, nj
         do i=1, ni
            grad = sqrt(((q(ivar,i+1,j,k,b) - q(ivar,i-1,j,k,b))/(2*dx))**2 + &
                        ((q(ivar,i,j+1,k,b) - q(ivar,i,j-1,k,b))/(2*dy))**2 + &
                        ((q(ivar,i,j,k+1,b) - q(ivar,i,j,k-1,b))/(2*dz))**2)
            grad = grad/(abs(q(ivar,i,j,k,b))+tol)
            gradient = max(gradient, grad)
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine compute_q_gradient

   pure subroutine compute_roe_average(ngc, b, i, j, k, ip, jp, kp, g, q_aux, &
                                       uu, vv, ww, h, qq, c, ci, b1, b2)
   !< Compute Roe averaged quantities.
   integer(I4P), intent(in)  :: ngc                               !< Number of ghost cells.
   integer(I4P), intent(in)  :: b, i, j, k, ip, jp, kp            !< Left/right cells indexes.
   real(R8P),    intent(in)  :: g                                 !< Specific heats ratio.
   real(R8P),    intent(in)  :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(out) :: uu, vv, ww, h, qq, c, ci, b1, b2  !< Roe state average variables.
   real(R8P)                 :: ri, up, vp, wp, hp, r, rp1, cc    !< Local varbiables.

   ! left state (node i)
   ri = 1._R8P/q_aux(VAR_R,i,j,k,b)
   uu = q_aux(VAR_U,i,j,k,b)
   vv = q_aux(VAR_V,i,j,k,b)
   ww = q_aux(VAR_W,i,j,k,b)
   h  = q_aux(VAR_H,i,j,k,b)
   ! right state (node i+1)
   up = q_aux(VAR_U,ip,jp,kp,b)
   vp = q_aux(VAR_V,ip,jp,kp,b)
   wp = q_aux(VAR_W,ip,jp,kp,b)
   hp = q_aux(VAR_H,ip,jp,kp,b)
   ! Roe average state
   r   = sqrt(q_aux(VAR_R,ip,jp,kp,b)*ri)
   rp1 = 1._R8P/(r+1._R8P)
   uu  = (r*up+uu)*rp1
   vv  = (r*vp+vv)*rp1
   ww  = (r*wp+ww)*rp1
   h   = (r*hp+h)*rp1
   qq  = 0.5_R8P * (uu*uu+vv*vv+ww*ww)
   cc  = (g-1._R8P) * (h - qq)
   c   = sqrt(cc)
   ci  = 1._R8P/c
   b2  = (g-1._R8P)/cc  ! alias 1/(cp*theta)
   b1  = b2 * qq        ! alias q/(cp*theta)
   endsubroutine compute_roe_average

   pure subroutine decompose_fluxes_convective(weno_s,b,i,j,k,ngc,nv,si,sir,el,q_aux,fmpc)
   !< Decompose convective fluxes.
   !< Flux vector splitting by local-Lax-Friedrics (Rusanov) with projection in pseudo-characteristics space.
   integer(I4P), intent(in)    :: weno_s                            !< Weno stencils number/dimension.
   integer(I4P), intent(in)    :: b, i, j, k                        !< Counter.
   integer(I4P), intent(in)    :: ngc                               !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                                !< Number of conservative varibales.
   integer(I4P), intent(in)    :: si(3)                             !< Stencil increment.
   real(R8P),    intent(in)    :: sir(3)                            !< Stencil increment, real cast.
   real(R8P),    intent(in)    :: el(1:,1:)                         !< Left eigeinvectors.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: fmpc(1:,1-weno_s:,1:)             !< Fluxes -+ decomposition in characteristics space.
   real(R8P)                   :: fmp(2)                            !< Fluxes -+ decomposition in each cell stencils.
   real(R8P)                   :: evmax                             !< Maximum eigenvalue.
   real(R8P)                   :: q(nv), f(nv)                      !< Conservative variables and fluxes.
   real(R8P)                   :: gc, wc                            !< Increments for fluxes decomposition.
   integer(I4P)                :: v, vv, s, is, js, ks              !< Counter.

   call compute_max_eigenvalues(si=si,weno_s=weno_s,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv,q_aux=q_aux,evmax=evmax)
   do s=1-weno_s, weno_s
      is = i + (s) * si(1) ; js = j + (s) * si(2) ; ks = k + (s) * si(3)
      call compute_conservative(auxiliary=q_aux(:,is,js,ks,b),conservative=q)
      call compute_fluxes_conservative(sir=sir,auxiliary=q_aux(:,is,js,ks,b),fluxes=f)
      do v=1, nv
         wc = 0._R8P
         gc = 0._R8P
         do vv=1, nv
            wc = wc + el(vv,v) * q(vv)
            gc = gc + el(vv,v) * f(vv)
         enddo
         fmp(2) = 0.5_R8P * (gc + evmax * wc)
         fmp(1) = gc - fmp(2)
         if (s<weno_s)   fmpc(2,s  ,v) = fmp(2)
         if (s>1-weno_s) fmpc(1,s-1,v) = fmp(1)
      enddo
   enddo
   endsubroutine decompose_fluxes_convective
endmodule adam_chase_cpu_object
