!< ADAM, Navier-Stokes equations system class definition, CPU backend.
module adam_nasto_cpu_object
!< ADAM, Navier-Stokes equations system class definition, CPU backend.

use adam_common_library
use adam_nasto_common_library
use penf
use mpi

implicit none
private
public :: nasto_cpu_object

type, extends(nasto_common_object) :: nasto_cpu_object
   !< Navier-Stokes equations system class definition, CPU backend.
   real(R8P), allocatable :: fl(:,:,:,:,:)       !< Residuals.
   real(R8P), allocatable :: flx(:,:,:,:,:)      !< Fluxes along x.
   real(R8P), allocatable :: fly(:,:,:,:,:)      !< Fluxes along y.
   real(R8P), allocatable :: flz(:,:,:,:,:)      !< Fluxes along z.
   real(R8P), allocatable :: gplus_x(:,:,:,:,:)  !< Positive fluxes for weno-x.
   real(R8P), allocatable :: gminus_x(:,:,:,:,:) !< Negative fluxes for weno-x.
   real(R8P), allocatable :: gplus_y(:,:,:,:,:)  !< Positive fluxes for weno-y.
   real(R8P), allocatable :: gminus_y(:,:,:,:,:) !< Negative fluxes for weno-y.
   real(R8P), allocatable :: gplus_z(:,:,:,:,:)  !< Positive fluxes for weno-z.
   real(R8P), allocatable :: gminus_z(:,:,:,:,:) !< Negative fluxes for weno-z.
   real(R8P), allocatable :: q_old(:,:,:,:,:)    !< Field cell centered variables (old iteration).
   contains
      ! auxiliary methods
      procedure, pass(self) :: allocate_cpu !< Allocate CPU data.
      procedure, pass(self) :: destroy      !< Destroy the equation.
      procedure, pass(self) :: initialize   !< Initialize the equation.
      ! AMR methods
      procedure, pass(self) :: amr_update       !< Do AMR update.
      procedure, pass(self) :: compute_phi      !< Compute phi, distance from IB solid.
      procedure, pass(self) :: mark_by_geo      !< Mark blocks to be refined/derefined by a geometric constrain.
      procedure, pass(self) :: mark_by_grad_var !< Mark blocks to be refined/derefined by a `grad(var)` value.
      procedure, pass(self) :: move_phi         !< Move phi.
      procedure, pass(self) :: refine_uniform   !< Refine all blocks uniformly.
      ! IB methods
      procedure, pass(self) :: integrate_eikonal_q !< Integrate eikonal equation over q.
      procedure, pass(self) :: invert_eikonal_q    !< Invert momentum eikonal equation over q.
      ! IO methods
      procedure, pass(self) :: load_restart_files   !< Load restart files.
      procedure, pass(self) :: save_hdf5            !< Save simulation data in HDF5 format.
      procedure, pass(self) :: save_residuals       !< Save residuals history.
      procedure, pass(self) :: save_restart_files   !< Save restart files.
      procedure, pass(self) :: save_simulation_data !< Save all simulation data.
      ! IC/BC
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions of equation.
      procedure, pass(self) :: set_initial_conditions  !< Set initial conditions of equation.
      procedure, pass(self) :: update_ghost            !< Update ghost cells and set boundary conditions.
      ! numerical methods
      procedure, pass(self) :: compute_dt        !< Compute time step.
      procedure, pass(self) :: compute_q_aux     !< Compute auxiliary variables.
      procedure, pass(self) :: compute_residuals !< Compute residuals.
      procedure, pass(self) :: compute_rk_q      !< Compute RK approximation over q.
      procedure, pass(self) :: integrate         !< Perform one step integration.
      procedure, pass(self) :: simulate          !< Perform the simulation.
endtype nasto_cpu_object
contains
   ! auxiliary methods
   subroutine allocate_cpu(self)
   !< Allocate CPU data.
   class(nasto_cpu_object), intent(inout) :: self !< The equation.
   character(:), allocatable              :: msg_ !< Allocating message base.
   character(:), allocatable              :: msg  !< Allocating message.
   integer(I4P)                           :: sn1  !< Solids number + 1.

   call self%mpih%print_message('nasto_cpu_object%allocate_cpu start')
   msg_ = self%mpih%myrankstr//'nasto_cpu_object%allocate_cpu '
   associate(nv=>self%nv, ns=>self%ns, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             nb=>self%nb, nv_aux=>self%nv_aux, iweno=>self%schemes%iweno, solids_number=>self%ib%solids_number)
   msg = msg_//' fl '
   call alloc_var_cpu(var=self%fl,       ulb=reshape([1,nv,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,5]),msg=msg)
   self%fl = 0._R8P
   msg = msg_//' flx '
   call alloc_var_cpu(var=self%flx,      ulb=reshape([1,nv,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,5]),msg=msg)
   self%flx = 0._R8P
   msg = msg_//' fly '
   call alloc_var_cpu(var=self%fly,      ulb=reshape([1,nv,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,5]),msg=msg)
   self%fly = 0._R8P
   msg = msg_//' flz '
   call alloc_var_cpu(var=self%flz,      ulb=reshape([1,nv,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,5]),msg=msg)
   self%flz = 0._R8P
   msg = msg_//' gplus_x '
   call alloc_var_cpu(var=self%gplus_x , ulb=reshape([1,nv,1,2*iweno,   1,nj,        1,nk,        1,nb],[2,5]),msg=msg)
   self%gplus_x = 0._R8P
   msg = msg_//' gminus_x '
   call alloc_var_cpu(var=self%gminus_x, ulb=reshape([1,nv,1,2*iweno,   1,nj,        1,nk,        1,nb],[2,5]),msg=msg)
   self%gminus_x = 0._R8P
   msg = msg_//' gplus_y '
   call alloc_var_cpu(var=self%gplus_y , ulb=reshape([1,nv,1,2*iweno,   1,ni,        1,nk,        1,nb],[2,5]),msg=msg)
   self%gplus_y = 0._R8P
   msg = msg_//' gminus_y '
   call alloc_var_cpu(var=self%gminus_y, ulb=reshape([1,nv,1,2*iweno,   1,ni,        1,nk,        1,nb],[2,5]),msg=msg)
   self%gminus_y = 0._R8P
   msg = msg_//' gplus_z '
   call alloc_var_cpu(var=self%gplus_z , ulb=reshape([1,nv,1,2*iweno,   1,ni,        1,nj,        1,nb],[2,5]),msg=msg)
   self%gplus_z = 0._R8P
   msg = msg_//' gminus_z '
   call alloc_var_cpu(var=self%gminus_z, ulb=reshape([1,nv,1,2*iweno,   1,ni,        1,nj,        1,nb],[2,5]),msg=msg)
   self%gminus_z = 0._R8P
   msg = msg_//' q_old '
   call alloc_var_cpu(var=self%q_old,    ulb=reshape([1,nv,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,5]),msg=msg)
   self%q_old = 0._R8P
   endassociate
   call self%mpih%print_message('nasto_cpu_object%allocate_cpu finish')
   endsubroutine allocate_cpu

   subroutine destroy(self)
   !< Destroy the equation.
   class(nasto_cpu_object), intent(inout) :: self  !< The equation.
   type(nasto_cpu_object)                 :: fresh !< Fresh equation.

   ! TODO to be implemented
   endsubroutine destroy

   subroutine initialize(self, filename)
   !< Initialize the equation.
   class(nasto_cpu_object), intent(inout) :: self     !< The equation.
   character(*),            intent(in)    :: filename !< Input file name.

   call self%mpih%initialize(do_mpi_init=.true.)
   call self%mpih%print_message('nasto_cpu_object%initialize start')
   call self%initialize_common(filename=filename, memory_avail=self%mpih%memory_avail)
   call self%allocate_cpu
   call self%mpih%print_message('nasto_cpu_object%initialize finish')
   endsubroutine initialize

   ! AMR methods
   subroutine amr_update(self)
   !< Do AMR update.
   class(nasto_cpu_object), intent(inout) :: self                 !< The equation.
   integer(I4P)                           :: iterations_          !< Number of AMR iterations, local var.
   logical                                :: is_grid_changed      !< Flag to check grid changes for each marker.
   logical                                :: is_grid_changed_all  !< Flag to check grid changes for each iter.
   integer(I4P)                           :: b, i, j, k, i_marker !< Counter.
   type(amr_marker_object)                :: amr_marker           !< Current amr marker.

   amr: do i=1, self%amr%iters
      is_grid_changed_all = .false.
      do i_marker=1, self%amr%markers_number
         amr_marker = self%amr%markers(i_marker)
         call self%update_ghost(q=self%field%q)
         select case(amr_marker%mode)
         case(AMR_GEO)
            call self%mark_by_geo(delta_fine=amr_marker%delta_fine, delta_coarse=amr_marker%delta_coarse)
         case(AMR_GRAD)
            select case(amr_marker%field)
            case(1)
               call self%mark_by_grad_var(grad_tol=amr_marker%tol, delta_fine=amr_marker%delta_fine, &
                                          delta_coarse=amr_marker%delta_coarse, ivar=amr_marker%ivar)
            case(2)
               call self%mark_by_grad_var(grad_tol=amr_marker%tol, delta_fine=amr_marker%delta_fine, &
                                          delta_coarse=amr_marker%delta_coarse, ivar=amr_marker%ivar)
            endselect
         endselect
         call self%adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false., is_grid_changed=is_grid_changed)
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
   class(nasto_cpu_object), intent(inout) :: self !< The equation.

   if (self%ib%solids_number>0) then
      call self%mpih%print_message('compute IB distance start')
      call self%ib%compute_phi
      call self%mpih%print_message('compute IB distance finish')
   endif
   endsubroutine compute_phi

   subroutine mark_by_geo(self, delta_fine, delta_coarse, threshold, do_init)
   !< Mark blocks to be refined/derefined by a geometric constrain.
   class(nasto_cpu_object), intent(inout)        :: self           !< The equation.
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
             blocks_number=>self%blocks_number, ns=>self%ns, dxyz=>self%field%dxyz, phi=>self%ib%phi)
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

   subroutine mark_by_grad_var(self, grad_tol, delta_fine, delta_coarse, ivar, threshold, do_init)
   !< Mark blocks to be refined/derefined by a `grad(var)` value.
   class(nasto_cpu_object), intent(inout)        :: self           !< The equation.
   real(R8P),               intent(in)           :: grad_tol       !< Gradiend tolerance value.
   real(R8P),               intent(in)           :: delta_fine     !< Maximum cell delta in fine grids.
   real(R8P),               intent(in)           :: delta_coarse   !< Minimum cell delta in coarse grids.
   integer(I4P),            intent(in), optional :: ivar           !< Variable for marking.
   real(R8P),               intent(in), optional :: threshold      !< Threshold for sphere proximity.
   logical,                 intent(in), optional :: do_init        !< Re-initialize refinements queries.
   integer(I4P)                                  :: ivar_          !< Variable for marking (local var).
   logical                                       :: do_init_       !< Re-initialize refinements queries, local var.
   real(R8P)                                     :: threshold_     !< Threshold for sphere proximity, local var.
   real(R8P)                                     :: max_cell_delta !< Maximum cell delta.
   real(R8P)                                     :: grad_var       !< Value (max) of gradient of var.
   integer(I4P)                                  :: b              !< Counter.

   ivar_     = 1_R4P    ; if (present(ivar)) ivar_ = ivar
   do_init_ = .true.    ; if (present(do_init)) do_init_ = do_init
   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   if (do_init_) self%field%refinements_needed = [(TO_BE_DEREFINED,b=1,self%blocks_number)]
   associate (ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
              blocks_number=>self%blocks_number, ns=>self%ns, dxyz=>self%field%dxyz)
      call self%update_ghost(q=self%field%q)
      call self%compute_q_aux(q=self%field%q, q_aux=self%q_aux)
      do b=1, blocks_number
         call compute_q_gradient(b=b, ni=ni, nj=nj, nk=nk, ngc=ngc, &
                                 dx=dxyz(1,b), dy=dxyz(2,b), dz=dxyz(3,b), q=self%q_aux, ivar=ivar_, gradient=grad_var)
         max_cell_delta = max_cell_delta_grad(grad=grad_var)
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
   !< Move phi and the actual ptree representation.
   class(nasto_cpu_object), intent(inout) :: self        !< The equation.
   real(R8P),               intent(in)    :: velocity(3) !< Velocity of the movement.
   integer(I4P),            intent(in)    :: s           !< Solid index.

   if (self%ib%solids_number>0) then
      call self%mpih%print_message('compute IB distance start')
      call self%ib%move_phi(velocity=velocity, s=s)
      call self%mpih%print_message('compute IB distance finish')
   endif
   endsubroutine move_phi

   subroutine refine_uniform(self, refinement_levels)
   !< Refine all blocks uniformly.
   class(nasto_cpu_object), intent(inout) :: self              !< The equation.
   integer(I4P),            intent(in)    :: refinement_levels !< Number of refinement to be performed.
   integer(I4P)                           :: l                 !< Counter.

   do l=1, refinement_levels
      call self%adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
      call self%adam%amr_update(do_blocks_reorder=.false., do_mpi_redistribute=.true.)
   enddo
   endsubroutine

   ! IB methods
   subroutine integrate_eikonal_q(self)
   !< Integrate eikonal equation over q.
   class(nasto_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: ib   !< Counter.

   associate(blocks_number=>self%blocks_number, solids_number=>self%ib%solids_number, q=>self%field%q)
   if (blocks_number > 0) then
      if (solids_number > 0) call self%ib%evolve_eikonal_q(q=q)
   endif
   endassociate
   endsubroutine integrate_eikonal_q

   subroutine invert_eikonal_q(self)
   !< Invert momentum eikonal equation over q.
   class(nasto_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: ib   !< Counter.

   associate(blocks_number=>self%blocks_number, solids_number=>self%ib%solids_number, q=>self%field%q)
   if (blocks_number > 0) then
      if (solids_number > 0) call self%ib%invert_eikonal_q(q=q)
   endif
   endassociate
   endsubroutine invert_eikonal_q

   ! IO methods
   subroutine load_restart_files(self, t, time)
   !< Save restart files.
   class(nasto_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P),            intent(out)   :: t    !< Time iteration.
   real(R8P),               intent(out)   :: time !< Time.

   call self%adam%load_restart_files(basename=self%io%restart_basename, t=t, time=time)
   call self%adam%make_comm_local_maps_ghost_bc
   endsubroutine load_restart_files

   subroutine save_hdf5(self, output_basename)
   !< Save simulation data in HDF5 format.
   class(nasto_cpu_object), intent(inout)        :: self             !< The equation.
   character(*),            intent(in), optional :: output_basename  !< Output basename.
   character(:), allocatable                     :: output_basename_ !< Output basename, local var.

   call self%mpih%barrier(tictoc=.true.)
   call self%mpih%print_message('save HDF5 files t: '//trim(str(self%time%it,.true.))//', time: '//&
                                trim(str(self%time%time,.true.)))
   output_basename_ = trim(self%io%output_basename)//'-'//trim(strz(self%time%it,9))
   if (present(output_basename)) output_basename_ = trim(output_basename)
   if (self%ib%solids_number>0) then
      call self%adam%save_hdf5(basename=trim(output_basename_),                                                       &
                               q=self%field%q,                                                                        &
                               q_aux=self%q_aux,                                                                      &
                               q_name=['rho','rhu','rhv','rhw','rhe'],                                                &
                               q_aux_name=['rhob ','u    ','v    ','w    ','ya   ','tem  ','pres ','ental','csp  '],  &
                               with_cell_morton=.true., phi=self%ib%phi)
   else
      call self%adam%save_hdf5(basename=trim(output_basename_),                                                       &
                               q=self%field%q,                                                                        &
                               q_aux=self%q_aux,                                                                      &
                               q_name=['rho','rhu','rhv','rhw','rhe'],                                                &
                               q_aux_name=['rhob ','u    ','v    ','w    ','ya   ','tem  ','pres ','ental','csp  '],  &
                               with_cell_morton=.true.)
   endif
   call self%mpih%barrier(tictoc=.true.)
   endsubroutine save_hdf5

   subroutine save_residuals(self)
   !< Save residuals history.
   class(nasto_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: v    !< Counter.

   if (self%time%is_to_save(it_save=self%io%residuals_save)) then
      call self%field%compute_normL2_residuals(dq=self%fl, norm=self%field%residuals)
      do v=1, self%nv
         call MPI_ALLREDUCE(MPI_IN_PLACE, self%field%residuals(v), 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, self%mpih%error)
         self%field%residuals(v) = sqrt(self%field%residuals(v))
      enddo
      if (self%mpih%myrank==0) call self%io%save_residuals(it=self%time%it, time=self%time%time, &
                                                           blocks_number=self%blocks_number, residuals=self%field%residuals)
   endif
   endsubroutine save_residuals

   subroutine save_restart_files(self)
   !< Save restart files.
   class(nasto_cpu_object), intent(inout) :: self !< The equation.

   call self%mpih%barrier(tictoc=.true.)
   call self%mpih%print_message('save restart files t: '//trim(str(self%time%it,.true.))//', time: '//&
                                trim(str(self%time%time,.true.)))
   call self%adam%save_restart_files(basename=self%io%restart_basename, t=self%time%it, time=self%time%time)
   call self%save_hdf5(output_basename=self%io%restart_basename)
   call self%mpih%barrier(tictoc=.true.)
   endsubroutine save_restart_files

   subroutine save_simulation_data(self)
   !< Save all simulation data.
   class(nasto_cpu_object), intent(inout) :: self !< The equation.

   if ((self%time%is_to_save(it_save=self%io%it_save)).or.      &
       (self%time%is_to_save(it_save=self%io%restart_save)).or. &
       (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max))) then
      call self%update_ghost(q=self%field%q)

      if (self%time%is_to_save(it_save=self%io%it_save)) call self%save_hdf5
      if (mod(self%time%it,self%io%restart_save)==0) call self%save_restart_files
      if (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max))&
         call self%slices%save_mat(basename=self%io%output_basename, &
                                   it=self%time%it,                  &
                                   it_max=self%time%it_max,          &
                                   time=self%time%time,              &
                                   time_max=self%time%time_max,      &
                                   adam=self%adam,                   &
                                   q=self%field%q,                   &
                                   q_name=['rho','rhu','rhv','rhw','rhe'])
   endif
   endsubroutine save_simulation_data

   ! IC/BC
   subroutine set_boundary_conditions(self, q)
   !< Set boundary conditions of equation.
   class(nasto_cpu_object), intent(in)    :: self              !< The equation.
   real(R8P),               intent(inout) :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,1:) !< Conservative variables.
   integer(I4P)                        :: b, c, i, j, k, v     !< Counter.
   integer(I4P)                        :: idelta,jdelta,kdelta !< IJK delta step for extrapolation.
   integer(I4P)                        :: bc_type              !< Boundary condition type.
   integer(I4P)                        :: crown                !< Crown counter.
   integer(I4P)                        :: fec                  !< Boundary fec (1 to 26).
   integer(I4P)                        :: fec_1_6              !< Boundary fec (1 to 6).

   associate(local_map_bc_crown=>self%field%maps%local_map_bc_crown, &
             nv=>self%nv, ngc=>self%ngc, cv=>self%physics%eos(1)%cv, R=>self%physics%eos(1)%R, q_bc_vars=>self%bc%q)
   if (allocated(self%field%maps%local_map_bc_crown)) then
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
   !< Set initial conditions of field.
   class(nasto_cpu_object), intent(inout) :: self !< The equation.

   call self%ic%set_initial_conditions(physics=self%physics, field=self%field)
   endsubroutine set_initial_conditions

   subroutine update_ghost(self, q, step)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(nasto_cpu_object), intent(inout)        :: self            !< The equation.
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
   class(nasto_cpu_object), intent(inout) :: self                            !< The equation.
   real(R8P)                              :: umax                            !< Maximum speed of waves propagation.
   real(R8P)                              :: ss                              !< Speed of sound.
   real(R8P)                              :: dx_locale, dy_locale, dz_locale !< Local space steps.
   integer(I4P)                           :: b, i, j, k                      !< Counter.

   call self%compute_q_aux(q=self%field%q, q_aux=self%q_aux)
   self%time%dt = huge(1._R8P)
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, blocks_number=>self%blocks_number, mu=>self%physics%eos(1)%mu, &
             dx=>self%field%dxyz(1,:), dy=>self%field%dxyz(2,:), dz=>self%field%dxyz(3,:), q_aux=>self%q_aux)
   umax = 0._R8P
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               dx_locale = dx(b)*0.5_R8P
               dy_locale = dy(b)*0.5_R8P
               dz_locale = dz(b)*0.5_R8P
               ss = q_aux(9,i,j,k,b)
               umax = max(umax, (abs(q_aux(2,i,j,k,b)) + ss)/dx_locale + 2._R8P*mu/(q_aux(1,i,j,k,b))/dx_locale**2 + &
                                (abs(q_aux(3,i,j,k,b)) + ss)/dy_locale + 2._R8P*mu/(q_aux(1,i,j,k,b))/dy_locale**2 + &
                                (abs(q_aux(4,i,j,k,b)) + ss)/dz_locale + 2._R8P*mu/(q_aux(1,i,j,k,b))/dz_locale**2)
            enddo
         enddo
      enddo
   enddo
   endassociate
   self%time%dt = min(self%time%dt, self%time%CFL / umax)
   call MPI_ALLREDUCE(MPI_IN_PLACE, self%time%dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, self%mpih%error)
   endsubroutine compute_dt

   subroutine compute_q_aux(self, q, q_aux)
   !< Compute auxiliary variables.
   class(nasto_cpu_object), intent(in)  :: self      !< The equation.
   real(R8P),               intent(in)  :: q(1:,         &
                                             1-self%ngc:,&
                                             1-self%ngc:,&
                                             1-self%ngc:,&
                                             1:)     !< Conservative variables.
   real(R8P),               intent(out) :: q_aux(1:,         &
                                                 1-self%ngc:,&
                                                 1-self%ngc:,&
                                                 1-self%ngc:,&
                                                 1:) !< Auxiliary variables.
   integer(I4P)                         :: b, i, j, k, s                          !< Counter.
   real(R8P)                            :: rho, uuu, vvv, www, rhe, rya, yya, tem !< State variables.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, ns=>self%ns, blocks_number=>self%blocks_number, &
             R=>self%physics%eos(1)%R, cv=>self%physics%eos(1)%cv, g=>self%physics%eos(1)%g, dha=>self%physics%eos(1)%dha)
   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               rho = q(1,i,j,k,b)
               uuu = q(2,i,j,k,b)/rho
               vvv = q(3,i,j,k,b)/rho
               www = q(4,i,j,k,b)/rho
               rhe = q(5,i,j,k,b)
               if (ns==2) then
                   rya = q(ns+4,i,j,k,b)
               else
                   rya = 0._R8P
               endif
               yya = rya/rho
               tem = ((rhe-rya*dha)/rho-0.5*(uuu**2+vvv**2+www**2))/cv

               q_aux(1,i,j,k,b) = rho           ! density
               q_aux(2,i,j,k,b) = uuu           ! velocity x
               q_aux(3,i,j,k,b) = vvv           ! velocity y
               q_aux(4,i,j,k,b) = www           ! velocity z
               q_aux(5,i,j,k,b) = yya           ! mass fraction
               q_aux(6,i,j,k,b) = tem           ! temperature
               q_aux(7,i,j,k,b) = R*rho*tem     ! pressure
               q_aux(8,i,j,k,b) = rhe/rho+R*tem ! entalpy
               q_aux(9,i,j,k,b) = sqrt(g*R*tem) ! sound speed
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine compute_q_aux

   subroutine compute_residuals(self)
   !< Compute residuals of equation.
   class(nasto_cpu_object), intent(inout) :: self         !< The equation.
   ! type(dim3)                             :: grid, tBlock !< CUDA grid and block.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk,                                                             &
             ngc=>self%ngc, nv=>self%nv, blocks_number=>self%blocks_number,                                     &
             dx=>self%field%dxyz(1,:), dy=>self%field%dxyz(2,:), dz=>self%field%dxyz(3,:),                      &
             q_aux=>self%q_aux, phi=>self%ib%phi, fl=>self%fl, flx=>self%flx, fly=>self%fly, flz=>self%flz,     &
             cell_scheme=>self%schemes%cell_scheme, ror_stats=>self%schemes%ror_stats,                          &
             fc_coeff=>self%schemes%fc_coeff,                                                                   &
             gminus_x=>self%gminus_x, gminus_y=>self%gminus_y, gminus_z=>self%gminus_z,                         &
             gplus_x=>self%gplus_x, gplus_y=>self%gplus_y, gplus_z=>self%gplus_z,                               &
             ror_schemes=>self%schemes%ror_schemes, ror_ivar=>self%schemes%ror_ivar,                            &
             ror_threshold=>self%schemes%ror_threshold, enable_ror_stats=>self%schemes%enable_ror_stats,        &
             lmax=>self%schemes%lmax, iweno=>self%schemes%iweno,                                                &
             cv=>self%physics%eos(1)%cv, g=>self%physics%eos(1)%g, R=>self%physics%eos(1)%R,                    &
             mu=>self%physics%eos(1)%mu, kd=>self%physics%eos(1)%kd, dha=>self%physics%eos(1)%dha)

   if (blocks_number > 0) then
      select case(self%schemes%fluxes_convective)
      case(SCHEME_FCONV_WENO_CENTRAL_2,SCHEME_FCONV_WENO_CENTRAL_4,SCHEME_FCONV_WENO_CENTRAL_6)
         ! to be implemented
      case(SCHEME_FCONV_WENO_UPWIND)
         call compute_flux_conv_x(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, &
                                  iweno=iweno, dha=dha, g=g, R=R, cv=cv,                            &
                                  ror_threshold=ror_threshold, enable_ror_stats=enable_ror_stats,   &
                                  cell_scheme=cell_scheme, ror_ivar=ror_ivar,                       &
                                  ror_schemes=ror_schemes, q_aux=q_aux, ror_stats=ror_stats,        &
                                  gplus=gplus_x, gminus=gminus_x, flx=flx)

         call compute_flux_conv_y(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, &
                                  iweno=iweno, dha=dha, g=g, R=R, cv=cv,                            &
                                  ror_threshold=ror_threshold, enable_ror_stats=enable_ror_stats,   &
                                  cell_scheme=cell_scheme, ror_ivar=ror_ivar,                       &
                                  ror_schemes=ror_schemes, q_aux=q_aux, ror_stats=ror_stats,        &
                                  gplus=gplus_y, gminus=gminus_y, fly=fly)

         call compute_flux_conv_z(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, &
                                  iweno=iweno, dha=dha, g=g, R=R, cv=cv,                            &
                                  ror_threshold=ror_threshold, enable_ror_stats=enable_ror_stats,   &
                                  cell_scheme=cell_scheme, ror_ivar=ror_ivar,                       &
                                  ror_schemes=ror_schemes, q_aux=q_aux, ror_stats=ror_stats,        &
                                  gplus=gplus_z, gminus=gminus_z, flz=flz)
      endselect
   endif

   if (mu > 0.) call compute_fluxes_diffusive(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, &
                                              mu=mu, kd=kd, q_aux=q_aux, dx=dx, dy=dy, dz=dz, flx=flx, fly=fly, flz=flz)

   call compute_fluxes_difference(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, ib_eps=1.e-12_R8P, &
                                  dx=dx, dy=dy, dz=dz, flx=flx, fly=fly, flz=flz, phi=phi, fl=fl)
   endassociate
   endsubroutine compute_residuals

   subroutine compute_rk_q(self, s)
   !< Compute RK approximation over q.
   class(nasto_cpu_object), intent(inout) :: self          !< The equation.
   integer(I4P),            intent(in)    :: s             !< Current RK stage.
   integer(I4P)                           :: all_solids    !< Last phi index, all solids summary.
   integer(I4P)                           :: i, j, k, b, v !< Counter.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, nv=>self%nv, blocks_number=>self%blocks_number, &
             dt=>self%time%dt, q=>self%field%q, q_old=>self%q_old, fl=>self%fl, phi=>self%ib%phi,   &
             solids_number=>self%ib%solids_number, ark=>self%schemes%ark(s), brk=>self%schemes%brk(s), crk=>self%schemes%crk(s))
   if (solids_number >0) then
      all_solids = ubound(phi, dim=1)
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do v=1, nv
                     if (phi(all_solids,i,j,k,b) < 0.) then
                        q(v,i,j,k,b) = ark * q_old(v,i,j,k,b) + brk * q(v,i,j,k,b) + dt * crk * fl(v,i,j,k,b)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
   else
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do v=1, nv
                     q(v,i,j,k,b) = ark * q_old(v,i,j,k,b) + brk * q(v,i,j,k,b) + dt * crk * fl(v,i,j,k,b)
                  enddo
               enddo
            enddo
         enddo
      enddo
   endif
   endassociate
   endsubroutine compute_rk_q

   subroutine integrate(self, t, do_ghost_syncro, residual)
   !< Perform one step integration.
   class(nasto_cpu_object), intent(inout)         :: self                  !< The equation.
   real(R8P),               intent(in)            :: t                     !< Time.
   logical,                 intent(in),  optional :: do_ghost_syncro       !< Flag to do syncrous ghost update.
   real(R8P),               intent(out), optional :: residual              !< Global residual.
   logical                                        :: do_ghost_syncro_      !< Flag to do syncrous ghost update, local var.
   integer(I4P)                                   :: s                     !< Counter.
   integer(I4P)                                   :: i_eikonal             !< Counter.
   integer(I4P), parameter                        :: n_eikonal=2           !< Counter.

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   self%q_old = self%field%q ! store previous conservative variables for RK integration
   do s=1, self%schemes%nrk
      if (self%ib%solids_number > 0) then ! integrate eikonal equation over q inside solids
         call self%update_ghost(q=self%field%q)
         do i_eikonal=1, n_eikonal
            call MPI_Barrier(MPI_COMM_WORLD, self%mpih%error)
            call self%integrate_eikonal_q
            call self%update_ghost(q=self%field%q)
         enddo
         call self%invert_eikonal_q
      endif
      call MPI_Barrier(MPI_COMM_WORLD, self%mpih%error)
      call self%compute_q_aux(q=self%field%q, q_aux=self%q_aux)
      call self%compute_residuals
      if (s==1) call self%save_residuals
      call self%compute_rk_q(s=s)
   enddo
   endsubroutine integrate

   subroutine simulate(self, filename)
   !< Perform the simulation.
   class(nasto_cpu_object), intent(inout) :: self             !< The equation.
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
      do i=1, 10
         call self%set_initial_conditions
         if (self%ib%solids_number > 0) call self%compute_phi()
         call self%amr_update()
      enddo
      call self%set_initial_conditions
      self%time%time = 0._R8P
      self%time%it = 0
   endif
   if (self%ib%solids_number > 0) call self%compute_phi()
   call self%amr_update()
   call self%save_simulation_data
   if (self%mpih%myrank==0) call self%io%open_file_residuals(nv=self%nv)

   ! integration
   call self%mpih%barrier(tictoc=.true., timing=timing(1), single=.true.)
   integration: do
      call self%mpih%barrier(tictoc=.true., timing=timing_step(1), single=.true.)
      self%time%it = self%time%it + 1

      if (self%io%save_memory_status) then
         call save_memory_cpu_status(file_name='memory_cpu-'//self%mpih%myrankstr//'.dat', tag=str(self%time%it,.true.))
      endif

      if (mod(self%time%it,self%amr%frequency)==0) then
         call self%mpih%barrier(tictoc=.true.)
         call self%amr_update()
         call self%mpih%barrier(tictoc=.true.)
      endif

      call self%compute_dt()
      if ((self%time%it_max <= 0).and.(self%time%time+self%time%dt > self%time%time_max)) &
         self%time%dt=self%time%time_max-self%time%time

      call self%integrate(t=self%time%time)

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
   endsubroutine simulate

   ! non TBP
   subroutine compute_flux_conv_x(blocks_number, ni, nj, nk, ngc, nv, iweno, dha, g, R, cv, &
                                  ror_threshold, enable_ror_stats, cell_scheme, ror_ivar,   &
                                  ror_schemes, q_aux, ror_stats, gplus, gminus, flx)
   !< Compute convective fluxes by means of upwind WENO reconstruction, x axis direction.
   integer(I4P), intent(in)    :: blocks_number, ni, nj, nk, ngc, nv, iweno
   real(R8P),    intent(in)    :: dha, g, R, cv
   real(R8P),    intent(in)    :: ror_threshold
   logical,      intent(in)    :: enable_ror_stats
   integer(I4P), intent(in)    :: cell_scheme(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P), intent(in)    :: ror_ivar(1:)
   integer(I4P), intent(in)    :: ror_schemes(1:)
   real(R8P),    intent(in)    :: q_aux(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   integer(I4P), intent(inout) :: ror_stats(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) ::  gplus(1:, 1:, 1:, 1:, 1:)
   real(R8P),    intent(inout) :: gminus(1:, 1:, 1:, 1:, 1:)
   real(R8P),    intent(inout) :: flx(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   integer(I4P)                :: b, i, j, k, l, ll, m, mm, v
   real(R8P)                   :: er(5,5), el(5,5), ev(5), evmax(5), ghat(5), gl(5), gr(5), fi(5), vi(5)
   real(R8P)                   :: uu, vv, ww, h, ya, qq, c, ci, b1, b2
   real(R8P)                   :: gc, wc
   integer(I4P)                :: wenorec_scheme, index_var
   logical                     :: ror_to_recompute

   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=0, ni ! loop on faces
               ! compute Roe average
               call compute_roe_average(q_aux=q_aux, dha=dha, g=g, ngc=ngc, b=b, i=i, j=j, k=k, ip=i+1, jp=j, kp=k, &
                                        uu=uu, vv=vv, ww=ww, h=h, ya=ya, qq=qq, c=c, ci=ci, b1=b1, b2=b2)
               ! compute right and left eigenvectors matrices (at Roe state)
               er(1,1)=1._R8P ; er(1,2)=uu-c   ; er(1,3)=vv     ; er(1,4)=ww     ; er(1,5)=h-uu*c
               er(2,1)=1._R8P ; er(2,2)=uu     ; er(2,3)=vv     ; er(2,4)=ww     ; er(2,5)=qq
               er(3,1)=1._R8P ; er(3,2)=uu+c   ; er(3,3)=vv     ; er(3,4)=ww     ; er(3,5)=h+uu*c
               er(4,1)=0._R8P ; er(4,2)=0._R8P ; er(4,3)=1._R8P ; er(4,4)=0._R8P ; er(4,5)=vv
               er(5,1)=0._R8P ; er(5,2)=0._R8P ; er(5,3)=0._R8P ; er(5,4)=1._R8P ; er(5,5)=ww

               el(1,1)= 0.5_R8P*(b1+uu*ci) ; el(1,2)=1._R8P-b1 ; el(1,3)= 0.5_R8P*(b1-uu*ci) ; el(1,4)=-vv    ; el(1,5)=-ww
               el(2,1)=-0.5_R8P*(b2*uu+ci) ; el(2,2)=b2*uu     ; el(2,3)=-0.5_R8P*(b2*uu-ci) ; el(2,4)=0._R8P ; el(2,5)=0._R8P
               el(3,1)=-0.5_R8P*(b2*vv   ) ; el(3,2)=b2*vv     ; el(3,3)=-0.5_R8P*(b2*vv   ) ; el(3,4)=1._R8P ; el(3,5)=0._R8P
               el(4,1)=-0.5_R8P*(b2*ww   ) ; el(4,2)=b2*ww     ; el(4,3)=-0.5_R8P*(b2*ww   ) ; el(4,4)=0._R8P ; el(4,5)=1._R8P
               el(5,1)= 0.5_R8P*b2         ; el(5,2)=-b2       ; el(5,3)= 0.5_R8P*b2         ; el(5,4)=0._R8P ; el(5,5)=0._R8P
               ! Find max eigenvalues on the stencil
               do m=1,nv  ! loop on characteristic fields
                  evmax(m) = -1._R8P
               enddo
               do l=1,2*iweno ! LLF
                  ll = i + l - iweno
                  uu = q_aux(2,ll,j,k,b)
                  c  = q_aux(9,ll,j,k,b)
                  ev(1) = abs(uu-c) ; ev(2) = abs(uu) ; ev(3) = abs(uu+c) ; ev(4) = ev(2) ; ev(5) = ev(2)
                  do m=1,nv
                     evmax(m) = max(ev(m),evmax(m))
                  enddo
               enddo
               ! Decompose fluxes as + and -
               do l=1,2*iweno ! loop over the stencil centered at face i
                  ll = i + l - iweno
                  vi(1) = q_aux(1,ll,j,k,b)
                  vi(2) = vi(1)*q_aux(2,ll,j,k,b)
                  vi(3) = vi(1)*q_aux(3,ll,j,k,b)
                  vi(4) = vi(1)*q_aux(4,ll,j,k,b)
                  vi(5) = vi(1)*(cv*q_aux(6,ll,j,k,b)+                                               &
                          0.5_R8P*(q_aux(2,ll,j,k,b)**2+q_aux(3,ll,j,k,b)**2+q_aux(4,ll,j,k,b)**2) + &
                          q_aux(5,ll,j,k,b)*dha)
                  fi(1) = vi(2)
                  fi(2) = fi(1) * q_aux(2,ll,j,k,b) + q_aux(7,ll,j,k,b)
                  fi(3) = fi(1) * q_aux(3,ll,j,k,b)
                  fi(4) = fi(1) * q_aux(4,ll,j,k,b)
                  fi(5) = fi(1) * vi(5) / vi(1) + q_aux(7,ll,j,k,b)*q_aux(2,ll,j,k,b)
                  do m=1,nv
                     wc = 0._R8P
                     gc = 0._R8P
                     do mm=1,nv
                        wc = wc + el(mm,m) * vi(mm)
                        gc = gc + el(mm,m) * fi(mm)
                     enddo
                     gplus (m,l,j,k,b) = 0.5_R8P * (gc + evmax(m) * wc)
                     gminus(m,l,j,k,b) = gc - gplus(m,l,j,k,b)
                  enddo
               enddo
               ! Reconstruction of the + and - fluxes
               wenorec_scheme = cell_scheme(1,i,j,k,b)
               call weno_reconstruction(nvar=nv, vp=gplus(1:,1:,j,k,b), vm=gminus(1:,1:,j,k,b), &
                                        vminus=gl, vplus=gr, iweno=iweno, wenorec_ord=wenorec_scheme)
               ror_x: do m = 2, size(ror_schemes)
                  ror_to_recompute = .false.
                  do mm = 1,size(ror_ivar)
                      index_var = ror_ivar(mm)
                      if ((abs(gl(index_var)-gplus(index_var,iweno,j,k,b))    > &
                           ror_threshold*abs(gplus(index_var,iweno,j,k,b))).or. &
                          (abs(gr(index_var)-gminus(index_var,iweno+1,j,k,b)) > &
                           ror_threshold*abs(gminus(index_var,iweno+1,j,k,b)))) then
                         ror_to_recompute = .true.
                      endif
                  enddo
                  if(ror_to_recompute) then
                     wenorec_scheme = ror_schemes(m)
                     call weno_reconstruction(nvar=nv, vp=gplus(1:,1:,j,k,b), vm=gminus(1:,1:,j,k,b), &
                                              vminus=gl, vplus=gr, iweno=iweno, wenorec_ord=wenorec_scheme)
                  else
                     exit ror_x
                  endif
               enddo ror_x
               if (enable_ror_stats) ror_stats(1,i,j,k,b) = wenorec_scheme
               ! Reassemble + and - characteristic fluxes
               do m=1,nv
                  ghat(m) = gl(m) + gr(m)
               enddo
               ! Return to conservative fluxes
               do m=1,nv
                  flx(m,i,j,k,b) = 0._R8P
                  do mm=1,nv
                     flx(m,i,j,k,b) = flx(m,i,j,k,b) + er(mm,m) * ghat(mm)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine compute_flux_conv_x

   subroutine compute_flux_conv_y(blocks_number, ni, nj, nk, ngc, nv, iweno, dha, g, R, cv, &
                                  ror_threshold, enable_ror_stats, cell_scheme, ror_ivar,   &
                                  ror_schemes, q_aux, ror_stats, gplus, gminus, fly)
   !< Compute convective fluxes by means of upwind WENO reconstruction, y axis direction.
   integer,      intent(in)    :: blocks_number, ni, nj, nk, ngc, nv, iweno
   real(R8P),    intent(in)    :: dha, g, R, cv
   real(R8P),    intent(in)    :: ror_threshold
   logical,      intent(in)    :: enable_ror_stats
   integer(I4P), intent(in)    :: cell_scheme(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P), intent(in)    :: ror_ivar(1:)
   integer(I4P), intent(in)    :: ror_schemes(1:)
   real(R8P),    intent(in)    :: q_aux(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   integer(I4P), intent(inout) :: ror_stats(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) ::  gplus(1:, 1:, 1:, 1:, 1:)
   real(R8P),    intent(inout) :: gminus(1:, 1:, 1:, 1:, 1:)
   real(R8P),    intent(inout) :: fly(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   integer                     :: b, i, j, k, l, ll, m, mm, v
   real(R8P)                   :: er(5,5), el(5,5), ev(5), evmax(5), ghat(5), gl(5), gr(5), fi(5), vi(5)
   real(R8P)                   :: uu, vv, ww, h, ya, qq, c, ci, b1, b2
   real(R8P)                   :: gc, wc
   integer                     :: wenorec_scheme, index_var
   logical                     :: ror_to_recompute

   do b=1, blocks_number
      do k=1,nk
         do j=0,nj ! loop on faces
            do i=1, ni
               ! Compute Roe average
               call compute_roe_average(q_aux=q_aux, dha=dha, g=g, ngc=ngc, b=b, i=i, j=j, k=k, ip=i, jp=j+1, kp=k, &
                                        uu=uu, vv=vv, ww=ww, h=h, ya=ya, qq=qq, c=c, ci=ci, b1=b1, b2=b2)
               ! Compute right and left eigenvectors matrices (at Roe state)
               er(1,1)=1._R8P ; er(1,2)=uu     ; er(1,3)=vv-c   ; er(1,4)=ww     ; er(1,5)=h-vv*c
               er(2,1)=1._R8P ; er(2,2)=uu     ; er(2,3)=vv     ; er(2,4)=ww     ; er(2,5)=qq
               er(3,1)=1._R8P ; er(3,2)=uu     ; er(3,3)=vv+c   ; er(3,4)=ww     ; er(3,5)=h+vv*c
               er(4,1)=0._R8P ; er(4,2)=1._R8P ; er(4,3)=0._R8P ; er(4,4)=0._R8P ; er(4,5)=ww
               er(5,1)=0._R8P ; er(5,2)=0._R8P ; er(5,3)=0._R8P ; er(5,4)=1._R8P ; er(5,5)=-uu

               el(1,1)= 0.5_R8P*(b1+vv*ci) ; el(1,2)=1._R8P-b1 ; el(1,3)=0.5_R8P*(b1-vv*ci)  ; el(1,4)=-ww    ; el(1,5)=uu
               el(2,1)=-0.5_R8P*(b2*uu)    ; el(2,2)=b2*uu     ; el(2,3)=-0.5_R8P*(b2*uu)    ; el(2,4)=0._R8P ; el(2,5)=-1._R8P
               el(3,1)=-0.5_R8P*(b2*vv+ci) ; el(3,2)=b2*vv     ; el(3,3)=-0.5_R8P*(b2*vv-ci) ; el(3,4)=0._R8P ; el(3,5)=0._R8P
               el(4,1)=-0.5_R8P*(b2*ww)    ; el(4,2)=b2*ww     ; el(4,3)=-0.5_R8P*(b2*ww)    ; el(4,4)=1._R8P ; el(4,5)=0._R8P
               el(5,1)= 0.5_R8P*b2         ; el(5,2)=-b2       ; el(5,3)=0.5_R8P*b2          ; el(5,4)=0._R8P ; el(5,5)=0._R8P
               ! Find max eigenvalues on the stencil
               do m=1,nv  ! loop on characteristic fields
                  evmax(m) = -1._R8P
               enddo
               do l=1,2*iweno ! LLF
                  ll = j + l - iweno
                  uu = q_aux(3,i,ll,k,b)
                  c  = q_aux(9,i,ll,k,b)
                  ev(1) = abs(uu-c) ; ev(2) = abs(uu) ; ev(3) = abs(uu+c) ; ev(4) = ev(2) ; ev(5) = ev(2) ;
                  do m=1,nv
                     evmax(m) = max(ev(m),evmax(m))
                  enddo
               enddo
               ! Decompose fluxes as + and -
               do l=1,2*iweno ! loop over the stencil centered at face i
                  ll = j + l - iweno
                  vi(1) = q_aux(1,i,ll,k,b)
                  vi(2) = vi(1)*q_aux(2,i,ll,k,b)
                  vi(3) = vi(1)*q_aux(3,i,ll,k,b)
                  vi(4) = vi(1)*q_aux(4,i,ll,k,b)
                  vi(5) = vi(1)*(cv*q_aux(6,i,ll,k,b)+                                           &
                      0.5_R8P*(q_aux(2,i,ll,k,b)**2+q_aux(3,i,ll,k,b)**2+q_aux(4,i,ll,k,b)**2) + &
                      q_aux(5,i,ll,k,b)*dha)
                  fi(1) = vi(3)
                  fi(2) = fi(1) * q_aux(2,i,ll,k,b)
                  fi(3) = fi(1) * q_aux(3,i,ll,k,b) + q_aux(7,i,ll,k,b)
                  fi(4) = fi(1) * q_aux(4,i,ll,k,b)
                  fi(5) = fi(1) * vi(5) / vi(1) + q_aux(7,i,ll,k,b)*q_aux(3,i,ll,k,b)
                  do m=1,nv
                     wc = 0._R8P
                     gc = 0._R8P
                     do mm=1,nv
                        wc = wc + el(mm,m) * vi(mm)
                        gc = gc + el(mm,m) * fi(mm)
                     enddo
                     gplus (m,l,i,k,b) = 0.5_R8P * (gc + evmax(m) * wc)
                     gminus(m,l,i,k,b) = gc - gplus(m,l,i,k,b)
                  enddo
               enddo
               ! Reconstruction of the + and - fluxes
               wenorec_scheme = cell_scheme(2,i,j,k,b)
               call weno_reconstruction(nvar=nv, vp=gplus(1:,1:,i,k,b), vm=gminus(1:,1:,i,k,b), &
                                        vminus=gl, vplus=gr, iweno=iweno, wenorec_ord=wenorec_scheme)
               ror_y: do m = 2, size(ror_schemes)
                  ror_to_recompute = .false.
                  do mm = 1,size(ror_ivar)
                      index_var = ror_ivar(mm)
                      if ((abs(gl(index_var)-gplus(index_var,iweno,i,k,b))    > &
                           ror_threshold*abs(gplus(index_var,iweno,i,k,b))).or. &
                          (abs(gr(index_var)-gminus(index_var,iweno+1,i,k,b)) > &
                           ror_threshold*abs(gminus(index_var,iweno+1,i,k,b)))) then
                         ror_to_recompute = .true.
                      endif
                  enddo
                  if (ror_to_recompute) then
                     wenorec_scheme = ror_schemes(m)
                     call weno_reconstruction(nvar=nv, vp=gplus(1:,1:,i,k,b), vm=gminus(1:,1:,i,k,b), &
                                              vminus=gl, vplus=gr, iweno=iweno, wenorec_ord=wenorec_scheme)
                  else
                     exit ror_y
                  endif
               enddo ror_y
               if (enable_ror_stats) ror_stats(2,i,j,k,b) = wenorec_scheme
               ! Reassemble + and - characteristic fluxes
               do m=1,nv
                  ghat(m) = gl(m) + gr(m)
               enddo
               ! Return to conservative fluxes
               do m=1,nv
                  fly(m,i,j,k,b) = 0._R8P
                  do mm=1,nv
                     fly(m,i,j,k,b) = fly(m,i,j,k,b) + er(mm,m) * ghat(mm)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine compute_flux_conv_y

   subroutine compute_flux_conv_z(blocks_number, ni, nj, nk, ngc, nv, iweno, dha, g, R, cv, &
                                  ror_threshold, enable_ror_stats, cell_scheme, ror_ivar,   &
                                  ror_schemes, q_aux, ror_stats, gplus, gminus, flz)
   !< Compute convective fluxes by means of upwind WENO reconstruction, z axis direction.
   integer,      intent(in)    :: blocks_number, ni, nj, nk, ngc, nv, iweno
   real(R8P),    intent(in)    :: dha, g, R, cv
   real(R8P),    intent(in)    :: ror_threshold
   logical,      intent(in)    :: enable_ror_stats
   integer(I4P), intent(in)    :: cell_scheme(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P), intent(in)    :: ror_ivar(1:)
   integer(I4P), intent(in)    :: ror_schemes(1:)
   real(R8P),    intent(in)    :: q_aux(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   integer(I4P), intent(inout) :: ror_stats(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout) ::  gplus(1:, 1:, 1:, 1:, 1:)
   real(R8P),    intent(inout) :: gminus(1:, 1:, 1:, 1:, 1:)
   real(R8P),    intent(inout) :: flz(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   integer                     :: b, i, j, k, l, ll, m, mm, v
   real(R8P)                   :: er(5,5), el(5,5), ev(5), evmax(5), ghat(5), gl(5), gr(5), fi(5), vi(5)
   real(R8P)                   :: uu, vv, ww, h, ya, qq, c, ci, b1, b2
   real(R8P)                   :: gc, wc
   integer                     :: wenorec_scheme, index_var
   logical                     :: ror_to_recompute

   do b=1, blocks_number
      do k=0,nk ! loop on faces
         do j=1,nj
            do i=1,ni
               ! Compute Roe average
               call compute_roe_average(q_aux=q_aux, dha=dha, g=g, ngc=ngc, b=b, i=i, j=j, k=k, ip=i, jp=j, kp=k+1, &
                                        uu=uu, vv=vv, ww=ww, h=h, ya=ya, qq=qq, c=c, ci=ci, b1=b1, b2=b2)
               ! Compute right and left eigenvectors matrices (at Roe state)
               er(1,1)=1._R8P ; er(1,2)=uu     ; er(1,3)=vv     ; er(1,4)=ww-c   ; er(1,5)=h-ww*c
               er(2,1)=1._R8P ; er(2,2)=uu     ; er(2,3)=vv     ; er(2,4)=ww     ; er(2,5)=qq
               er(3,1)=1._R8P ; er(3,2)=uu     ; er(3,3)=vv     ; er(3,4)=ww+c   ; er(3,5)=h+ww*c
               er(4,1)=0._R8P ; er(4,2)=1._R8P ; er(4,3)=0._R8P ; er(4,4)=0._R8P ; er(4,5)=uu
               er(5,1)=0._R8P ; er(5,2)=0._R8P ; er(5,3)=1._R8P ; er(5,4)=0._R8P ; er(5,5)=vv

               el(1,1)=0.5_R8P*(b1+ww*ci)  ; el(1,2)=1._R8P-b1 ; el(1,3)=0.5_R8P*(b1-ww*ci)  ; el(1,4)=-uu    ; el(1,5)=-vv
               el(2,1)=-0.5_R8P*(b2*uu)    ; el(2,2)=b2*uu     ; el(2,3)=-0.5_R8P*(b2*uu)    ; el(2,4)=1._R8P ; el(2,5)=0._R8P
               el(3,1)=-0.5_R8P*(b2*vv)    ; el(3,2)=b2*vv     ; el(3,3)=-0.5_R8P*(b2*vv)    ; el(3,4)=0._R8P ; el(3,5)=1._R8P
               el(4,1)=-0.5_R8P*(b2*ww+ci) ; el(4,2)=b2*ww     ; el(4,3)=-0.5_R8P*(b2*ww-ci) ; el(4,4)=0._R8P ; el(4,5)=0._R8P
               el(5,1)=0.5_R8P*b2          ; el(5,2)=-b2       ; el(5,3)=0.5_R8P*b2          ; el(5,4)=0._R8P ; el(5,5)=0._R8P
               ! Find max eigenvalues on the stencil
               do m=1,nv  ! loop on characteristic fields
                  evmax(m) = -1._R8P
               enddo
               do l=1,2*iweno ! LLF
                  ll = k + l - iweno
                  uu = q_aux(4,i,j,ll,b)
                  c  = q_aux(9,i,j,ll,b)
                  ev(1) = abs(uu-c) ; ev(2) = abs(uu) ; ev(3) = abs(uu+c) ; ev(4) = ev(2) ; ev(5) = ev(2)
                  do m=1,nv
                     evmax(m) = max(ev(m),evmax(m))
                  enddo
               enddo
               ! Decompose fluxes as + and -
               do l=1,2*iweno ! loop over the stencil centered at face i
                  ll = k + l - iweno
                  vi(1) = q_aux(1,i,j,ll,b)
                  vi(2) = vi(1)*q_aux(2,i,j,ll,b)
                  vi(3) = vi(1)*q_aux(3,i,j,ll,b)
                  vi(4) = vi(1)*q_aux(4,i,j,ll,b)
                  vi(5) = vi(1)*(cv*q_aux(6,i,j,ll,b)+                                           &
                      0.5_R8P*(q_aux(2,i,j,ll,b)**2+q_aux(3,i,j,ll,b)**2+q_aux(4,i,j,ll,b)**2) + &
                      q_aux(5,i,j,ll,b)*dha)
                  fi(1) = vi(4)
                  fi(2) = fi(1) * q_aux(2,i,j,ll,b)
                  fi(3) = fi(1) * q_aux(3,i,j,ll,b)
                  fi(4) = fi(1) * q_aux(4,i,j,ll,b) + q_aux(7,i,j,ll,b)
                  fi(5) = fi(1) * vi(5) / vi(1) + q_aux(7,i,j,ll,b)*q_aux(4,i,j,ll,b)
                  do m=1,nv
                     wc = 0._R8P
                     gc = 0._R8P
                     do mm=1,nv
                        wc = wc + el(mm,m) * vi(mm)
                        gc = gc + el(mm,m) * fi(mm)
                     enddo
                     gplus (m,l,i,j,b) = 0.5_R8P * (gc + evmax(m) * wc)
                     gminus(m,l,i,j,b) = gc - gplus(m,l,i,j,b)
                  enddo
               enddo
               ! Reconstruction of the + and - fluxes
               wenorec_scheme = cell_scheme(3,i,j,k,b)
               call weno_reconstruction(nvar=nv, vp=gplus(1:,1:,i,j,b), vm=gminus(1:,1:,i,j,b), &
                                        vminus=gl, vplus=gr, iweno=iweno, wenorec_ord=wenorec_scheme)
               ror_z: do m = 2, size(ror_schemes)
                  ror_to_recompute = .false.
                  do mm = 1,size(ror_ivar)
                      index_var = ror_ivar(mm)
                      if ((abs(gl(index_var)-gplus(index_var,iweno,i,j,b))    > &
                           ror_threshold*abs(gplus(index_var,iweno,i,j,b))).or. &
                          (abs(gr(index_var)-gminus(index_var,iweno+1,i,j,b)) > &
                           ror_threshold*abs(gminus(index_var,iweno+1,i,j,b)))) then
                         ror_to_recompute = .true.
                      endif
                  enddo
                  if(ror_to_recompute) then
                     wenorec_scheme = ror_schemes(m)
                     call weno_reconstruction(nvar=nv, vp=gplus(1:,1:,i,j,b), vm=gminus(1:,1:,i,j,b), &
                                              vminus=gl, vplus=gr, iweno=iweno, wenorec_ord=wenorec_scheme)
                  else
                     exit ror_z
                  endif
               enddo ror_z
               if (enable_ror_stats) ror_stats(3,i,j,k,b) = wenorec_scheme
               ! Reassemble + and - characteristic fluxes
               do m=1,nv
                  ghat(m) = gl(m) + gr(m)
               enddo
               ! Return to conservative fluxes
               do m=1,nv
                  flz(m,i,j,k,b) = 0._R8P
                  do mm=1,nv
                     flz(m,i,j,k,b) = flz(m,i,j,k,b) + er(mm,m) * ghat(mm)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine compute_flux_conv_z

   subroutine compute_fluxes_difference(blocks_number, ni, nj, nk, ngc, nv, ib_eps, dx, dy, dz, flx, fly, flz, phi, fl)
   !< Compute fluxes difference.
   integer(I4P), intent(in)    :: blocks_number                   !< Number of blocks.
   integer(I4P), intent(in)    :: ni                              !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                              !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                              !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                             !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                              !< Number of conservative varibales.
   real(R8P),    intent(in)    :: ib_eps                          !< Tolerance IB delta ratio.
   real(R8P),    intent(in)    :: dx(1:), dy(1:), dz(1:)          !< Space steps.
   real(R8P),    intent(in)    :: flx(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< X direction fluxes.
   real(R8P),    intent(in)    :: fly(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Y direction fluxes.
   real(R8P),    intent(in)    :: flz(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Z direction fluxes.
   real(R8P),    intent(in)    :: phi(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< IB distance function.
   real(R8P),    intent(inout) ::  fl(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Fluxes differences.
   real(R8P)                   :: delta_x, delta_y, delta_z       !< Space steps.
   real(R8P)                   :: dx_locale, dy_locale, dz_locale !< Local space steps.
   integer(I4P)                :: b, i, j, k, v                   !< Counter.
   integer(I4P)                :: all_solids                      !< Last phi index, all solids summary.

   all_solids = ubound(phi, dim=1)
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
         fl(v,i,j,k,b) = - (flx(v,i,j,k,b)-flx(v,i-1,j,k,b))/dx_locale &
                         - (fly(v,i,j,k,b)-fly(v,i,j-1,k,b))/dy_locale &
                         - (flz(v,i,j,k,b)-flz(v,i,j,k-1,b))/dz_locale
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_fluxes_difference

   subroutine compute_fluxes_diffusive(blocks_number, ni, nj, nk, ngc, nv, mu, kd, q_aux, dx, dy, dz, flx, fly, flz)
   !< Compute diffusive fluxes.
   integer(I4P), intent(in)    :: blocks_number, ni, nj, nk, ngc, nv
   real(R8P),    intent(in)    :: mu, kd
   real(R8P),    intent(in)    :: dx(1:), dy(1:), dz(1:)
   real(R8P),    intent(in)    :: q_aux(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P),    intent(inout) ::   flx(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P),    intent(inout) ::   fly(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P),    intent(inout) ::   flz(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   integer(I4P)                :: b, i, j, k, v
   real(R8P)                   :: du_dx, dv_dx, dw_dx, du_dy, dv_dy, dw_dy, du_dz, dv_dz, dw_dz
   real(R8P)                   :: dx_locale, dy_locale, dz_locale
   real(R8P)                   :: delta_x, delta_y, delta_z
   real(R8P)                   :: sigq, sigl
   real(R8P)                   :: tau_1_1, tau_2_1, tau_3_1, dT_dx
   real(R8P)                   :: tau_1_2, tau_2_2, tau_3_2, dT_dy
   real(R8P)                   :: tau_1_3, tau_2_3, tau_3_3, dT_dz
   real(R8P)                   :: vel_u, vel_v, vel_w
   real(R8P), parameter        :: ib_eps=1.e-12_R8P

   do b=1,blocks_number
      do k=1,nk
         do j=1,nj
            do i=0,ni ! loop on faces
                du_dx = (q_aux(2,i+1,j,k,b)-q_aux(2,i,j,k,b))/dx(b)
                dv_dx = (q_aux(3,i+1,j,k,b)-q_aux(3,i,j,k,b))/dx(b)
                dw_dx = (q_aux(4,i+1,j,k,b)-q_aux(4,i,j,k,b))/dx(b)

                du_dy = (q_aux(2,i+1,j+1,k,b) - q_aux(2,i+1,j-1,k,b)+ &
                         q_aux(2,i,  j+1,k,b) - q_aux(2,i,  j-1,k,b))*0.25_R8P/dy(b)
                dv_dy = (q_aux(3,i+1,j+1,k,b) - q_aux(3,i+1,j-1,k,b)+ &
                         q_aux(3,i,  j+1,k,b) - q_aux(3,i,  j-1,k,b))*0.25_R8P/dy(b)

                du_dz = (q_aux(2,i+1,j,k+1,b) - q_aux(2,i+1,j,k-1,b)+ &
                         q_aux(2,i,  j,k+1,b) - q_aux(2,i,  j,k-1,b))*0.25_R8P/dz(b)
                dw_dz = (q_aux(4,i+1,j,k+1,b) - q_aux(4,i+1,j,k-1,b)+ &
                         q_aux(4,i,  j,k+1,b) - q_aux(4,i,  j,k-1,b))*0.25_R8P/dz(b)

                vel_u = 0.5*(q_aux(2,i,j,k,b) + q_aux(2,i+1,j,k,b))
                vel_v = 0.5*(q_aux(3,i,j,k,b) + q_aux(3,i+1,j,k,b))
                vel_w = 0.5*(q_aux(4,i,j,k,b) + q_aux(4,i+1,j,k,b))

                tau_1_1 = 2.0*mu*(du_dx-1./3.*(du_dx+dv_dy+dw_dz))
                tau_2_1 = mu*(dv_dx+du_dy)
                tau_3_1 = mu*(dw_dx+du_dz)

                dT_dx = (q_aux(6,i+1,j,k,b)-q_aux(6,i,j,k,b))/dx(b)

                sigq = kd*dT_dx
                sigl = vel_u*tau_1_1+vel_v*tau_2_1+vel_w*tau_3_1

                flx(2,i,j,k,b) = flx(2,i,j,k,b) - tau_1_1
                flx(3,i,j,k,b) = flx(3,i,j,k,b) - tau_2_1
                flx(4,i,j,k,b) = flx(4,i,j,k,b) - tau_3_1
                flx(5,i,j,k,b) = flx(5,i,j,k,b) - sigq + sigl
            enddo
         enddo
      enddo
   enddo

   do b=1,blocks_number
      do k=1,nk
         do j=0,nj ! loop on faces
            do i=1,ni
                du_dy = (q_aux(2,i,j+1,k,b)-q_aux(2,i,j,k,b))/dy(b)
                dv_dy = (q_aux(3,i,j+1,k,b)-q_aux(3,i,j,k,b))/dy(b)
                dw_dy = (q_aux(4,i,j+1,k,b)-q_aux(4,i,j,k,b))/dy(b)

                du_dx = (q_aux(2,i+1,j+1,k,b) - q_aux(2,i-1,j+1,k,b)+ &
                         q_aux(2,i+1,j,  k,b) - q_aux(2,i-1,j,  k,b))*0.25_R8P/dx(b)
                dv_dx = (q_aux(3,i+1,j+1,k,b) - q_aux(3,i-1,j+1,k,b)+ &
                         q_aux(3,i+1,j,  k,b) - q_aux(3,i-1,j,  k,b))*0.25_R8P/dx(b)

                dv_dz = (q_aux(3,i+1,j,k+1,b) - q_aux(3,i-1,j,k+1,b)+ &
                         q_aux(3,i+1,j,k,  b) - q_aux(3,i-1,j,k,  b))*0.25_R8P/dz(b)
                dw_dz = (q_aux(4,i+1,j,k+1,b) - q_aux(4,i-1,j,k+1,b)+ &
                         q_aux(4,i+1,j,k,  b) - q_aux(4,i-1,j,k,  b))*0.25_R8P/dz(b)

                vel_u = 0.5*(q_aux(2,i,j,k,b) + q_aux(2,i,j+1,k,b))
                vel_v = 0.5*(q_aux(3,i,j,k,b) + q_aux(3,i,j+1,k,b))
                vel_w = 0.5*(q_aux(4,i,j,k,b) + q_aux(4,i,j+1,k,b))

                tau_1_2 = mu*(du_dy+dv_dx)
                tau_2_2 = 2.0*mu*(dv_dy-1./3.*(du_dx+dv_dy+dw_dz))
                tau_3_2 = mu*(dw_dy+dv_dz)

                dT_dy = (q_aux(6,i,j+1,k,b)-q_aux(6,i,j,k,b))/dy(b)

                sigq = kd*dT_dy
                sigl = vel_u*tau_1_2+vel_v*tau_2_2+vel_w*tau_3_2

                fly(2,i,j,k,b) = fly(2,i,j,k,b) - tau_1_2
                fly(3,i,j,k,b) = fly(3,i,j,k,b) - tau_2_2
                fly(4,i,j,k,b) = fly(4,i,j,k,b) - tau_3_2
                fly(5,i,j,k,b) = fly(5,i,j,k,b) - sigq + sigl
            enddo
         enddo
      enddo
   enddo

   do b=1,blocks_number
      do k=0,nk ! loop on faces
         do j=1,nj
            do i=1,ni
                du_dz = (q_aux(2,i,j,k+1,b)-q_aux(2,i,j,k,b))/dz(b)
                dv_dz = (q_aux(3,i,j,k+1,b)-q_aux(3,i,j,k,b))/dz(b)
                dw_dz = (q_aux(4,i,j,k+1,b)-q_aux(4,i,j,k,b))/dz(b)

                du_dx = (q_aux(2,i+1,j,k+1,b) - q_aux(2,i-1,j,k+1,b)+ &
                         q_aux(2,i+1,j,k,  b) - q_aux(2,i-1,j,k,  b))*0.25_R8P/dx(b)
                dw_dx = (q_aux(4,i+1,j,k+1,b) - q_aux(4,i-1,j,k+1,b)+ &
                         q_aux(4,i+1,j,k,  b) - q_aux(4,i-1,j,k,  b))*0.25_R8P/dx(b)

                dv_dy = (q_aux(3,i,j+1,k+1,b) - q_aux(3,i,j-1,k+1,b)+ &
                         q_aux(3,i,j+1,k,  b) - q_aux(3,i,j-1,k,  b))*0.25_R8P/dy(b)
                dw_dy = (q_aux(4,i,j+1,k+1,b) - q_aux(4,i,j-1,k+1,b)+ &
                         q_aux(4,i,j+1,k,  b) - q_aux(4,i,j-1,k,  b))*0.25_R8P/dy(b)

                vel_u = 0.5*(q_aux(2,i,j,k,b) + q_aux(2,i,j,k+1,b))
                vel_v = 0.5*(q_aux(3,i,j,k,b) + q_aux(3,i,j,k+1,b))
                vel_w = 0.5*(q_aux(4,i,j,k,b) + q_aux(4,i,j,k+1,b))

                tau_1_3 = mu*(du_dz+dw_dx)
                tau_2_3 = mu*(dv_dz+dw_dy)
                tau_3_3 = 2.0*mu*(dw_dz-1./3.*(du_dx+dv_dy+dw_dz))

                dT_dz = (q_aux(6,i,j,k+1,b)-q_aux(6,i,j,k,b))/dz(b)

                sigq = kd*dT_dz
                sigl = vel_u*tau_1_3+vel_v*tau_2_3+vel_w*tau_3_3

                flz(2,i,j,k,b) = flz(2,i,j,k,b) - tau_1_3
                flz(3,i,j,k,b) = flz(3,i,j,k,b) - tau_2_3
                flz(4,i,j,k,b) = flz(4,i,j,k,b) - tau_3_3
                flz(5,i,j,k,b) = flz(5,i,j,k,b) - sigq + sigl
            enddo
         enddo
      enddo
   enddo
   endsubroutine compute_fluxes_diffusive

   subroutine compute_roe_average(q_aux, dha, g, ngc, b, i, j, k, ip, jp, kp, uu, vv, ww, h, ya, qq, c, ci, b1, b2)
   !< Compute Roe averaged quantities.
   real(R8P),    intent(in)  :: q_aux(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P),    intent(in)  :: dha, g
   integer(I4P), intent(in)  :: ngc, b, i, j, k, ip, jp, kp
   real(R8P),    intent(out) :: uu, vv, ww, h, ya, qq, c, ci, b1, b2
   real(R8P)                 :: ri, up, vp, wp, hp, yap, r, rp1, cc
   ! Left state (node i)
   ri        =  1._R8P/q_aux(1,i,j,k,b)
   uu        =  q_aux(2,i,j,k,b)
   vv        =  q_aux(3,i,j,k,b)
   ww        =  q_aux(4,i,j,k,b)
   h         =  q_aux(8,i,j,k,b)
   ya        =  q_aux(5,i,j,k,b)
   ! Right state (node i+1)
   up        =  q_aux(2,ip,jp,kp,b)
   vp        =  q_aux(3,ip,jp,kp,b)
   wp        =  q_aux(4,ip,jp,kp,b)
   hp        =  q_aux(8,ip,jp,kp,b)
   yap       =  q_aux(5,ip,jp,kp,b)
   ! Average state
   r         =  sqrt(q_aux(1,ip,jp,kp,b)*ri)
   rp1       =  1._R8P/(r+1._R8P)
   uu        =  (r*up+uu)*rp1
   vv        =  (r*vp+vv)*rp1
   ww        =  (r*wp+ww)*rp1
   h         =  (r*hp+h)*rp1
   ya        =  (r*yap+ya)*rp1
   qq        =  0.5_R8P * (uu*uu+vv*vv+ww*ww)
   cc        =  (g-1._R8P) * (h - qq - ya*dha)
   !ERRATODIREIcc        =  g * (g-1._R8P) * (h - qq - ya*dha)
   c         =  sqrt(cc)
   ci        =  1._R8P/c
   b2        = (g-1)/cc  ! alias 1/(cp*theta)
   b1        = b2 * qq   ! alias q/(cp*theta)
   endsubroutine compute_roe_average

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
   integer(I4P), intent(in)  :: ivar                          !< Ghost cells number.
   real(R8P),    intent(out) :: gradient                      !< Maximum gradient of q.
   real(R8P)                 :: grad                          !< Current gradient of q.
   integer(I4P)              :: i, j, k                       !< Counter.
   real(R8P), parameter      :: tol=1.e-12                    !< Gradient denominator tolerance.

   gradient = 0._R8P
   do k=1, nk
      do j=1, nj
         do i=1, ni
            grad = sqrt(((q(ivar,i+1,j,k,b) - q(ivar,i-1,j,k,b))/(2*dx))**2 + &
                        ((q(ivar,i,j+1,k,b) - q(ivar,i,j-1,k,b))/(2*dy))**2 + &
                        ((q(ivar,i,j,k+1,b) - q(ivar,i,j,k-1,b))/(2*dz))**2)
            grad = grad/(abs(q(b,i,j,k,ivar))+tol)
            gradient = max(gradient, grad)
         enddo
      enddo
   enddo
   endsubroutine compute_q_gradient

   subroutine weno_reconstruction(nvar, vp, vm, vminus, vplus, iweno, wenorec_ord)
   !< Compute WENO reconstruction.
   integer, intent(in)                     :: nvar, iweno, wenorec_ord
   !real(R8P), dimension(nvar,2*iweno), intent(in)  :: vm,vp
   real(R8P), dimension(1:,1:)        :: vm,vp
   real(R8P), dimension(nvar), intent(out) :: vminus,vplus
   real(R8P), dimension(-1:4)              :: dwe               ! linear weights
   real(R8P), dimension(-1:4)              :: alfp,alfm         ! alpha_l
   real(R8P), dimension(-1:4)              :: alfp_map,alfm_map ! alpha_l
   real(R8P), dimension(-1:4)              :: betap,betam       ! beta_l
   real(R8P), dimension(-1:4)              :: omp,omm           ! WENO weights
   integer                                 :: r,i,j,k,l,m
   real(R8P)                               :: c0,c1,c2,c3,c4,d0,d1,d2,d3,d4,summ,sump
   real(R8P)                               :: x,y,y2

   if (wenorec_ord==1) then ! Godunov

       i = iweno ! index of intermediate node to perform reconstruction

       vminus(1:nvar) = vp(1:nvar,i)
       vplus (1:nvar) = vm(1:nvar,i+1)

   elseif (wenorec_ord==2) then ! WENO-3

       i = iweno ! index of intermediate node to perform reconstruction

       dwe(1)   = 2._R8P/3._R8P
       dwe(0)   = 1._R8P/3._R8P

       do m=1,nvar

           betap(0)  = (vp(m,i  )-vp(m,i-1))**2
           betap(1)  = (vp(m,i+1)-vp(m,i  ))**2
           betam(0)  = (vm(m,i+2)-vm(m,i+1))**2
           betam(1)  = (vm(m,i+1)-vm(m,i  ))**2

           sump = 0._R8P
           summ = 0._R8P
           do l=0,1
               alfp(l) = dwe(l)/(0.000001_R8P+betap(l))**2
               alfm(l) = dwe(l)/(0.000001_R8P+betam(l))**2
               sump = sump + alfp(l)
               summ = summ + alfm(l)
           enddo
           do l=0,1
               omp(l) = alfp(l)/sump
               omm(l) = alfm(l)/summ
           enddo

           vminus(m) = omp(0) *(-vp(m,i-1)+3*vp(m,i  )) + omp(1) *( vp(m,i  )+ vp(m,i+1))
           vplus(m)  = omm(0) *(-vm(m,i+2)+3*vm(m,i+1)) + omm(1) *( vm(m,i  )+ vm(m,i+1))

       enddo

       do m=1,nvar
           vminus(m) = 0.5_R8P*vminus(m)
           vplus(m)  = 0.5_R8P*vplus(m)
       enddo

     elseif (wenorec_ord==3) then ! WENO-5
!
      i = iweno ! index of intermediate node to perform reconstruction
!
      dwe( 0) = 1._R8P/10._R8P
      dwe( 1) = 6._R8P/10._R8P
      dwe( 2) = 3._R8P/10._R8P
!     JS
      d0 = 13._R8P/12._R8P
      d1 = 1._R8P/4._R8P
!     Weights for polynomial reconstructions
      c0 = 1._R8P/3._R8P
      c1 = 5._R8P/6._R8P
      c2 =-1._R8P/6._R8P
      c3 =-7._R8P/6._R8P
      c4 =11._R8P/6._R8P
!
      do m=1,nvar
!
       betap(2) = d0*(     vp(m,i)-2._R8P*vp(m,i+1)+vp(m,i+2))**2+d1*(3._R8P*vp(m,i)-4._R8P*vp(m,i+1)+vp(m,i+2))**2
       betap(1) = d0*(     vp(m,i-1)-2._R8P*vp(m,i)+vp(m,i+1))**2+d1*(     vp(m,i-1)-vp(m,i+1) )**2
       betap(0) = d0*(     vp(m,i)-2._R8P*vp(m,i-1)+vp(m,i-2))**2+d1*(3._R8P*vp(m,i)-4._R8P*vp(m,i-1)+vp(m,i-2))**2
!
       betam(2) = d0*(     vm(m,i+1)-2._R8P*vm(m,i)+vm(m,i-1))**2+d1*(3._R8P*vm(m,i+1)-4._R8P*vm(m,i)+vm(m,i-1))**2
       betam(1) = d0*(     vm(m,i+2)-2._R8P*vm(m,i+1)+vm(m,i))**2+d1*(     vm(m,i+2)-vm(m,i) )**2
       betam(0) = d0*(     vm(m,i+1)-2._R8P*vm(m,i+2)+vm(m,i+3))**2+d1*(3._R8P*vm(m,i+1)-4._R8P*vm(m,i+2)+vm(m,i+3))**2
!
       sump = 0._R8P
       summ = 0._R8P
       do l=0,2
        alfp(l) = dwe(  l)/(0.000001_R8P+betap(l))**2
        alfm(l) = dwe(  l)/(0.000001_R8P+betam(l))**2
        sump = sump + alfp(l)
        summ = summ + alfm(l)
       enddo
       do l=0,2
        omp(l) = alfp(l)/sump
        omm(l) = alfm(l)/summ
       enddo
!
       vminus(m)   = omp(2)*(c0*vp(m,i  )+c1*vp(m,i+1)+c2*vp(m,i+2)) + &
         & omp(1)*(c2*vp(m,i-1)+c1*vp(m,i  )+c0*vp(m,i+1)) + omp(0)*(c0*vp(m,i-2)+c3*vp(m,i-1)+c4*vp(m,i  ))
       vplus(m)   = omm(2)*(c0*vm(m,i+1)+c1*vm(m,i  )+c2*vm(m,i-1)) +  &
         & omm(1)*(c2*vm(m,i+2)+c1*vm(m,i+1)+c0*vm(m,i  )) + omm(0)*(c0*vm(m,i+3)+c3*vm(m,i+2)+c4*vm(m,i+1))
!
      enddo ! end of m-loop
!
   elseif (wenorec_ord==4) then ! WENO-7
!
      i = iweno ! index of intermediate node to perform reconstruction
!
      dwe( 0) = 1._R8P/35._R8P
      dwe( 1) = 12._R8P/35._R8P
      dwe( 2) = 18._R8P/35._R8P
      dwe( 3) = 4._R8P/35._R8P
!
!     JS weights
      d1 = 1._R8P/36._R8P
      d2 = 13._R8P/12._R8P
      d3 = 781._R8P/720._R8P
!
      do m=1,nvar
!
       betap(3)= d1*(-11*vp(m,  i)+18*vp(m,i+1)- 9*vp(m,i+2)+ 2*vp(m,i+3))**2+&
       &  d2*(  2*vp(m,  i)- 5*vp(m,i+1)+ 4*vp(m,i+2)-   vp(m,i+3))**2+ &
       & d3*(   -vp(m,  i)+ 3*vp(m,i+1)- 3*vp(m,i+2)+   vp(m,i+3))**2
       betap(2)= d1*(- 2*vp(m,i-1)- 3*vp(m,i  )+ 6*vp(m,i+1)-   vp(m,i+2))**2+&
       &  d2*(    vp(m,i-1)- 2*vp(m,i  )+   vp(m,i+1)             )**2+&
       &  d3*(   -vp(m,i-1)+ 3*vp(m,i  )- 3*vp(m,i+1)+   vp(m,i+2))**2
       betap(1)= d1*(    vp(m,i-2)- 6*vp(m,i-1)+ 3*vp(m,i  )+ 2*vp(m,i+1))**2+&
       &  d2*( vp(m,i-1)- 2*vp(m,i  )+   vp(m,i+1))**2+ &
       &  d3*(   -vp(m,i-2)+ 3*vp(m,i-1)- 3*vp(m,i  )+   vp(m,i+1))**2
       betap(0)= d1*(- 2*vp(m,i-3)+ 9*vp(m,i-2)-18*vp(m,i-1)+11*vp(m,i  ))**2+&
       &  d2*(-   vp(m,i-3)+ 4*vp(m,i-2)- 5*vp(m,i-1)+ 2*vp(m,i  ))**2+&
       &  d3*(   -vp(m,i-3)+ 3*vp(m,i-2)- 3*vp(m,i-1)+   vp(m,i  ))**2
!
       betam(3)= d1*(-11*vm(m,i+1)+18*vm(m,i  )- 9*vm(m,i-1)+ 2*vm(m,i-2))**2+&
       &  d2*(  2*vm(m,i+1)- 5*vm(m,i  )+ 4*vm(m,i-1)-   vm(m,i-2))**2+&
       &  d3*(   -vm(m,i+1)+ 3*vm(m,i  )- 3*vm(m,i-1)+   vm(m,i-2))**2
       betam(2)= d1*(- 2*vm(m,i+2)- 3*vm(m,i+1)+ 6*vm(m,i  )-   vm(m,i-1))**2+&
       &  d2*(    vm(m,i+2)- 2*vm(m,i+1)+   vm(m,i  )             )**2+&
       &  d3*(   -vm(m,i+2)+ 3*vm(m,i+1)- 3*vm(m,i  )+   vm(m,i-1))**2
       betam(1)= d1*(    vm(m,i+3)- 6*vm(m,i+2)+ 3*vm(m,i+1)+ 2*vm(m,i  ))**2+&
       &  d2*(                 vm(m,i+2)- 2*vm(m,i+1)+   vm(m,i  ))**2+&
       &  d3*(   -vm(m,i+3)+ 3*vm(m,i+2)- 3*vm(m,i+1)+   vm(m,i  ))**2
       betam(0)= d1*(- 2*vm(m,i+4)+ 9*vm(m,i+3)-18*vm(m,i+2)+11*vm(m,i+1))**2+&
       &  d2*(-   vm(m,i+4)+ 4*vm(m,i+3)- 5*vm(m,i+2)+ 2*vm(m,i+1))**2+&
       &  d3*(   -vm(m,i+4)+ 3*vm(m,i+3)- 3*vm(m,i+2)+   vm(m,i+1))**2
!
       sump = 0._R8P
       summ = 0._R8P
       do l=0,3
        alfp(l) = dwe(  l)/(0.000001_R8P+betap(l))**2
        alfm(l) = dwe(  l)/(0.000001_R8P+betam(l))**2
        sump = sump + alfp(l)
        summ = summ + alfm(l)
       enddo
       do l=0,3
        omp(l) = alfp(l)/sump
        omm(l) = alfm(l)/summ
       enddo
!
       vminus(m)   = omp(3)*( 6*vp(m,i  )+26*vp(m,i+1)-10*vp(m,i+2)+ 2*vp(m,i+3))+&
        omp(2)*(-2*vp(m,i-1)+14*vp(m,i  )+14*vp(m,i+1)- 2*vp(m,i+2))+&
        omp(1)*( 2*vp(m,i-2)-10*vp(m,i-1)+26*vp(m,i  )+ 6*vp(m,i+1))+&
        omp(0)*(-6*vp(m,i-3)+26*vp(m,i-2)-46*vp(m,i-1)+50*vp(m,i  ))
       vplus(m)   =  omm(3)*( 6*vm(m,i+1)+26*vm(m,i  )-10*vm(m,i-1)+ 2*vm(m,i-2))+&
        omm(2)*(-2*vm(m,i+2)+14*vm(m,i+1)+14*vm(m,i  )- 2*vm(m,i-1))+&
        omm(1)*( 2*vm(m,i+3)-10*vm(m,i+2)+26*vm(m,i+1)+ 6*vm(m,i  ))+&
        omm(0)*(-6*vm(m,i+4)+26*vm(m,i+3)-46*vm(m,i+2)+50*vm(m,i+1))
!
      enddo ! end of m-loop
!
      vminus = vminus/24._R8P
      vplus  = vplus /24._R8P
!
   else
      write(*,*) 'Error! WENO scheme not implemented'
      stop
   endif

   endsubroutine weno_reconstruction
endmodule adam_nasto_cpu_object
