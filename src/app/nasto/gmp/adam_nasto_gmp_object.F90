!< ADAM, Navier-Stokes equations system class definition, GPU (GMP) backend.
module adam_nasto_gmp_object
!< ADAM, Navier-Stokes equations system class definition, GPU (GMP) backend.

use adam_common_library
use adam_gmp_library
use adam_nasto_common_library
use adam_nasto_gmp_kernels
use penf
use mpi

implicit none
private
public :: nasto_gmp_object

type, extends(nasto_common_object) :: nasto_gmp_object
   !< Navier-Stokes equations system class definition, GPU (GMP) backend.
   ! ADAM library objects
   type(mpih_gmp_object)  :: mpih_gpu  !< MPI handler, GMP backend.
   type(field_gmp_object) :: field_gpu !< The field, GMP backend.
   type(ib_gmp_object)    :: ib_gpu    !< IB handler, GMP backend.
   type(rk_gmp_object)    :: rk_gpu    !< RK integrator, GMP backend.
   type(weno_gmp_object)  :: weno_gpu  !< WENO reconstructor, GMP backend.
   ! device data
   real(R8P), pointer, contiguous :: q_gpu(:,:,:,:,:)=>null()     !< Field cell centered variables.
   real(R8P), pointer, contiguous :: q_aux_gpu(:,:,:,:,:)=>null() !< Auxiliary cell centered variables.
   real(R8P), pointer, contiguous :: dq_gpu(:,:,:,:,:)=>null()    !< Eikonal right hand side.
   real(R8P), pointer, contiguous :: flx_gpu(:,:,:,:,:)=>null()   !< Fluxes along x.
   real(R8P), pointer, contiguous :: fly_gpu(:,:,:,:,:)=>null()   !< Fluxes along y.
   real(R8P), pointer, contiguous :: flz_gpu(:,:,:,:,:)=>null()   !< Fluxes along z.
   real(R8P), pointer, contiguous :: q_bc_vars_gpu(:,:)=>null()   !< Variables array for boundary conditions on GPU.
   contains
      ! auxiliary methods
      procedure, pass(self) :: allocate_gpu !< Allocate GPU data.
      procedure, pass(self) :: copy_cpu_gpu !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu !< Copy data from GPU to CPU.
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
      procedure, pass(self) :: save_hdf5            !< Save simulation data in HDF5 format.
      procedure, pass(self) :: save_residuals       !< Save residuals history.
      procedure, pass(self) :: save_restart_files   !< Save restart files.
      procedure, pass(self) :: save_simulation_data !< Save all simulation data.
      ! IC/BC
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions of equation.
      procedure, pass(self) :: set_initial_conditions  !< Set initial conditions of equation.
      procedure, pass(self) :: update_ghost            !< Update ghost cells and set boundary conditions.
      ! numerical methods
      procedure, pass(self) :: compute_dt          !< Compute time step.
      procedure, pass(self) :: compute_q_auxiliary !< Compute auxiliary variables.
      procedure, pass(self) :: compute_residuals   !< Compute residuals.
      procedure, pass(self) :: integrate           !< Perform one step integration.
      procedure, pass(self) :: simulate            !< Perform the simulation.
endtype nasto_gmp_object

contains
   ! auxiliary methods
   subroutine allocate_gpu(self, q_gpu)
   !< Allocate GPU data.
   class(nasto_gmp_object), intent(inout)         :: self      !< The equation.
   real(R8P),               intent(inout), target :: q_gpu(1:,         &
                                                           1-self%ngc:,&
                                                           1-self%ngc:,&
                                                           1-self%ngc:,&
                                                           1:) !< Conservative variables.
   character(:), allocatable                      :: msg_      !< Allocating message base.
   character(:), allocatable                      :: ms        !< Allocating message.
   integer(I4P)                                   :: sn1       !< Solids number + 1.

   call self%mpih_gpu%print_message('nasto_gmp_object%allocate_gpu start')
   self%q_gpu => q_gpu
   msg_ = self%mpih%myrankstr//'nasto_gmp_object%allocate_gpu '
   associate(nv=>self%nv, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, nrk=>self%rk%nrk, nb=>self%nb, nv_aux=>self%nv_aux,&
             mydev=>self%mpih_gpu%mydev)
   call assign_allocatable_gpu(lhs=self%q_bc_vars_gpu, rhs=self%bc%q, omp_dev=mydev, msg=msg_//' q_bc_vars_gpu ')
   ms = msg_//' q_aux_gpu '
   call alloc_var_gpu(var=self%q_aux_gpu,ulb=reshape([1,nb,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nv_aux],[2,5]),&
                      omp_dev=mydev,init_val=0._R8P,msg=ms)
   ms = msg_//' dq_gpu '
   call alloc_var_gpu(var=self%dq_gpu,   ulb=reshape([1,nb,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nv    ],[2,5]),&
                      omp_dev=mydev,init_val=0._R8P,msg=ms)
   ms = msg_//' flx_gpu '
   call alloc_var_gpu(var=self%flx_gpu,  ulb=reshape([1,nb,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nv    ],[2,5]),&
                      omp_dev=mydev,init_val=0._R8P,msg=ms)
   ms = msg_//' fly_gpu '
   call alloc_var_gpu(var=self%fly_gpu,  ulb=reshape([1,nb,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nv    ],[2,5]),&
                      omp_dev=mydev,init_val=0._R8P,msg=ms)
   ms = msg_//' flz_gpu '
   call alloc_var_gpu(var=self%flz_gpu,  ulb=reshape([1,nb,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nv    ],[2,5]),&
                      omp_dev=mydev,init_val=0._R8P,msg=ms)
   endassociate
   call self%mpih_gpu%print_message('nasto_gmp_object%allocate_gpu finish')
   endsubroutine allocate_gpu

   subroutine copy_cpu_gpu(self)
   !< Copy data from CPU to GPU.
   class(nasto_gmp_object), intent(inout) :: self !< The equation.

   call self%field_gpu%copy_transpose_cpu_gpu(nv=self%nv, q_cpu=self%field%q, q_gpu=self%q_gpu)
   call self%field_gpu%copy_cpu_gpu(verbose=.false.)
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self, compute_copy_q_aux, copy_phi)
   !< Copy data from GPU to CPU.
   class(nasto_gmp_object), intent(inout)        :: self               !< The equation.
   logical,                 intent(in), optional :: compute_copy_q_aux !< Flag to compute auxiliary variables.
   logical,                 intent(in), optional :: copy_phi           !< Copy also phi.

   call self%field_gpu%copy_transpose_gpu_cpu(nv=self%nv, q_gpu=self%q_gpu, q_cpu=self%field%q)
   if (present(compute_copy_q_aux)) then
      if (compute_copy_q_aux) then
         call self%compute_q_auxiliary(q_gpu=self%q_gpu, q_aux_gpu=self%q_aux_gpu)
         call self%field_gpu%copy_transpose_gpu_cpu(nv=self%nv_aux, q_gpu=self%q_aux_gpu, q_cpu=self%q_aux)
      endif
   endif
   if (present(copy_phi)) then
      if (copy_phi) then
         if (self%ib%solids_number>0) call self%field_gpu%copy_transpose_gpu_cpu(nv=self%ib%solids_number+1, &
                                                                                 q_gpu=self%ib_gpu%phi_gpu, q_cpu=self%ib%phi)
      endif
   endif
   endsubroutine copy_gpu_cpu

   subroutine initialize(self, filename)
   !< Initialize the equation.
   class(nasto_gmp_object), intent(inout) :: self     !< The equation.
   character(*),            intent(in)    :: filename !< Input file name.

   call self%mpih_gpu%initialize(do_mpi_init=.true., do_device_init=.true.)
   call self%mpih_gpu%print_message('nasto_gmp_object%initialize start')
   call self%initialize_common(filename=filename, memory_avail=self%mpih_gpu%memory_avail)
   call self%field_gpu%initialize(mpih=self%mpih_gpu, field=self%adam%field, nv_aux=self%nv_aux, verbose=.false.)
   call self%ib_gpu%initialize(mpih=self%mpih_gpu, ib=self%ib, field_gpu=self%field_gpu)
   call self%rk_gpu%initialize(mpih=self%mpih_gpu,rk=self%rk,nb=self%nb,ngc=self%ngc,ni=self%ni,nj=self%nj,nk=self%nk,nv=self%nv)
   call self%weno_gpu%initialize(mpih=self%mpih_gpu, weno=self%weno)
   call self%allocate_gpu(q_gpu=self%field_gpu%q_gpu)
   print '(A)', self%mpih_gpu%description()
   call self%mpih_gpu%print_message('nasto_gmp_object%initialize finish')
   endsubroutine initialize

   ! AMR methods
   subroutine amr_update(self)
   !< Do AMR update.
   class(nasto_gmp_object), intent(inout) :: self                !< The equation.
   logical                                :: is_grid_changed     !< Flag to check grid changes for each marker.
   logical                                :: is_grid_changed_all !< Flag to check grid changes for each iter.
   integer(I4P)                           :: i, i_marker         !< Counter.
   type(amr_marker_object)                :: amr_marker          !< Current amr marker.

   amr: do i=1, self%amr%iters
      is_grid_changed_all = .false.
      do i_marker=1, self%amr%markers_number
         amr_marker = self%amr%markers(i_marker)
         call self%update_ghost(q_gpu=self%q_gpu)
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
         call self%copy_gpu_cpu ! needed for adam%amr_update
         call self%adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false., is_grid_changed=is_grid_changed)
         call self%copy_cpu_gpu
         if (self%ib%solids_number > 0) call self%compute_phi()
         is_grid_changed_all = is_grid_changed_all.or.is_grid_changed
      enddo
      if (.not.is_grid_changed_all) then
          call self%mpih_gpu%print_message('AMR Grid stabilized after : '//trim(str(i))//' AMR iterations')
          exit amr
       elseif (i==self%amr%iters) then
          call self%mpih_gpu%print_message('AMR Grid is NOT stabilized after : '//trim(str(i))//' AMR iterations')
      endif
   enddo amr
   endsubroutine amr_update

   subroutine compute_phi(self)
   !< Compute phi, distance from IB solid.
   class(nasto_gmp_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: ib   !< Counter.

   associate(blocks_number=>self%blocks_number, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc,                             &
             x_cell_gpu=>self%field_gpu%x_cell_gpu, y_cell_gpu=>self%field_gpu%y_cell_gpu, z_cell_gpu=>self%field_gpu%z_cell_gpu, &
             phi_gpu=>self%ib_gpu%phi_gpu, solids_number=>self%ib%solids_number, definition=>self%ib%definition,                  &
             sphere=>self%ib%sphere, ib_reduction_extent=>self%weno%ib_reduction_extent,                                          &
             ib_reduced_order=>self%weno%ib_reduced_order, iweno=>self%weno%S, cell_scheme_gpu=>self%weno_gpu%cell_scheme_gpu)
   if (solids_number>0) then
      call self%mpih_gpu%print_message('compute IB distance start')
      do ib=1, solids_number
         ! compute phi
         select case(trim(adjustl(definition(ib))))
         case(trim(IB_ANALYTICAL_SPHERE))
            call compute_phi_analytical_sphere_gmp(ib=ib, ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number,    &
                                                   sphere=self%ib%sphere_to_array(ib=ib),                               &
                                                   x_cell_gpu=x_cell_gpu, y_cell_gpu=y_cell_gpu, z_cell_gpu=z_cell_gpu, &
                                                   phi_gpu=phi_gpu)
         endselect
         ! reduce local order of spatial operator close to solids if requested
         if (ib_reduction_extent > 0) call reduce_cell_order_phi_gmp(ib=ib,ni=ni,nj=nj,nk=nk,ngc=ngc,blocks_number=blocks_number, &
                                                                     iweno=iweno, ib_reduced_order=ib_reduced_order,              &
                                                                     ib_reduction_extent=ib_reduction_extent, phi_gpu=phi_gpu,    &
                                                                     cell_scheme_gpu=cell_scheme_gpu)
      enddo
      call compute_phi_all_solids_gmp(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number,phi_gpu=phi_gpu)
      call self%mpih_gpu%print_message('compute IB distance finish')
   endif
   endassociate
   endsubroutine compute_phi

   subroutine mark_by_geo(self, delta_fine, delta_coarse, threshold, do_init)
   !< Mark blocks to be refined/derefined by a geometric constrain.
   class(nasto_gmp_object), intent(inout)        :: self           !< The equation.
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

   subroutine mark_by_grad_var(self, grad_tol, delta_fine, delta_coarse, ivar, threshold, do_init)
   !< Mark blocks to be refined/derefined by a `grad(var)` value.
   class(nasto_gmp_object), intent(inout)        :: self           !< The equation.
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
   associate (blocks_number=>self%blocks_number, dxyz=>self%field%dxyz)
      call self%update_ghost(q_gpu=self%q_gpu)
      call self%compute_q_auxiliary(q_gpu=self%q_gpu, q_aux_gpu=self%q_aux_gpu)
      do b=1, blocks_number
         call self%field_gpu%compute_q_gradient(b=b, ivar=ivar_, q_gpu=self%q_aux_gpu, gradient=grad_var)
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

   subroutine move_phi(self, velocity)
   !< Move phi.
   class(nasto_gmp_object), intent(inout) :: self        !< The equation.
   real(R8P),               intent(in)    :: velocity(3) !< Velocity of the movement.

   if (self%ib%solids_number>0) then
      call self%mpih_gpu%print_message('move IB distance start')
      call move_phi_gmp(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                        velocity=velocity, phi_gpu=self%ib_gpu%phi_gpu, dphi_gpu=self%dq_gpu)
      call self%mpih_gpu%print_message('move IB distance finish')
   endif
   endsubroutine move_phi

   subroutine refine_uniform(self, refinement_levels)
   !< Refine all blocks uniformly.
   class(nasto_gmp_object), intent(inout) :: self              !< The equation.
   integer(I4P),            intent(in)    :: refinement_levels !< Number of refinement to be performed.
   integer(I4P)                           :: l                 !< Counter.

   do l=1, refinement_levels
      call self%adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
      call self%adam%amr_update(do_blocks_reorder=.false., do_mpi_redistribute=.true.)
   enddo
   endsubroutine

   ! IB methods
   subroutine integrate_eikonal(self, q_gpu)
   !< Integrate eikonal equation.
   class(nasto_gmp_object), intent(inout) :: self      !< The equation.
   real(R8P),               intent(inout) :: q_gpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1:) !< Conservative variables.
   integer(I4P)                           :: i_eikonal !< Counter.

   associate(blocks_number=>self%blocks_number, solids_number=>self%ib%solids_number, dq_gpu=>self%dq_gpu)
   if (blocks_number > 0) then
      if (solids_number > 0) then
         call self%update_ghost(q_gpu=q_gpu)
         do i_eikonal=1, self%ib%n_eikonal
            call self%mpih_gpu%barrier
            call self%ib_gpu%evolve_eikonal(dq_gpu=dq_gpu, q_gpu=q_gpu)
            call self%update_ghost(q_gpu=q_gpu)
         enddo
         call self%ib_gpu%invert_eikonal(q_gpu=q_gpu)
         call self%mpih_gpu%barrier
      endif
   endif
   endassociate
   endsubroutine integrate_eikonal

   ! IO methods
   subroutine load_restart_files(self, t, time)
   !< Save restart files.
   class(nasto_gmp_object), intent(inout) :: self !< The equation.
   integer(I4P),            intent(out)   :: t    !< Time iteration.
   real(R8P),               intent(out)   :: time !< Time.

   call self%adam%load_restart_files(basename=self%io%restart_basename, t=t, time=time)
   call self%adam%make_comm_local_maps_ghost_bc
   call self%copy_cpu_gpu
   endsubroutine load_restart_files

   subroutine save_hdf5(self, output_basename)
   !< Save simulation data in HDF5 format.
   class(nasto_gmp_object), intent(inout)        :: self             !< The equation.
   character(*),            intent(in), optional :: output_basename  !< Output basename.
   character(:), allocatable                     :: output_basename_ !< Output basename, local var.

   call self%mpih_gpu%barrier(tictoc=.true.)
   call self%mpih_gpu%print_message('save HDF5 files t: '//trim(str(self%time%it,.true.))//', time: '//&
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
   call self%mpih_gpu%barrier(tictoc=.true.)
   endsubroutine save_hdf5

   subroutine save_residuals(self)
   !< Save residuals history.
   class(nasto_gmp_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: v    !< Counter.

   if (self%time%is_to_save(it_save=self%io%residuals_save)) then
      call compute_normL2_residuals_gmp(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, nv=self%nv, &
                                        blocks_number=self%blocks_number, dq_gpu=self%dq_gpu, norm=self%field%residuals)
      do v=1, self%nv
         call MPI_ALLREDUCE(MPI_IN_PLACE, self%field%residuals(v), 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, self%mpih_gpu%error)
         self%field%residuals(v) = sqrt(self%field%residuals(v))
      enddo
      if (self%mpih_gpu%myrank==0) call self%io%save_residuals(it=self%time%it, time=self%time%time, &
                                                               blocks_number=self%blocks_number, residuals=self%field%residuals)
   endif
   endsubroutine save_residuals

   subroutine save_restart_files(self)
   !< Save restart files.
   class(nasto_gmp_object), intent(inout) :: self !< The equation.

   call self%mpih_gpu%barrier(tictoc=.true.)
   call self%mpih_gpu%print_message('save restart files t: '//trim(str(self%time%it,.true.))//', time: '//&
                                    trim(str(self%time%time,.true.)))
   call self%adam%save_restart_files(basename=self%io%restart_basename, t=self%time%it, time=self%time%time)
   call self%save_hdf5(output_basename=self%io%restart_basename)
   call self%mpih_gpu%barrier(tictoc=.true.)
   endsubroutine save_restart_files

   subroutine save_simulation_data(self)
   !< Save all simulation data.
   class(nasto_gmp_object), intent(inout) :: self !< The equation.

   if ((self%time%is_to_save(it_save=self%io%it_save)).or.      &
       (self%time%is_to_save(it_save=self%io%restart_save)).or. &
       (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max))) then
      call self%update_ghost(q_gpu=self%q_gpu)
      call self%copy_gpu_cpu(compute_copy_q_aux=.true., copy_phi=.true.)

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
   subroutine set_boundary_conditions(self, q_gpu)
   !< Set boundary conditions of equation.
   class(nasto_gmp_object), intent(in)    :: self                  !< The equation.
   real(R8P),               intent(inout) :: q_gpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,1:) !< Conservative variables.

   if (associated(self%field_gpu%maps%local_map_bc_crown_gpu)) &
      call set_bc_q_gpu_gmp(BC_EXTRAPOLATION=BC_EXTRAPOLATION, BC_INFLOW=BC_INFLOW,                       &
                            nv=self%nv, ngc=self%ngc, cv=self%physics%eos(1)%cv, R=self%physics%eos(1)%R, &
                            local_map_bc_gpu=self%field_gpu%maps%local_map_bc_crown_gpu,                  &
                            fec_1_6_array_gpu=self%field_gpu%fec_1_6_array_gpu,                           &
                            q_bc_vars_gpu=self%q_bc_vars_gpu,                                             &
                            q_gpu=q_gpu)
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self)
   !< Set initial conditions of field.
   class(nasto_gmp_object), intent(inout) :: self !< The equation.

   call self%ic%set_initial_conditions(physics=self%physics, field=self%field)
   call self%copy_cpu_gpu
   endsubroutine set_initial_conditions

   subroutine update_ghost(self, q_gpu, step)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(nasto_gmp_object), intent(inout)        :: self            !< The equation.
   real(R8P),               intent(inout)        :: q_gpu(1:,         &
                                                          1-self%ngc:,&
                                                          1-self%ngc:,&
                                                          1-self%ngc:,&
                                                          1:)       !< Conservative variables.
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

   if (do_local_update) call self%field_gpu%update_ghost_local(q_gpu=q_gpu)
                        call self%field_gpu%update_ghost_mpi(q_gpu=q_gpu, step=step)
   if (do_set_bc)       call self%set_boundary_conditions(q_gpu=q_gpu)
   endsubroutine update_ghost

   ! numerical methods
   subroutine compute_dt(self)
   !< Compute maximum time step accordingly to CFL stabilty criterion.
   class(nasto_gmp_object), intent(inout) :: self !< The equation.
   real(R8P)                              :: umax !< Maximum speed of waves propagation.
   integer(I4P)                           :: b    !< Counter.

   call self%compute_q_auxiliary(q_gpu=self%q_gpu, q_aux_gpu=self%q_aux_gpu)
   self%time%dt = huge(1._R8P)
   call compute_umax_gmp(ni=self%ni,nj=self%nj,nk=self%nk,ngc=self%ngc,blocks_number=self%blocks_number,mu=self%physics%eos(1)%mu,&
                         dx_gpu=self%field_gpu%dxyz_gpu(:,1),                                                                     &
                         dy_gpu=self%field_gpu%dxyz_gpu(:,2),                                                                     &
                         dz_gpu=self%field_gpu%dxyz_gpu(:,3), q_aux_gpu=self%q_aux_gpu,umax=umax)
   self%time%dt = min(self%time%dt, self%time%CFL / umax)
   call MPI_ALLREDUCE(MPI_IN_PLACE, self%time%dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, self%mpih_gpu%error)
   endsubroutine compute_dt

   subroutine compute_q_auxiliary(self, q_gpu, q_aux_gpu)
   !< Compute auxiliary variables.
   class(nasto_gmp_object), intent(in)    :: self          !< The equation.
   real(R8P),               intent(in)    :: q_gpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1:)     !< Conservative variables.
   real(R8P),               intent(inout) :: q_aux_gpu(1:,         &
                                                       1-self%ngc:,&
                                                       1-self%ngc:,&
                                                       1-self%ngc:,&
                                                       1:) !< Auxiliary variables.

   call compute_q_aux_gmp(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number, &
                          R=self%physics%eos(1)%R, cv=self%physics%eos(1)%cv, g=self%physics%eos(1)%g,        &
                          q_gpu=q_gpu, q_aux_gpu=q_aux_gpu)
   endsubroutine compute_q_auxiliary

   subroutine compute_residuals(self, q_gpu, dq_gpu)
   !< Compute residuals of equation.
   class(nasto_gmp_object), intent(inout) :: self       !< The equation.
   real(R8P),               intent(inout) :: q_gpu(1:,       &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1:)  !< Conservative variables.
   real(R8P),               intent(inout) :: dq_gpu(1:,         &
                                                    1-self%ngc:,&
                                                    1-self%ngc:,&
                                                    1-self%ngc:,&
                                                    1:) !< Residuals.

   call self%update_ghost(q_gpu=q_gpu)
   call self%integrate_eikonal(q_gpu=q_gpu)
   call self%compute_q_auxiliary(q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu)
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv, blocks_number=>self%blocks_number,      &
             dx_gpu=>self%field_gpu%dxyz_gpu(:,1),                                                                      &
             dy_gpu=>self%field_gpu%dxyz_gpu(:,2),                                                                      &
             dz_gpu=>self%field_gpu%dxyz_gpu(:,3),                                                                      &
             q_aux_gpu=>self%q_aux_gpu, phi_gpu=>self%ib_gpu%phi_gpu,                                                   &
             flx_gpu=>self%flx_gpu, fly_gpu=>self%fly_gpu, flz_gpu=>self%flz_gpu,                                       &
             weno_s=>self%weno%S,                                                                                       &
             weno_a_gpu=>self%weno_gpu%a_gpu, weno_p_gpu=>self%weno_gpu%p_gpu, weno_d_gpu=>self%weno_gpu%d_gpu,         &
             ror_number=>self%weno%ror_number,                                                                          &
             ror_schemes_gpu=>self%weno_gpu%ror_schemes_gpu, ror_ivar_gpu=>self%weno_gpu%ror_ivar_gpu,                  &
             ror_threshold=>self%weno%ror_threshold, enable_ror_stats=>self%weno%enable_ror_stats,                      &
             cell_scheme_gpu=>self%weno_gpu%cell_scheme_gpu, ror_stats_gpu=>self%weno_gpu%ror_stats_gpu,                &
             weno_zeps=>self%weno%zeps,                                                                                 &
             solids_number=>self%ib%solids_number,                                                                      &
             cv=>self%physics%eos(1)%cv, g=>self%physics%eos(1)%g, R=>self%physics%eos(1)%R,                            &
             mu=>self%physics%eos(1)%mu, kd=>self%physics%eos(1)%kd, mpih=>self%mpih_gpu)
   if (blocks_number > 0) then
      call compute_fluxes_convective_gmp(dir=1,                                                            &
                                         blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, &
                                         weno_s=weno_s, weno_a_gpu=weno_a_gpu, weno_p_gpu=weno_p_gpu,      &
                                         weno_d_gpu=weno_d_gpu, weno_zeps=weno_zeps,                       &
                                         ror_number=ror_number, ror_schemes_gpu=ror_schemes_gpu,           &
                                         ror_threshold=ror_threshold, ror_ivar_gpu=ror_ivar_gpu,           &
                                         ror_stats_gpu=ror_stats_gpu,                                      &
                                         g=g, q_aux_gpu=q_aux_gpu, fluxes_gpu=flx_gpu)
      call compute_fluxes_convective_gmp(dir=2,                                                            &
                                         blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, &
                                         weno_s=weno_s, weno_a_gpu=weno_a_gpu, weno_p_gpu=weno_p_gpu,      &
                                         weno_d_gpu=weno_d_gpu, weno_zeps=weno_zeps,                       &
                                         ror_number=ror_number, ror_schemes_gpu=ror_schemes_gpu,           &
                                         ror_threshold=ror_threshold, ror_ivar_gpu=ror_ivar_gpu,           &
                                         ror_stats_gpu=ror_stats_gpu,                                      &
                                         g=g, q_aux_gpu=q_aux_gpu, fluxes_gpu=fly_gpu)
      call compute_fluxes_convective_gmp(dir=3,                                                            &
                                         blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, &
                                         weno_s=weno_s, weno_a_gpu=weno_a_gpu, weno_p_gpu=weno_p_gpu,      &
                                         weno_d_gpu=weno_d_gpu, weno_zeps=weno_zeps,                       &
                                         ror_number=ror_number, ror_schemes_gpu=ror_schemes_gpu,           &
                                         ror_threshold=ror_threshold, ror_ivar_gpu=ror_ivar_gpu,           &
                                         ror_stats_gpu=ror_stats_gpu,                                      &
                                         g=g, q_aux_gpu=q_aux_gpu, fluxes_gpu=flz_gpu)
      if (mu > 0.) call compute_fluxes_diffusive_gmp(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, &
                                                     mu=mu, kd=kd, q_aux_gpu=q_aux_gpu,                                &
                                                     dx_gpu=dx_gpu, dy_gpu=dy_gpu, dz_gpu=dz_gpu,                      &
                                                     flx_gpu=flx_gpu, fly_gpu=fly_gpu, flz_gpu=flz_gpu)
      if (solids_number>0) then
         call compute_fluxes_difference_gmp(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, ib_eps=1.e-12_R8P, &
                                            dx_gpu=dx_gpu, dy_gpu=dy_gpu, dz_gpu=dz_gpu,                                         &
                                            flx_gpu=flx_gpu, fly_gpu=fly_gpu, flz_gpu=flz_gpu, phi_gpu=phi_gpu, dq_gpu=dq_gpu)
      else
         call compute_fluxes_difference_gmp(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, ib_eps=1.e-12_R8P, &
                                            dx_gpu=dx_gpu, dy_gpu=dy_gpu, dz_gpu=dz_gpu,                                         &
                                            flx_gpu=flx_gpu, fly_gpu=fly_gpu, flz_gpu=flz_gpu, dq_gpu=dq_gpu)
      endif
   endif
   endassociate
   endsubroutine compute_residuals

   subroutine integrate(self, do_ghost_syncro)
   !< Perform one step integration.
   class(nasto_gmp_object), intent(inout)         :: self                  !< The equation.
   logical,                 intent(in),  optional :: do_ghost_syncro       !< Flag to do syncrous ghost update.
   logical                                        :: do_ghost_syncro_      !< Flag to do syncrous ghost update, local var.
   integer(I4P)                                   :: s                     !< Counter.

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   call self%rk_gpu%initialize_stages(q_gpu=self%q_gpu)
   select case(self%rk%scheme)
   case(RK_1, RK_2, RK_3)
      ! low storage RK working on q_rk_gpu(:,:,:,:,:,1)/q_gpu as stages, update q_gpu in place
      do s=1, self%rk%nrk
         call self%compute_residuals(q_gpu=self%q_gpu, dq_gpu=self%dq_gpu)
         if (s==1) call self%save_residuals
         if (self%ib%solids_number>0) then
            call self%rk_gpu%compute_stage_ls(s=s,dt=self%time%dt,phi_gpu=self%ib_gpu%phi_gpu,dq_gpu=self%dq_gpu,q_gpu=self%q_gpu)
         else
            call self%rk_gpu%compute_stage_ls(s=s,dt=self%time%dt,dq_gpu=self%dq_gpu,q_gpu=self%q_gpu)
         endif
      enddo
   case(RK_SSP_22, RK_SSP_33, RK_SSP_54)
      ! RK working on q_rk_gpu as stages
      do s=1, self%rk%nrk
         if (self%ib%solids_number>0) then
            call self%rk_gpu%compute_stage(s=s, dt=self%time%dt, phi_gpu=self%ib_gpu%phi_gpu)
         else
            call self%rk_gpu%compute_stage(s=s, dt=self%time%dt)
         endif
         call self%compute_residuals(q_gpu=self%rk_gpu%q_rk_gpu(:,:,:,:,:,s), dq_gpu=self%dq_gpu)
         if (s==1) call self%save_residuals
         if (self%ib%solids_number>0) then
            call self%rk_gpu%assign_stage(s=s, q_gpu=self%dq_gpu, phi_gpu=self%ib_gpu%phi_gpu)
         else
            call self%rk_gpu%assign_stage(s=s, q_gpu=self%dq_gpu)
         endif
      enddo
      if (self%ib%solids_number>0) then
         call self%rk_gpu%update_q(dt=self%time%dt, phi_gpu=self%ib_gpu%phi_gpu, q_gpu=self%q_gpu)
      else
         call self%rk_gpu%update_q(dt=self%time%dt, q_gpu=self%q_gpu)
      endif
   endselect
   endsubroutine integrate

   subroutine simulate(self, filename)
   !< Perform the simulation.
   class(nasto_gmp_object), intent(inout) :: self             !< The equation.
   character(*),            intent(in)    :: filename         !< Input file name.
   real(R8P)                              :: timing(1:2)      !< Tic toc timing.
   real(R8P)                              :: timing_step(1:2) !< Tic toc timing.
   integer(I4P)                           :: i                !< Counter.

   ! initialization
   call self%initialize(filename=filename)
   if (self%io%restart) then
      call self%mpih_gpu%print_message('restart simulation from "'//trim(self%io%restart_basename)//'" files')
      call self%load_restart_files(t=self%time%it, time=self%time%time)
      call self%mpih_gpu%print_message('restart [t, time]: '//trim(str(self%time%it))//', '//trim(str(self%time%time)))
   else
      call self%mpih_gpu%print_message('impose initial conditions start')
      do i=1, self%ic%amr_iterations
         call self%mpih_gpu%print_message('  AMR/set IC iteration:'//trim(str(i,.true.)))
         call self%set_initial_conditions
         if (self%ib%solids_number > 0) call self%compute_phi
         call self%amr_update
      enddo
      call self%set_initial_conditions
      self%time%time = 0._R8P
      self%time%it = 0
      call self%mpih_gpu%print_message('impose initial conditions finish')
   endif
   if (self%ib%solids_number > 0) call self%compute_phi
   call self%amr_update
   call self%save_simulation_data
   if (self%mpih_gpu%myrank==0) call self%io%open_file_residuals(nv=self%nv)

   ! integration
   call self%mpih_gpu%barrier(tictoc=.true., timing=timing(1), single=.true.)
   integration: do
      call self%mpih_gpu%barrier(tictoc=.true., timing=timing_step(1), single=.true.)
      self%time%it = self%time%it + 1

      if (self%io%save_memory_status) then
         call save_memory_cpu_status(file_name='memory_cpu-'//self%mpih_gpu%myrankstr//'.dat', tag=str(self%time%it,.true.))
         call save_memory_gpu_status(file_name='memory_gpu-'//self%mpih_gpu%myrankstr//'.dat', tag=str(self%time%it,.true.))
      endif

      if (mod(self%time%it,self%amr%frequency)==0) then
         call self%mpih_gpu%barrier(tictoc=.true.)
         call self%amr_update
         call self%mpih_gpu%barrier(tictoc=.true.)
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

      call self%mpih_gpu%barrier(tictoc=.true., timing=timing_step(2), single=.true.)
   enddo integration
   call self%mpih_gpu%barrier(tictoc=.true., timing=timing(2), single=.true.)
   call self%save_simulation_data
   if (self%mpih_gpu%myrank==0) call self%io%close_file_residuals
   endsubroutine simulate
endmodule adam_nasto_gmp_object
