!< ADAM, PRISM (Plasma Research usIng Simulation Methods) equations system class definition, GPU (FNL) backend.

#include "fundal.H"

module adam_prism_fnl_object
!< ADAM, PRISM (Plasma Research usIng Simulation Methods) equations system class definition, GPU (FNL) backend.

use :: adam_prism_common_library
use :: adam_prism_fnl_library
use :: fundal, save_memory_status_gpu=>save_memory_status
use :: penf,   save_memory_status_cpu=>save_memory_status
use :: mpi

implicit none
private
public :: prism_fnl_object

type, extends(prism_common_object) :: prism_fnl_object
   !< PRISM equations system class definition, GPU (FNL) backend.
   ! ADAM library objects
   type(mpih_fnl_object)  :: mpih_gpu  !< MPI handler, FNL backend.
   type(field_fnl_object) :: field_gpu !< The field, FNL backend.
   type(rk_fnl_object)    :: rk_gpu    !< RK integrator, FNL backend.
   type(weno_fnl_object)  :: weno_gpu  !< WENO reconstructor, FNL backend.
   ! device data
   real(R8P), pointer :: q_gpu(  :,:,:,:,:)=>null() !< Field cell centered variables.
   real(R8P), pointer :: dq_gpu( :,:,:,:,:)=>null() !< Residuals right hand side.
   real(R8P), pointer :: flx_gpu(:,:,:,:,:)=>null() !< Fluxes along x.
   real(R8P), pointer :: fly_gpu(:,:,:,:,:)=>null() !< Fluxes along y.
   real(R8P), pointer :: flz_gpu(:,:,:,:,:)=>null() !< Fluxes along z.
   contains
      ! auxiliary methods
      procedure, pass(self) :: allocate_gpu !< Allocate GPU data.
      procedure, pass(self) :: copy_cpu_gpu !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu !< Copy data from GPU to CPU.
      procedure, pass(self) :: initialize   !< Initialize the equation.
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
      procedure, pass(self) :: compute_residuals   !< Compute residuals.
      procedure, pass(self) :: integrate           !< Perform one step integration.
      procedure, pass(self) :: simulate            !< Perform the simulation.
endtype prism_fnl_object

contains
   ! auxiliary methods
   subroutine allocate_gpu(self, q_gpu)
   !< Allocate GPU data.
   class(prism_fnl_object), intent(inout)         :: self      !< The equation.
   real(R8P),               intent(inout), target :: q_gpu(1:,         &
                                                           1-self%ngc:,&
                                                           1-self%ngc:,&
                                                           1-self%ngc:,&
                                                           1:) !< Conservative variables.
   integer(I4P)                                   :: ierr      !< Error status.

   call self%mpih%print_message('prism_fnl_object%allocate_gpu start')
   self%q_gpu => q_gpu ! q_gpu is allocated by field_gpu%initialize
   associate(nv=>self%nv, nb=>self%nb, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk)
   call dev_alloc(fptr_dev=self%dq_gpu, &
                  ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   call dev_alloc(fptr_dev=self%flx_gpu, &
                  ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   call dev_alloc(fptr_dev=self%fly_gpu, &
                  ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   call dev_alloc(fptr_dev=self%flz_gpu, &
                  ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   endassociate
   call self%mpih%print_message('prism_fnl_object%allocate_gpu finish')
   endsubroutine allocate_gpu

   subroutine copy_cpu_gpu(self)
   !< Copy data from CPU to GPU.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   call self%field_gpu%copy_transpose_cpu_gpu(nv=self%nv, q_cpu=self%field%q, q_gpu=self%q_gpu)
   call self%field_gpu%copy_cpu_gpu(verbose=.false.)
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self, compute_copy_q_aux, copy_phi)
   !< Copy data from GPU to CPU.
   class(prism_fnl_object), intent(inout)        :: self               !< The equation.
   logical,                 intent(in), optional :: compute_copy_q_aux !< Flag to compute auxiliary variables.
   logical,                 intent(in), optional :: copy_phi           !< Copy also phi.

   call self%field_gpu%copy_transpose_gpu_cpu(nv=self%nv, q_gpu=self%q_gpu, q_cpu=self%field%q)
   endsubroutine copy_gpu_cpu

   subroutine initialize(self, filename)
   !< Initialize the equation.
   class(prism_fnl_object), intent(inout) :: self     !< The equation.
   character(*),            intent(in)    :: filename !< Input file name.

   call self%mpih_gpu%initialize(do_mpi_init=.true., do_device_init=.true., verbose=.true.)
   call self%mpih_gpu%print_message('prism_fnl_object%initialize start')
   call self%initialize_common(field = self%adam%field, filename=filename, memory_avail=real(self%mpih_gpu%dev_memory_avail,R8P), &
                               verbose=.true.)
   ! call self%field_gpu%initialize(field=self%adam%field, nv_aux=self%nv_aux, verbose=.false.)
   ! call self%ib_gpu%initialize(ib=self%ib, field_gpu=self%field_gpu)
   ! call self%rk_gpu%initialize(rk=self%rk, nb=self%nb, ngc=self%ngc, ni=self%ni, nj=self%nj, nk=self%nk, nv=self%nv)
   ! call self%weno_gpu%initialize(weno=self%weno)
   call self%allocate_gpu(q_gpu=self%field_gpu%q_gpu)
   call self%mpih%print_message(self%mpih_gpu%description()//new_line('a')//'prism_fnl_object%initialize finish')
   endsubroutine initialize

   ! IO methods
   subroutine load_restart_files(self, t, time)
   !< Save restart files.
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P),            intent(out)   :: t    !< Time iteration.
   real(R8P),               intent(out)   :: time !< Time.

   call self%adam%load_restart_files(basename=self%io%restart_basename, t=t, time=time, q=self%q)
   call self%adam%make_comm_local_maps_ghost_bc
   call self%copy_cpu_gpu
   endsubroutine load_restart_files

   subroutine save_hdf5(self, output_basename)
   !< Save simulation data in HDF5 format.
   class(prism_fnl_object), intent(inout)        :: self             !< The equation.
   character(*),            intent(in), optional :: output_basename  !< Output basename.
   character(:), allocatable                     :: output_basename_ !< Output basename, local var.
   character(3), allocatable                     :: q_name(:)        !< Fields name.
   character(6), allocatable                     :: q_aux_name(:)    !< Auxiliary fields name.

   if ((.not.self%physics%D_divergence_cleaner).and.(.not.self%physics%B_divergence_cleaner)) then
      q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ']
   elseif ((self%physics%D_divergence_cleaner).and.(.not.self%physics%B_divergence_cleaner)) then
      q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ', 'phi']
   elseif ((self%physics%D_divergence_cleaner).and.(self%physics%B_divergence_cleaner)) then
      q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ', 'phi', 'psi']
   endif
   q_aux_name = ['DivD_d','DivB_d','DivD_f','DivB_f']

   call self%mpih_gpu%barrier(tictoc=.true.)
   call self%mpih_gpu%print_message('save HDF5 files t: '//trim(str(self%time%it,.true.))//', time: '//&
                                    trim(str(self%time%time,.true.)))
   output_basename_ = trim(self%io%output_basename)//'-'//trim(strz(self%time%it,9))
   if (present(output_basename)) output_basename_ = trim(output_basename)
   if (self%ib%solids_number>0) then
      call self%adam%save_hdf5(basename         = trim(output_basename_), &
                               q                = self%field%q,           &
                               q_aux            = self%field_div,         &
                               q_name           = q_name,                 &
                               q_aux_name       = q_aux_name,             &
                               with_cell_morton = .true.,                 &
                               phi=self%ib%phi)
   else
      call self%adam%save_hdf5(basename         = trim(output_basename_), &
                               q                = self%field%q,           &
                               q_aux            = self%field_div,         &
                               q_name           = q_name,                 &
                               q_aux_name       = q_aux_name,             &
                               with_cell_morton = .true.)
   endif
   call self%mpih_gpu%barrier(tictoc=.true.)
   endsubroutine save_hdf5

   subroutine save_residuals(self)
   !< Save residuals history.
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: v    !< Counter.

   if (self%time%is_to_save(it_save=self%io%residuals_save)) then
      call compute_normL2_residuals_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, nv=self%nv, &
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
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   call self%mpih_gpu%barrier(tictoc=.true.)
   call self%mpih_gpu%print_message('save restart files t: '//trim(str(self%time%it,.true.))//', time: '//&
                                    trim(str(self%time%time,.true.)))
   call self%adam%save_restart_files(basename=self%io%restart_basename, t=self%time%it, time=self%time%time, q=self%q)
   call self%save_hdf5(output_basename=self%io%restart_basename)
   call self%mpih_gpu%barrier(tictoc=.true.)
   endsubroutine save_restart_files

   subroutine save_simulation_data(self)
   !< Save all simulation data.
   class(prism_fnl_object), intent(inout) :: self      !< The equation.
   character(3), allocatable              :: q_name(:) !< Fields name.

   if ((.not.self%physics%D_divergence_cleaner).and.(.not.self%physics%B_divergence_cleaner)) then
      q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ']
   elseif ((self%physics%D_divergence_cleaner).and.(.not.self%physics%B_divergence_cleaner)) then
      q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ', 'phi']
   elseif ((self%physics%D_divergence_cleaner).and.(self%physics%B_divergence_cleaner)) then
      q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ', 'phi', 'psi']
   endif

   if ((self%time%is_to_save(it_save=self%io%it_save)).or.      &
       (self%time%is_to_save(it_save=self%io%restart_save)).or. &
       (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max))) then
      ! call self%update_ghost(q_gpu=self%q_gpu)
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
                                   q_name=q_name)
   endif
   endsubroutine save_simulation_data

   ! IC/BC
   subroutine set_boundary_conditions(self, q_gpu)
   !< Set boundary conditions of equation.
   class(prism_fnl_object), intent(in)    :: self                  !< The equation.
   real(R8P),               intent(inout) :: q_gpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,1:) !< Conservative variables.

   if (associated(self%field_gpu%maps%local_map_bc_crown_gpu)) &
      call set_bc_q_gpu_dev(BC_EXTRAPOLATION=BC_EXTRAPOLATION, BC_fWLAYER=BC_fWLAYER,     &
                            nv=self%nv, ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, &
                            D_divergence_cleaner=self%physics%D_divergence_cleaner,       &
                            B_divergence_cleaner=self%physics%B_divergence_cleaner,       &
                            dxyz_gpu=self%field_gpu%dxyz_gpu,                             &
                            l_map_bc_gpu=self%field_gpu%maps%local_map_bc_crown_gpu,      &
                            fec_1_6_array_gpu=self%field_gpu%fec_1_6_array_gpu,           &
                            q_gpu=q_gpu)
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self)
   !< Set initial conditions of field.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   call self%ic%set_initial_conditions(physics=self%physics, field=self%field)
   call self%coil%set_coils(physics=self%physics, field=self%field)
   call self%copy_cpu_gpu
   endsubroutine set_initial_conditions

   subroutine update_ghost(self, q_gpu, step)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(prism_fnl_object), intent(inout)        :: self            !< The equation.
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

   if (do_local_update) call self%field_gpu%update_ghost_local_gpu(q_gpu=q_gpu)
                        call self%field_gpu%update_ghost_mpi_gpu(q_gpu=q_gpu, step=step)
   if (do_set_bc)       call self%set_boundary_conditions(q_gpu=q_gpu)
   endsubroutine update_ghost

   ! numerical methods
   subroutine compute_dt(self)
   !< Compute maximum time step accordingly to CFL stabilty criterion.
   class(prism_fnl_object), intent(inout) :: self     !< The equation.
   real(R8P)                              :: umax     !< Maximum speed of waves propagation (light speed)
   real(R8P)                              :: dxyz_min !< Minimum space step.

   call compute_dxyz_min_dev(blocks_number=self%blocks_number, dxyz_gpu=self%field_gpu%dxyz_gpu, dxyz_min=dxyz_min)

   if (self%physics%D_divergence_cleaner) then
      umax = max(self%physics%chi*sqrt(1._R8P/(EPS0*MU0)), self%physics%eta*sqrt(1._R8P/(EPS0*MU0)))
   else
      umax = sqrt(1._R8P/(EPS0*MU0))
   endif
   self%time%dt = self%time%CFL*dxyz_min/umax
   call MPI_ALLREDUCE(MPI_IN_PLACE, self%time%dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, self%mpih_gpu%error)
   endsubroutine compute_dt

   subroutine compute_residuals(self, q_gpu, dq_gpu)
   !< Compute residuals of equation.
   class(prism_fnl_object), intent(inout) :: self       !< The equation.
   real(R8P),               intent(inout) :: q_gpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1:)  !< Conservative variables.
   real(R8P),               intent(inout) :: dq_gpu(1:,         &
                                                    1-self%ngc:,&
                                                    1-self%ngc:,&
                                                    1-self%ngc:,&
                                                    1:) !< Residuals.
   real(R8P)                              :: evmax      !< Maximum waves speeds estimation.

   call self%update_ghost(q_gpu=q_gpu)
   ! call self%integrate_eikonal(q_gpu=q_gpu)
   ! call self%compute_q_auxiliary(q_gpu=q_gpu, q_aux_gpu=self%q_aux_gpu)
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv, blocks_number=>self%blocks_number, &
             dx_gpu=>self%field_gpu%dxyz_gpu(:,1),                                                                 &
             dy_gpu=>self%field_gpu%dxyz_gpu(:,2),                                                                 &
             dz_gpu=>self%field_gpu%dxyz_gpu(:,3),                                                                 &
             q_gpu=>self%q_gpu,                                                                                    &
             flx_gpu=>self%flx_gpu, fly_gpu=>self%fly_gpu, flz_gpu=>self%flz_gpu,                                  &
             solids_number=>self%ib%solids_number,                                                                 &
             null_xyz=>self%grid%null_xyz,                                                                         &
             eta=>self%physics%eta, chi=>self%physics%chi, d_divergence_cleaner=>self%physics%d_divergence_cleaner,&
             b_divergence_cleaner=>self%physics%b_divergence_cleaner)
   if (blocks_number > 0) then
      if (d_divergence_cleaner) then
         evmax = chi*sqrt(1._R8P/(EPS0*MU0))
      else
         evmax = sqrt(1._R8P/(EPS0*MU0))
      endif
      if (.not.null_xyz(1)) then
         call compute_fluxes_convective_dev(dir=1,                                                            &
                                            blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, &
                                            evmax=evmax, q_gpu=q_gpu, fluxes_gpu=flx_gpu)
      else
         flx_gpu = 0._R8P
      endif
      if (.not.null_xyz(2)) then
         call compute_fluxes_convective_dev(dir=2,                                                            &
                                            blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, &
                                            evmax=evmax, q_gpu=q_gpu, fluxes_gpu=fly_gpu)
      else
         fly_gpu = 0._R8P
      endif
      if (.not.null_xyz(3)) then
         call compute_fluxes_convective_dev(dir=3,                                                            &
                                            blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, &
                                            evmax=evmax, q_gpu=q_gpu, fluxes_gpu=flz_gpu)
      else
         flz_gpu = 0._R8P
      endif
      if (solids_number>0) then
         call compute_fluxes_difference_dev(null_x=null_xyz(1), null_y=null_xyz(2), null_z=null_xyz(3),                          &
                                            blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, ib_eps=1.e-12_R8P, &
                                            dx_gpu=dx_gpu, dy_gpu=dy_gpu, dz_gpu=dz_gpu,                                         &
                                            flx_gpu=flx_gpu, fly_gpu=fly_gpu, flz_gpu=flz_gpu,                                   &
                                            ! phi_gpu=phi_gpu,                  &
                                            q_gpu=q_gpu, eta=eta, chi=chi, d_divergence_cleaner=d_divergence_cleaner,            &
                                            b_divergence_cleaner=b_divergence_cleaner, dq_gpu=dq_gpu)
      else
         call compute_fluxes_difference_dev(null_x=null_xyz(1), null_y=null_xyz(2), null_z=null_xyz(3),                          &
                                            blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, ib_eps=1.e-12_R8P, &
                                            dx_gpu=dx_gpu, dy_gpu=dy_gpu, dz_gpu=dz_gpu,                                         &
                                            flx_gpu=flx_gpu, fly_gpu=fly_gpu, flz_gpu=flz_gpu,                                   &
                                            q_gpu=q_gpu, eta=eta, chi=chi, d_divergence_cleaner=d_divergence_cleaner,            &
                                            b_divergence_cleaner=b_divergence_cleaner, dq_gpu=dq_gpu)
      endif
      !!Calcolo della divergenza tramite i flussi
      !do b=1, blocks_number
      !   do k=1, nk
      !      do j=1, nj
      !         do i=1, ni
      !            self%field_Div(3,i,j,k,b) = CFL/vmax*(dq(1,i,j,k,b)+dq(2,i,j,k,b)+dq(3,i,j,k,b))
      !            self%field_Div(4,i,j,k,b) = CFL/vmax*(dq(4,i,j,k,b)+dq(5,i,j,k,b)+dq(6,i,j,k,b))
      !         enddo
      !      enddo
      !   enddo
      !enddo
   endif
   endassociate
   endsubroutine compute_residuals

   subroutine integrate(self, do_ghost_syncro)
   !< Perform one step integration.
   class(prism_fnl_object), intent(inout)         :: self             !< The equation.
   logical,                 intent(in),  optional :: do_ghost_syncro  !< Flag to do syncrous ghost update.
   logical                                        :: do_ghost_syncro_ !< Flag to do syncrous ghost update, local var.
   integer(I4P)                                   :: s                !< Counter.

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   ! if (self%coil%total_coils_number >= 1_I4P) then
   !    call compute_coils_current(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, q=self%q, time=time, A=A, d=d, &
   !                               f=freq, phase=phase, coil_flag=coil_flag, td=td, J_vec=J_vec, dx=dx)
   ! endif
   ! call self%rk_gpu%initialize_stages(q_gpu=self%q_gpu)
   select case(self%rk%scheme)
   case(RK_1, RK_2, RK_3)
      ! low storage RK working on q_rk_gpu(:,:,:,:,:,1)/q_gpu as stages, update q_gpu in place
      do s=1, self%rk%nrk
         call self%compute_residuals(q_gpu=self%q_gpu, dq_gpu=self%dq_gpu)
         if (s==1) call self%save_residuals
         if (self%ib%solids_number>0) then
            ! call self%rk_gpu%compute_stage_ls(s=s,dt=self%time%dt,phi_gpu=self%ib_gpu%phi_gpu,dq_gpu=self%dq_gpu,q_gpu=self%q_gpu)
         else
            ! call self%rk_gpu%compute_stage_ls(s=s,dt=self%time%dt,dq_gpu=self%dq_gpu,q_gpu=self%q_gpu)
         endif
      enddo
   case(RK_SSP_22, RK_SSP_33, RK_SSP_54)
      ! RK working on q_rk_gpu as stages
      do s=1, self%rk%nrk
         if (self%ib%solids_number>0) then
            ! call self%rk_gpu%compute_stage(s=s, dt=self%time%dt, phi_gpu=self%ib_gpu%phi_gpu)
         else
            ! call self%rk_gpu%compute_stage(s=s, dt=self%time%dt)
         endif
         call self%compute_residuals(q_gpu=self%rk_gpu%q_rk_gpu(:,:,:,:,:,s), dq_gpu=self%dq_gpu)
         if (s==1) call self%save_residuals
         if (self%ib%solids_number>0) then
            ! call self%rk_gpu%assign_stage(s=s, q_gpu=self%dq_gpu, phi_gpu=self%ib_gpu%phi_gpu)
         else
            ! call self%rk_gpu%assign_stage(s=s, q_gpu=self%dq_gpu)
         endif
      enddo
      if (self%ib%solids_number>0) then
         ! call self%rk_gpu%update_q(dt=self%time%dt, phi_gpu=self%ib_gpu%phi_gpu, q_gpu=self%q_gpu)
      else
         ! call self%rk_gpu%update_q(dt=self%time%dt, q_gpu=self%q_gpu)
      endif
   endselect

   !!calcolo della divergenza tramite differenze finite
   !do b=1, blocks_number
   !   do k=1, nk
   !      do j=1, nj
   !         do i=1, ni
   !               !if (self%coil%coil_flag(i,j,k,b) > 0_I4P) then
   !               !   self%q(1:6,i,j,k,b) = 0._R8P ! azzero i campi dentro le spire
   !               !endif
   !            self%field_Div(1,i,j,k,b) = 0.5_R8P*((self%q(1,i+1,j,k,b) - self%q(1,i-1,j,k,b))/dx + &
   !                                       (self%q(2,i,j+1,k,b) - self%q(2,i,j-1,k,b))/dx + &
   !                                       (self%q(3,i,j,k+1,b) - self%q(3,i,j,k-1,b))/dx)
   !            self%field_Div(2,i,j,k,b) = 0.5_R8P*((self%q(4,i+1,j,k,b) - self%q(4,i-1,j,k,b))/dx + &
   !                                       (self%q(5,i,j+1,k,b) - self%q(5,i,j-1,k,b))/dx + &
   !                                       (self%q(6,i,j,k+1,b) - self%q(6,i,j,k-1,b))/dx)
   !         enddo
   !      enddo
   !   enddo
   !enddo
   endsubroutine integrate

   subroutine simulate(self, filename)
   !< Perform the simulation.
   class(prism_fnl_object), intent(inout) :: self             !< The equation.
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
         ! if (self%ib%solids_number > 0) call self%compute_phi
         ! call self%amr_update
      enddo
      call self%set_initial_conditions
      self%time%time = 0._R8P
      self%time%it = 0
      call self%mpih_gpu%print_message('impose initial conditions finish')
   endif
   ! if (self%ib%solids_number > 0) call self%compute_phi
   ! call self%amr_update
   call self%save_simulation_data
   if (self%mpih_gpu%myrank==0) call self%io%open_file_residuals(nv=self%nv)

   ! integration
   call self%mpih_gpu%barrier(tictoc=.true., timing=timing(1), single=.true.)
   integration: do
      call self%mpih_gpu%barrier(tictoc=.true., timing=timing_step(1), single=.true.)
      self%time%it = self%time%it + 1

      if (self%io%save_memory_status) then
         call save_memory_status_cpu(file_name='memory_cpu-'//self%mpih_gpu%myrankstr//'.dat', tag=str(self%time%it,.true.))
         call save_memory_status_gpu(file_name='memory_gpu-'//self%mpih_gpu%myrankstr//'.dat', tag=str(self%time%it,.true.))
      endif

      if (mod(self%time%it,self%amr%frequency)==0) then
         call self%mpih_gpu%barrier(tictoc=.true.)
         ! call self%amr_update
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
endmodule adam_prism_fnl_object
