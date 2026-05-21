!< ADAM, Maxwell equations system class definition, common data to all backends.
module adam_prism_common_object

! ADAM common library
use :: adam_common_library
! PRISM classes, libraries, parameters
use :: adam_prism_bc_object
use :: adam_prism_coil_object
use :: adam_prism_external_fields_object
use :: adam_prism_fWLayer_object
use :: adam_prism_ic_object
use :: adam_prism_leapfrog_pic_object
use :: adam_prism_numerics_object
use :: adam_prism_physics_object
use :: adam_prism_pic_object
use :: adam_prism_particle_injection_object
use :: adam_prism_rk_pic_object
use :: adam_prism_rk_bc_object
use :: adam_prism_time_object
! third party modules
use :: motion
use :: penf
use :: stringifor

implicit none
private
public :: prism_common_object

type, extends(realm_object) :: prism_common_object
   !< Maxwell equations system class definition, common data to all backends.
   ! physics data replica for easy handling
   integer(I4P), pointer :: nv_c=>null()  !< Number of conservative variables in q vector.
   integer(I4P), pointer :: nv_s=>null()  !< Number of source variables in q vector.
   integer(I4P), pointer :: nv_cl=>null() !< Number of divergence cleaning variables in q vector.
   ! fields data [1:nv,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb].
   real(R8P),    allocatable :: q(         :,:,:,:,:) !< Conservative cell centered variables.
   real(R8P),    allocatable :: dq(        :,:,:,:,:) !< Residuals right hand side.
   real(R8P),    allocatable :: q_pic(           :,:) !< PIC variables.
   real(R8P),    allocatable :: pic_fields(      :,:) !< Fields value at particle locations.
   real(R8P),    allocatable :: curl(      :,:,:,:,:) !< Curl fields.
   real(R8P),    allocatable :: divergence(:,:,:,:,:) !< Divergence fields.
   type(string), allocatable :: q_name(:)             !< Fields names [1:nv].
   type(string), allocatable :: dq_name(:)            !< Residuals names [1:nv].
   type(string), allocatable :: q_pic_name(:)         !< PIC Fields names.
   type(string), allocatable :: curl_name(:)          !< Curl fields names.
   type(string), allocatable :: div_name(:)           !< Divergence fields names.
   ! auxiliary data
   real(R8P), allocatable :: energy_D(:)                !< Energy of field D, time history.
   real(R8P), allocatable :: energy_B(:)                !< Energy of field B, time history.
   real(R8P), allocatable :: coil_power(:)              !< Power of coils, time history.
   real(R8P), allocatable :: Poynting_flux(:)           !< Total Poynting flux from boundary, time history.
   real(R8P)              :: rms_energy_error_D=0.0_R8P !< RMS energy error of D field.
   real(R8P)              :: rms_energy_error_B=0.0_R8P !< RMS energy error of B field.
   real(R8P)              :: max_divergence_D=0.0_R8P   !< Maximum of divergence of D field.
   real(R8P)              :: max_divergence_B=0.0_R8P   !< Maximum of divergence of B field.
   real(R8P)              :: max_divergence_J=0.0_R8P   !< Maximum of divergence of J field.
   ! PRISM C.3 closure (issue #13 D.4b follow-up): the 13 PRISM-app objects
   ! that previously lived as module-scope singletons are now per-realm value
   ! components. The 13 adam_prism_*_global modules became pointer shims that
   ! alias these components via adam_prism_common_bind. Without this, N>1
   ! runs (rmf-2realm) crashed on realm 2's prism_bc_object%initialize
   ! ("allocate already allocated") because realm 1 had already initialized
   ! the shared singleton.
   type(prism_bc_object)                 :: bc                 !< Boundary conditions.
   type(prism_coil_object)               :: coil               !< Coil source term.
   type(prism_external_fields_object)    :: external_fields    !< External fields.
   type(prism_fWLayer_object)            :: fWLayer            !< Far-field weighted layer.
   type(prism_ic_object)                 :: ic                 !< Initial conditions.
   type(prism_leapfrog_pic_object)       :: leapfrog_pic       !< Leapfrog PIC integrator.
   type(prism_numerics_object)           :: numerics           !< Numerics configuration.
   type(prism_particle_injection_object) :: particle_injection !< Particle injection.
   type(prism_physics_object)            :: physics            !< Physics configuration.
   type(prism_pic_object)                :: pic                !< PIC state.
   type(prism_rk_bc_object)              :: rk_bc              !< Runge-Kutta BC handler.
   type(prism_rk_pic_object)             :: rk_pic             !< Runge-Kutta PIC handler.
   type(prism_time_object)               :: time               !< Time integration state.
   contains
      procedure, pass(self) :: allocate_common                              !< Allocate common data.
      procedure, pass(self) :: compute_auxiliary_fields                     !< Compute auxiliary fields.
      procedure, pass(self) :: initialize                                   !< Initialize the equation common data.
      ! IO methods
      procedure, pass(self) :: load_restart_files      !< Load restart files.
      procedure, pass(self) :: save_energy_error       !< Save energy error history.
      procedure, pass(self) :: save_energy_history     !< Save energy history.
      procedure, pass(self) :: save_divergence_history !< Save divergence history.
      procedure, pass(self) :: save_restart_files      !< Save restart files.
      procedure, pass(self) :: save_xh5f               !< Save simulation data in XH5F format.
      ! coils initialization methods
      procedure, pass(self) :: initialize_coils       !< Initialize coils.
      procedure, pass(self) :: set_rectangular_coil_x !< Subroutine to set a rectangular coil source with +-x normal
      procedure, pass(self) :: set_rectangular_coil_y !< Subroutine to set a rectangular coil source with +-y normal
      procedure, pass(self) :: set_rectangular_coil_z !< Subroutine to set a rectangular coil source with +-z normal
      procedure, pass(self) :: set_circular_coil_x    !< Subroutine to set a circular coil source with +-x normal
      procedure, pass(self) :: set_circular_coil_y    !< Subroutine to set a circular coil source with +-y normal
      procedure, pass(self) :: set_circular_coil_z    !< Subroutine to set a circular coil source with +-z normal
      procedure, pass(self) :: set_solenoid_x         !< Subroutine to set a solenoid source with +-x normal
      procedure, pass(self) :: set_solenoid_y         !< Subroutine to set a solenoid source with +-y normal
      procedure, pass(self) :: set_solenoid_z         !< Subroutine to set a solenoid source with +-z normal
endtype prism_common_object

contains
   ! public methods
   subroutine allocate_common(self)
   !< Allocate common data.
   class(prism_common_object), intent(inout) :: self !< The equation.

   associate(nv=>self%nv, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, nb=>self%nb, &
             particle_number=>self%pic%particle_number)
   call allocate_variable(var=self%q,                &
                          ulb=reshape([1,nv,         &
                                       1-ngc,ni+ngc, &
                                       1-ngc,nj+ngc, &
                                       1-ngc,nk+ngc, &
                                       1,nb],[2,5]), &
                          msg=mpih%myrankstr//'prism_common_object%allocate_common(q) ', verbose=.true.)
   self%q = 0._R8P
   call allocate_variable(var=self%dq,               &
                          ulb=reshape([1,nv,         &
                                       1-ngc,ni+ngc, &
                                       1-ngc,nj+ngc, &
                                       1-ngc,nk+ngc, &
                                       1,nb],[2,5]), &
                          msg=mpih%myrankstr//'prism_common_object%allocate_common(dq) ', verbose=.true.)
   self%dq = 0._R8P
   call allocate_variable(var=self%divergence,        &
                          ulb=reshape([1,self%nv,    &
                                       1-ngc,ni+ngc, &
                                       1-ngc,nj+ngc, &
                                       1-ngc,nk+ngc, &
                                       1,nb],[2,5]), &
                          msg=mpih%myrankstr//'prism_common_object%allocate_common(divergence) ', verbose=.true.)
   self%divergence = 0._R8P

   call allocate_variable(var=self%curl,             &
                          ulb=reshape([1,self%nv,    &
                                       1-ngc,ni+ngc, &
                                       1-ngc,nj+ngc, &
                                       1-ngc,nk+ngc, &
                                       1,nb],[2,5]), &
                          msg=mpih%myrankstr//'prism_common_object%allocate_common(curl) ', verbose=.true.)
   self%curl = 0._R8P

   allocate(self%q_name(   self%nv))
   allocate(self%dq_name(  self%nv))
   allocate(self%div_name( self%nv))
   allocate(self%curl_name(self%nv))

   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
      ! allocate(q_pic_name(:))
      call allocate_variable(var=self%q_pic,                          &
                             ulb=reshape([1,8,                        &
                                          1,particle_number],[2,2]),  &
                             msg=mpih%myrankstr//'prism_common_object%allocate_common(q_pic) ', verbose=.true.)
      self%q_pic = 0._R8P
      call allocate_variable(var=self%pic_fields,                     &
                             ulb=reshape([1,6,                        &
                                          1,particle_number],[2,2]),  &
                             msg=mpih%myrankstr//'prism_common_object%allocate_common(pic_fields) ', verbose=.true.)
      self%pic_fields = 0._R8P
   endif
   endassociate
   endsubroutine allocate_common

   subroutine compute_auxiliary_fields(self)
   !< Compute auxiliary fields.
   class(prism_common_object), intent(inout) :: self !< The equation.

   associate(hs=>self%fdv_half_stencil,var_jx=>self%physics%var_jx)
   if (self%io%save_divergence_fields) then
      call self%compute_divergence(hs=hs,ivar=VAR_DX,q=self%q,divergence=self%divergence(1,:,:,:,:))
      call self%compute_divergence(hs=hs,ivar=VAR_BX,q=self%q,divergence=self%divergence(2,:,:,:,:))
      call self%compute_divergence(hs=hs,ivar=var_jx,q=self%q,divergence=self%divergence(3,:,:,:,:))
   endif
   if (self%io%save_curl_fields) then
      call self%compute_curl(hs=hs,ivar=VAR_DX,q=self%q,curl=self%curl(1:3,:,:,:,:))
      call self%compute_curl(hs=hs,ivar=VAR_BX,q=self%q,curl=self%curl(4:6,:,:,:,:))
      call self%compute_curl(hs=hs,ivar=var_jx,q=self%q,curl=self%curl(7:9,:,:,:,:))
   endif
   endassociate
   endsubroutine compute_auxiliary_fields

   subroutine initialize(self, filename, memory_avail, nv, verbose)
   !< Initialize the equation common data.
   class(prism_common_object), intent(inout), target :: self         !< The equation.
   character(*),               intent(in)            :: filename     !< Input file name.
   real(R8P),                  intent(in), value     :: memory_avail !< Memory available for single MPI process.
   integer(I4P),               intent(in), optional  :: nv           !< Number of field variables.
   logical,                    intent(in), optional  :: verbose      !< Trigger verbose output.
   logical                                           :: verbose_     !< Trigger verbose output, local variable.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   call mpih%initialize(verbose=verbose_)
   if (verbose_) call mpih%print_message('prism_common_object%initialize start')
   call self%io%initialize(filename=trim(filename),verbose=verbose_)
   associate(file_parameters=>self%io%file_parameters)
   call self%numerics%initialize(file_parameters=file_parameters)
   call self%physics%initialize(file_parameters=file_parameters,                              &
                                reconstruction_vars=self%numerics%reconstruction_vars,        &
                                div_corr_var=self%numerics%div_corr_var,                      &
                                constrained_transport_D=self%numerics%constrained_transport_D,&
                                constrained_transport_B=self%numerics%constrained_transport_B)
   self%nv_c   => self%physics%nv_c
   self%nv_s   => self%physics%nv_s
   self%nv_cl  => self%physics%nv_cl
   !self%nv_pic => self%physics%nv_pic
   call self%realm_object%initialize(filename=filename, memory_avail=memory_avail, nv=self%physics%nv, verbose=verbose_)
   call self%bc%initialize(file_parameters=file_parameters)
   call grid%set_bc_type(bc_type=self%bc%bc_type)
   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) &
      call self%pic%initialize(file_parameters=file_parameters)
   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) &
      call self%particle_injection%initialize(file_parameters=file_parameters, pic=self%pic)
   call self%time%initialize(file_parameters=file_parameters)
   call self%ic%initialize(file_parameters=file_parameters)
   call self%fWLayer%initialize(file_parameters=file_parameters, physics=self%physics)
   call self%coil%initialize(file_parameters=file_parameters)
   call self%external_fields%initialize(file_parameters=file_parameters)
   if (self%numerics%scheme_time==NUM_SCHEME_TIME_RUNGE_KUTTA) &
      call self%rk_bc%initialize(file_parameters=file_parameters, rk=rk, physics=self%physics)
   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
      if (self%pic%scheme_time==NUM_SCHEME_TIME_PIC_LEAPFROG) &
         call self%leapfrog_pic%initialize(file_parameters=file_parameters, pic=self%pic)
      if (self%pic%scheme_time==NUM_SCHEME_TIME_PIC_RUNGE_KUTTA) &
         call self%rk_pic%initialize(file_parameters=file_parameters, rk=rk, pic=self%pic)
   endif
   call check_ngc_number
   call self%allocate_common
   if (self%adam%tree%iu_ref_levels>0) then
      call adam%refine_uniform(refinement_levels=self%adam%tree%iu_ref_levels, do_mpi_redistribute=.true., &
                               do_blocks_reorder=.false., q=self%q)
   endif
   call io_initialize
   endassociate
   if (verbose_) call mpih%print_message('prism_common_object%initialize finish')
   contains
      subroutine check_ngc_number
      !< Check if the number of ghost cells is consistent with the numerical schemes used, if not an error is echoed and
      !< the simulation is stop.

      if (self%numerics%scheme_space==NUM_SCHEME_SPACE_WENO) then
         if (self%weno%S > grid%ngc) &
            call mpih%error_stop(msg=': ghost cells number (ngc) must be >= of weno stencil number (weno%S):'//&
                                      ' ngc='//trim(str(grid%ngc))//' weno%S='//trim(str(self%weno%S)))
      endif
      if (self%fdv_half_stencil > grid%ngc) &
         call mpih%error_stop(msg=': ghost cells number (ngc) must be >= of FDV half stencil number (fdv_hs):'//&
                                   ' ngc='//trim(str(grid%ngc))//' fdv_hs='//trim(str(self%fdv_half_stencil)))
      endsubroutine check_ngc_number

      subroutine io_initialize
      !< Initialize IO data.
      character(9),  allocatable :: q_name(:) !< Variables names buffer.
      logical                    :: add_rho   !< Flag to add rho for PIC.
      integer(I4P)               :: v         !< Counter.

      add_rho = self%physics%physical_model == PIC_PHYSICAL_MODEL
      associate(add_phi=>self%numerics%constrained_transport_D, &
                add_psi=>self%numerics%constrained_transport_B)
                   q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ']
      if (add_phi) q_name = [q_name, 'phi']
      if (add_psi) q_name = [q_name, 'psi']
                   q_name = [q_name, 'Jx ','Jy ','Jz ']
      if (add_rho) q_name = [q_name, 'rho']
      do v=1, size(q_name,dim=1)
         self%q_name(v) = trim(q_name(v))
      enddo
      if (allocated(q_name)) deallocate(q_name)

                   q_name = ['res_Dx ','res_Dy ','res_Dz ','res_Bx ','res_By ','res_Bz ']
      if (add_phi) q_name = [q_name, 'res_phi']
      if (add_psi) q_name = [q_name, 'res_psi']
                   q_name = [q_name, 'res_Jx ','res_Jy ','res_Jz ']
      if (add_rho) q_name = [q_name, 'res_rho']
      do v=1, size(q_name,dim=1)
         self%dq_name(v) = trim(q_name(v))
      enddo
      if (allocated(q_name)) deallocate(q_name)

                   q_name = ['div_D ','div_B ','div_J ','totvar','div_05','div_06']
      if (add_phi) q_name = [q_name, 'div_07']
      if (add_psi) q_name = [q_name, 'div_08']
                   q_name = [q_name, 'div_09','div_10','div_11']
      if (add_rho) q_name = [q_name, 'div_12']
      do v=1, size(q_name,dim=1)
         self%div_name(v) = trim(q_name(v))
      enddo
      if (allocated(q_name)) deallocate(q_name)

                   q_name = ['curl_Dx','curl_Dy','curl_Dz','curl_Bx','curl_By','curl_Bz']
      if (add_phi) q_name = [q_name, 'curl_07']
      if (add_psi) q_name = [q_name, 'curl_08']
                   q_name = [q_name, 'curl_Jx','curl_Jy','curl_Jz']
      if (add_rho) q_name = [q_name, 'curl_12']
      do v=1, size(q_name,dim=1)
         self%curl_name(v) = trim(q_name(v))
      enddo
      endassociate
      endsubroutine io_initialize
   endsubroutine initialize

   ! IO methods
   subroutine load_restart_files(self, t, time)
   !< Save restart files.
   class(prism_common_object), intent(inout) :: self !< The equation.
   integer(I4P),               intent(out)   :: t    !< Time iteration.
   real(R8P),                  intent(out)   :: time !< Time.

   call adam%load_restart_files(basename=self%io%restart_basename, t=t, time=time, q=self%q)
   call adam%make_comm_local_maps_ghost_bc
   endsubroutine load_restart_files

   subroutine save_divergence_history(self, is_to_open, is_to_close)
   !< Save divergence history.
   class(prism_common_object), intent(inout)        :: self        !< The equation.
   logical,                    intent(in), optional :: is_to_open  !< Flag to open  file before first saving.
   logical,                    intent(in), optional :: is_to_close !< Flag to close file after last saving.
   real(R8P)                                        :: max_div_D   !< Maximum of divergence of D field.
   real(R8P)                                        :: max_div_B   !< Maximum of divergence of B
   real(R8P)                                        :: max_div_J   !< Maximum of divergence of J field.
   !real(R8P)                                        :: r           !< Auxiliary variable to identify fWL presence

   !associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, C=>self%fWLayer%C, hs=>self%fdv_half_stencils(1))
   !r = nint(real(C)/(real(C)+1_I4P))
   if (self%time%is_to_save(it_save=self%io%divergence_history_save)) then
      !max_div_D = maxval(abs(self%divergence(1,1+r*(C+hs):ni-r*(C+hs-1_I4P),1+r*(C+hs):nj-r*(C+hs-1_I4P), &
      !                                       1+r*(C+hs):nk-r*(C+hs-1_I4P),:)))
      !max_div_B = maxval(abs(self%divergence(2,1+r*(C+hs):ni-r*(C+hs-1_I4P),1+r*(C+hs):nj-r*(C+hs-1_I4P), &
      !                                       1+r*(C+hs):nk-r*(C+hs-1_I4P),:)))
      !max_div_J = maxval(abs(self%divergence(3,1+r*(C+hs):ni-r*(C+hs-1_I4P),1+r*(C+hs):nj-r*(C+hs-1_I4P), &
      !                                       1+r*(C+hs):nk-r*(C+hs-1_I4P),:)))
      call self%io%save_divergence_history(it=self%time%it,time=self%time%time,blocks_number=self%blocks_number,                          &
                                           div_D=self%max_divergence_D,div_B=self%max_divergence_B,div_J=self%max_divergence_J, &
                                           is_to_open=is_to_open,is_to_close=is_to_close)
   endif
   !endassociate
   endsubroutine save_divergence_history

   subroutine save_energy_error(self, is_to_open, is_to_close)
   !< Save energy error history.
   class(prism_common_object), intent(inout)        :: self        !< The equation.
   logical,                    intent(in), optional :: is_to_open  !< Flag to open  file before first saving.
   logical,                    intent(in), optional :: is_to_close !< Flag to close file after last saving.

   if (self%time%is_to_save(it_save=self%io%energy_error_save)) then
      call self%io%save_energy_error(it=self%time%it,time=self%time%time,blocks_number=self%blocks_number,                  &
                                     energy_D=self%energy_D,energy_B=self%energy_B,                                         &
                                     rms_energy_error_D=self%rms_energy_error_D,rms_energy_error_B=self%rms_energy_error_B, &
                                     is_to_open=is_to_open,is_to_close=is_to_close)
   endif
   endsubroutine save_energy_error

   subroutine save_energy_history(self, is_to_open, is_to_close)
   !< Save energy history.
   class(prism_common_object), intent(inout)        :: self        !< The equation.
   logical,                    intent(in), optional :: is_to_open  !< Flag to open  file before first saving.
   logical,                    intent(in), optional :: is_to_close !< Flag to close file after last saving.

   if (self%time%is_to_save(it_save=self%io%energy_history_save)) then
      call self%io%save_energy_history(it=self%time%it,time=self%time%time,blocks_number=self%blocks_number, &
                                       energy_D=self%energy_D,energy_B=self%energy_B,                        &
                                       coil_power=self%coil_power,Poynting_flux=self%Poynting_flux,          &
                                       is_to_open=is_to_open,is_to_close=is_to_close)
   endif
   endsubroutine save_energy_history

   subroutine save_restart_files(self)
   !< Save restart files.
   class(prism_common_object), intent(inout) :: self !< The equation.

   call mpih%barrier(tictoc=.true.)
   call mpih%print_message('save restart files t: '//trim(str(self%time%it,.true.))//', time: '//&
                                    trim(str(self%time%time,.true.)))
   call adam%save_restart_files(basename=self%io%restart_basename, t=self%time%it, time=self%time%time, q=self%q)
   call self%save_xh5f(output_basename=self%io%restart_basename)
   call mpih%barrier(tictoc=.true.)
   endsubroutine save_restart_files

   subroutine save_xh5f(self, output_basename, with_ghost)
   !< Save simulation data in HDF5 format.
   class(prism_common_object), intent(inout)        :: self             !< The equation.
   character(*),               intent(in), optional :: output_basename  !< Output basename.
   logical,                    intent(in), optional :: with_ghost       !< Flag to save ghost cells.
   character(:), allocatable                        :: output_basename_ !< Output basename, local var.
   logical                                          :: with_ghost_      !< Flag to save ghost cells, local var.
   type(xh5f_file_object)                           :: xh5f             !< XH5F file handler.
   integer(I4P)                                     :: ngc              !< Ghost cells saved.
   integer(I4P)                                     :: ijk(2,3)         !< Blocks extents.
   integer(I8P)                                     :: nijk(3)          !< Blocks dimensions.
   character(:), allocatable                        :: bn               !< Block name.
   integer(I4P)                                     :: b, c, v          !< Counter.

   call mpih%barrier(tictoc=.true.)
   call mpih%print_message('save HDF5 files t: '//trim(str(self%time%it,.true.))//', time: '//&
                                trim(str(self%time%time,.true.)))

   output_basename_ = trim(self%io%output_basename)//'-'//trim(strz(self%time%it,9))
   with_ghost_      = .false. ; if (present(with_ghost      )) with_ghost_       = with_ghost

   if (present(output_basename)) output_basename_ = trim(output_basename)

   ngc = 0_I4P ; if (with_ghost_) ngc = grid%ngc
   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk)
   ijk(:,1) = [1-ngc,ni+ngc]
   ijk(:,2) = [1-ngc,nj+ngc]
   ijk(:,3) = [1-ngc,nk+ngc]
   nijk = [ijk(2,1)-ijk(1,1)+1, &
           ijk(2,2)-ijk(1,2)+1, &
           ijk(2,3)-ijk(1,3)+1]
   endassociate
   call self%open_file_xh5f(basename=trim(output_basename_), xh5f=xh5f)
   do b=1, field%blocks_number
      bn = 'block_'//trim(strz(b,9))//'-proc'//trim(strz(mpih%myrank,6))
      call self%open_block_xh5f(xh5f=xh5f, b=b, nijk=nijk, t=self%time%it, time=self%time%time)

      call self%io%save_field(xh5f=xh5f, block_name=bn, ijk=ijk, nijk=nijk, q=self%q(:,:,:,:,b), q_name=self%q_name)

      if (self%coil%total_coils_number>0) then
         do c=1, self%coil%total_coils_number
            call self%io%save_field(xh5f=xh5f, block_name=bn, ijk=ijk, nijk=nijk, &
                                    q=self%coil%j_vec(:,:,:,:,b,c), q_name=self%coil%j_vec_name(:,c))
         enddo
      endif

      if (self%fWLayer%C>0) &
         call self%io%save_field(xh5f=xh5f, block_name=bn, ijk=ijk, nijk=nijk, &
                                 q=self%fWLayer%f(:,:,:,:,b), q_name=self%fWLayer%f_name)

      if (self%io%save_residual_fields) &
         call self%io%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=self%dq(:,:,:,:,b), q_name=self%dq_name)

      if (self%io%save_curl_fields) &
         call self%io%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=self%curl(:,:,:,:,b), q_name=self%curl_name)

      if (self%io%save_divergence_fields) &
         call self%io%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=self%divergence(:,:,:,:,b), q_name=self%div_name)

      call self%close_block_xh5f(xh5f=xh5f)
   enddo
   call self%close_file_xh5f(xh5f=xh5f)

   call mpih%barrier(tictoc=.true.)
   endsubroutine save_xh5f

   ! coils initialization methods
   subroutine initialize_coils(self)
   !< Initialize coils.
   class(prism_common_object), intent(inout) :: self !< The equation.
   integer(I4P)                              :: n    !< Counter.

   do n=1, self%coil%total_coils_number
      selectcase(self%coil%coil_type(n))
      case(COIL_TYPE_RECTANGULAR)
         select case(self%coil%normal(n))
         case(NORMAL_P_X)
            call self%set_rectangular_coil_x(n=n, verse = 1._R8P)
         case(NORMAL_P_Y)
            call self%set_rectangular_coil_y(n=n, verse = 1._R8P)
         case(NORMAL_P_Z)
            call self%set_rectangular_coil_z(n=n, verse = 1._R8P)
         case(NORMAL_M_X)
            call self%set_rectangular_coil_x(n=n, verse = -1._R8P)
         case(NORMAL_M_Y)
            call self%set_rectangular_coil_y(n=n, verse = -1._R8P)
         case(NORMAL_M_Z)
            call self%set_rectangular_coil_z(n=n, verse = -1._R8P)
         endselect
      case(COIL_TYPE_CIRCULAR)
         select case(self%coil%normal(n))
         case(NORMAL_P_X)
            call self%set_circular_coil_x(n=n, verse = 1._R8P)
         case(NORMAL_P_Y)
            call self%set_circular_coil_y(n=n, verse = 1._R8P)
         case(NORMAL_P_Z)
            call self%set_circular_coil_z(n=n, verse = 1._R8P)
         case(NORMAL_M_X)
            call self%set_circular_coil_x(n=n, verse = -1._R8P)
         case(NORMAL_M_Y)
            call self%set_circular_coil_y(n=n, verse = -1._R8P)
         case(NORMAL_M_Z)
            call self%set_circular_coil_z(n=n, verse = -1._R8P)
         endselect
      case(COIL_TYPE_SOLENOID)
         select case(self%coil%normal(n))
         case(NORMAL_P_X)
            call self%set_solenoid_x(n=n, verse = 1._R8P)
         case(NORMAL_P_Y)
            call self%set_solenoid_y(n=n, verse = 1._R8P)
         case(NORMAL_P_Z)
            call self%set_solenoid_z(n=n, verse = 1._R8P)
         case(NORMAL_M_X)
            call self%set_solenoid_x(n=n, verse = -1._R8P)
         case(NORMAL_M_Y)
            call self%set_solenoid_y(n=n, verse = -1._R8P)
         case(NORMAL_M_Z)
            call self%set_solenoid_z(n=n, verse = -1._R8P)
         endselect
      endselect
      !call self%compute_divergence(hs=self%fdv_half_stencils(1),ivar=1_I4P,q=self%coil%J_vec(1:3,:,:,:,:,n),&
      !                             divergence=self%divergence(3,:,:,:,:))
      !print '(A)', mpih%myrankstr//'Divergenza J vec della spira: ' &
      !            //trim(str(n))//' pari a: '//trim(str(maxval(abs(self%divergence(3,:,:,:,:)))))
   enddo
   endsubroutine initialize_coils

   subroutine set_rectangular_coil_x(self, n, verse)
   !< Set rectangular coil with normal direction parallel to x.
   class(prism_common_object), intent(inout) :: self                    !< Cpu object.
   integer(I4P),               intent(in)    :: n                       !< Coil number.
   real(R8P),                  intent(in)    :: verse                   !< Coil normal direction, +1=+x, -1=-x.
   real(R8P), allocatable                    :: A(:,:,:,:,:)            !< Total coil vector potential field.
   real(R8P), allocatable                    :: J_vec_buffer(:,:,:,:,:) !< Buffer variable for self%coil%J_vec.
   real(R8P), allocatable                    :: A_gc(:,:,:,:)           !< Vector potential contribution from ghost cell reconstruction.
   real(R8P)                                 :: A_1                     !< Vector potential contribution from side 1.
   real(R8P)                                 :: A_2                     !< Vector potential contribution from side 2.
   real(R8P)                                 :: A_3                     !< Vector potential contribution from side 3.
   real(R8P)                                 :: A_4                     !< Vector potential contribution from side 4.
   real(R8P)                                 :: cell_coord(3)           !< Cell-center coordinate vector.
   real(R8P)                                 :: gc_coord(3)             !< Ghost cell center coordinate vector for ghost cell reconstruction.
   real(R8P)                                 :: y_1                     !< Lower y-boundary of the rectangular coil.
   real(R8P)                                 :: y_2                     !< Upper y-boundary of the rectangular coil.
   real(R8P)                                 :: z_1                     !< Lower z-boundary of the rectangular coil.
   real(R8P)                                 :: z_2                     !< Upper z-boundary of the rectangular coil.
   real(R8P)                                 :: F_n                     !< Normal error-function profile.
   real(R8P)                                 :: W_t                     !< Tangential window function.
   real(R8P)                                 :: W_x                     !< Coil thickness window function along x.
   integer(I4P)                              :: i_dir_n, i_dir_a, i_dir_b !< Normal and tangential directions indices for ghost cell reconstruction.
   integer(I4P), parameter                   :: n_faces = 6_I4P         !< Number of ghost-cell face slabs.
   integer(I4P)                              :: i1_f(n_faces), i2_f(n_faces)
   integer(I4P)                              :: j1_f(n_faces), j2_f(n_faces)
   integer(I4P)                              :: k1_f(n_faces), k2_f(n_faces)
   integer(I4P)                              :: b,i,j,k,f               !< Counters.

   ! Schema locale nel piano della spira (y,z):
   !
   !          z
   !          ^
   !          |
   !   z2  3 <-------- 2
   !       |            |
   !       |            |
   !       |            |
   !   z1  4 --------> 1
   !          y1       y2  ---> y
   !
   ! Scelta locale:
   !   a = y
   !   b = z
   !   a x b = y x z = x
   !
   ! Lati, per normale +x:
   !   1: z = z_1, tangente lungo +y
   !   2: y = y_2, tangente lungo +z
   !   3: z = z_2, tangente lungo -y
   !   4: y = y_1, tangente lungo -z

   ! Associate grid, field, and coil data.
   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,         &
            nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                     &
            lx=>self%coil%lx(n), ly=>self%coil%ly(n),                                        &
            y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), normal=>self%coil%normal(n),  &
            nb=>field%nb, x_cell=>field%x_cell, y_cell=>field%y_cell,              &
            z_cell=>field%z_cell, sigma=>self%coil%sigma(n), J_vec=>self%coil%J_vec,         &
            dxyz=>field%dxyz, hs=>self%fdv_half_stencil)

   ! Set rectangular coil boundaries with normal direction parallel to x
   ! and center located at (x_c, y_c, z_c).
   y_1 = -lx/2.0_R8P + y_c
   y_2 = +lx/2.0_R8P + y_c
   z_1 = -ly/2.0_R8P + z_c
   z_2 = +ly/2.0_R8P + z_c

   allocate(A(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   A(:,:,:,:,:) = 0.0_R8P

   allocate(A_gc(1:3, 1:2*hs+1, 1:2*hs+1, 1:2*hs+1))
   A_gc(:,:,:,:) = 0.0_R8P

   allocate(J_vec_buffer(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   J_vec_buffer(:,:,:,:,:) = 0.0_R8P

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]

               ! Finestra nella direzione normale, uguale per tutti i lati.
               W_x = tangential_window(s=cell_coord(1), smin=x_c-sigma, smax=x_c+sigma, sigma=sigma)

               ! Side 1: z = z_1, tangential direction along +y.
               F_n = erf_function(s=cell_coord(3), mu=z_1, sigma=sigma)
               W_t = tangential_window(s=cell_coord(2), smin=y_1, smax=y_2, sigma=sigma)
               A_1 = F_n*W_t*W_x

               ! Side 2: y = y_2, tangential direction along +z.
               F_n = erf_function(s=cell_coord(2), mu=y_2, sigma=sigma)
               W_t = tangential_window(s=cell_coord(3), smin=z_1, smax=z_2, sigma=sigma)
               A_2 = -F_n*W_t*W_x

               ! Side 3: z = z_2, tangential direction along -y.
               F_n = erf_function(s=cell_coord(3), mu=z_2, sigma=sigma)
               W_t = tangential_window(s=cell_coord(2), smin=y_1, smax=y_2, sigma=sigma)
               A_3 = -F_n*W_t*W_x

               ! Side 4: y = y_1, tangential direction along -z.
               F_n = erf_function(s=cell_coord(2), mu=y_1, sigma=sigma)
               W_t = tangential_window(s=cell_coord(3), smin=z_1, smax=z_2, sigma=sigma)
               A_4 = F_n*W_t*W_x

               ! Sum vector potential contributions.
               A(1,i,j,k,b) = A_1 + A_2 + A_3 + A_4
            enddo
         enddo
      enddo
   enddo

   call self%compute_curl(hs=hs, ivar=1_I4P, q=A, curl=J_vec_buffer)

   ! Parte relativa alle ghost cells, visto che compute_curl lavora solo
   ! sui punti interni del dominio.

   i_dir_n = 1_I4P
   i_dir_a = 2_I4P
   i_dir_b = 3_I4P

   ! Faccia -x
   i1_f(1) = 1_I4P - ngc
   i2_f(1) = 0_I4P
   j1_f(1) = 1_I4P
   j2_f(1) = nj
   k1_f(1) = 1_I4P
   k2_f(1) = nk

   ! Faccia +x
   i1_f(2) = ni + 1_I4P
   i2_f(2) = ni + ngc
   j1_f(2) = 1_I4P
   j2_f(2) = nj
   k1_f(2) = 1_I4P
   k2_f(2) = nk

   ! Faccia -y
   i1_f(3) = 1_I4P
   i2_f(3) = ni
   j1_f(3) = 1_I4P - ngc
   j2_f(3) = 0_I4P
   k1_f(3) = 1_I4P
   k2_f(3) = nk

   ! Faccia +y
   i1_f(4) = 1_I4P
   i2_f(4) = ni
   j1_f(4) = nj + 1_I4P
   j2_f(4) = nj + ngc
   k1_f(4) = 1_I4P
   k2_f(4) = nk

   ! Faccia -z
   i1_f(5) = 1_I4P
   i2_f(5) = ni
   j1_f(5) = 1_I4P
   j2_f(5) = nj
   k1_f(5) = 1_I4P - ngc
   k2_f(5) = 0_I4P

   ! Faccia +z
   i1_f(6) = 1_I4P
   i2_f(6) = ni
   j1_f(6) = 1_I4P
   j2_f(6) = nj
   k1_f(6) = nk + 1_I4P
   k2_f(6) = nk + ngc

   do f = 1_I4P, n_faces
      do b = 1, blocks_number
         do k = k1_f(f), k2_f(f)
            do j = j1_f(f), j2_f(f)
               do i = i1_f(f), i2_f(f)

                  gc_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]

                  A_gc = build_A_ghost_stencil(hs=hs, dxyz=dxyz(:,b), gc_coord=gc_coord,            &
                                               x_c=x_c, y_c=y_c, z_c=z_c, sigma=sigma,              &
                                               i_dir_n=i_dir_n, i_dir_a=i_dir_a, i_dir_b=i_dir_b,   &
                                               a1=y_1, a2=y_2, b1=z_1, b2=z_2)

                  call compute_curl_fd_centered(s=hs, dxyz=dxyz(:,b), q=A_gc,                       &
                                                curl=J_vec_buffer(:,i,j,k,b))

               enddo
            enddo
         enddo
      enddo
   enddo

   J_vec_buffer = J_vec_buffer * verse

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (maxval(abs(J_vec_buffer(:,i,j,k,b))) < 1.0e-12_R8P) then
                  J_vec_buffer(:,i,j,k,b) = 0.0_R8P
               endif
            enddo
         enddo
      enddo
   enddo

   !J_vec(1:3,:,:,:,:,n) = J_vec(1:3,:,:,:,:,n) + J_vec_buffer
   J_vec(1:3,:,:,:,:,n) = J_vec_buffer

   if (n == 1_I4P) then
      call compute_coil_current_density_flux_analytic_x(x_c=x_c, y_c=y_c, z_c=z_c, y_1=y_1, y_2=y_2, z_1=z_1,   &
                                                        z_2=z_2, A=self%coil%A(n), amplitude=self%coil%coil_amplitude(n), &
                                                        sigma=sigma, n=n, adjust_amplitude=.true.)
   else
      self%coil%coil_amplitude(n) = self%coil%coil_amplitude(1)
      call compute_coil_current_density_flux_analytic_x(x_c=x_c, y_c=y_c, z_c=z_c, y_1=y_1, y_2=y_2, z_1=z_1,   &
                                                        z_2=z_2, A=self%coil%A(n), amplitude=self%coil%coil_amplitude(n), &
                                                        sigma=sigma, n=n, adjust_amplitude=.false.)
   endif

   endassociate

   contains

      subroutine compute_coil_current_density_flux_analytic_x(x_c, y_c, z_c, y_1, y_2, z_1, z_2, sigma, &
                                                               A, amplitude, n, adjust_amplitude)
      !< Analytic amplitude correction for a rectangular coil normal to x.
      !< The current is evaluated as the flux of Jy through the lower side z = z_1.
      real(R8P),    intent(in)    :: x_c              !< Coil center x-coordinate.
      real(R8P),    intent(in)    :: y_c              !< Coil center y-coordinate.
      real(R8P),    intent(in)    :: z_c              !< Coil center z-coordinate.
      real(R8P),    intent(in)    :: y_1              !< Lower y-boundary of the rectangular coil.
      real(R8P),    intent(in)    :: y_2              !< Upper y-boundary of the rectangular coil.
      real(R8P),    intent(in)    :: z_1              !< Lower z-boundary of the rectangular coil.
      real(R8P),    intent(in)    :: z_2              !< Upper z-boundary of the rectangular coil.
      real(R8P),    intent(in)    :: sigma            !< Coil Gaussian smoothing width.
      real(R8P),    intent(in)    :: A                !< Coil amplitude, target.
      real(R8P),    intent(inout) :: amplitude        !< Corrected amplitude.
      integer(I4P), intent(in)    :: n                !< Coil number.
      logical,      intent(in)    :: adjust_amplitude !< If true, correct the coil amplitude.
      real(R8P), parameter        :: alpha = 3.5_R8P  !< Half-width of the integration interval in sigma units.
      real(R8P)                   :: x_l              !< Lower x-boundary of the coil thickness window.
      real(R8P)                   :: x_u              !< Upper x-boundary of the coil thickness window.
      real(R8P)                   :: x_s              !< x-coordinate of the flux-section center.
      real(R8P)                   :: y_s              !< y-coordinate of the flux-section plane.
      real(R8P)                   :: z_s              !< z-coordinate of the flux-section center.
      real(R8P)                   :: x_min            !< Lower x-integration boundary.
      real(R8P)                   :: x_max            !< Upper x-integration boundary.
      real(R8P)                   :: z_min            !< Lower z-integration boundary.
      real(R8P)                   :: z_max            !< Upper z-integration boundary.
      real(R8P)                   :: e_yd_ys          !< E_yd evaluated at y_s.
      real(R8P)                   :: e_yt_ys          !< E_yt evaluated at y_s.
      real(R8P)                   :: w_y_ys           !< W_y evaluated at y_s.
      real(R8P)                   :: dw_z             !< Difference W_z(z_max)-W_z(z_min).
      real(R8P)                   :: de_zb            !< Difference E_zb(z_max)-E_zb(z_min).
      real(R8P)                   :: de_zf            !< Difference E_zf(z_max)-E_zf(z_min).
      real(R8P)                   :: i_x              !< Integral of W_x over the x-integration interval.
      real(R8P)                   :: i_y_1            !< Flux contribution associated with A_1.
      real(R8P)                   :: i_y_2            !< Flux contribution associated with A_2.
      real(R8P)                   :: i_y_3            !< Flux contribution associated with A_3.
      real(R8P)                   :: i_y_4            !< Flux contribution associated with A_4.
      real(R8P)                   :: flux_unit        !< Analytic current flux for unit amplitude.
      real(R8P)                   :: flux             !< Analytic current flux before amplitude correction.
      real(R8P)                   :: correction       !< Amplitude correction factor.
      real(R8P)                   :: target_current   !< Target input current.

      target_current = A

      ! Set x-window boundaries:
      ! W_x = 1/2 * [E_x_l - E_x_u].
      x_l = x_c - sigma
      x_u = x_c + sigma

      ! Set flux section for the lower side, where the current is directed along +y.
      x_s = x_c
      y_s = y_c
      z_s = z_1

      x_min = x_s - alpha * sigma
      x_max = x_s + alpha * sigma
      z_min = z_s - alpha * sigma
      z_max = z_s + alpha * sigma

      ! Integral of W_x over x:
      !
      ! I_x = 1/2 * [ P_xl(xmax) - P_xl(xmin)
      !              -P_xu(xmax) + P_xu(xmin) ]
      i_x = 0.5_R8P * (                                              &
            erf_primitive_function(s=x_max, mu=x_l, sigma=sigma)     &
          - erf_primitive_function(s=x_min, mu=x_l, sigma=sigma)     &
          - erf_primitive_function(s=x_max, mu=x_u, sigma=sigma)     &
          + erf_primitive_function(s=x_min, mu=x_u, sigma=sigma))

      ! Evaluate y-dependent functions at the flux-section plane.
      e_yd_ys = erf_function(s=y_s, mu=y_1, sigma=sigma)
      e_yt_ys = erf_function(s=y_s, mu=y_2, sigma=sigma)

      w_y_ys = tangential_window(s     = y_s,     &
                                 smin  = y_1,     &
                                 smax  = y_2,     &
                                 sigma = sigma)

      ! Compute Delta W_z = W_z(z_max) - W_z(z_min).
      dw_z = tangential_window(s     = z_max,     &
                               smin  = z_1,       &
                               smax  = z_2,       &
                               sigma = sigma)     &
           - tangential_window(s     = z_min,     &
                               smin  = z_1,       &
                               smax  = z_2,       &
                               sigma = sigma)

      ! Compute Delta E_zb = E_zb(z_max) - E_zb(z_min).
      de_zb = erf_function(s=z_max, mu=z_1, sigma=sigma) &
            - erf_function(s=z_min, mu=z_1, sigma=sigma)

      ! Compute Delta E_zf = E_zf(z_max) - E_zf(z_min).
      de_zf = erf_function(s=z_max, mu=z_2, sigma=sigma) &
            - erf_function(s=z_min, mu=z_2, sigma=sigma)

      ! Contributions to Iy from the four vector-potential components:
      !
      ! A_1 =  E_z1(z) W_y(y) W_x(x)
      ! A_2 = -E_y2(y) W_z(z) W_x(x)
      ! A_3 = -E_z2(z) W_y(y) W_x(x)
      ! A_4 =  E_y1(y) W_z(z) W_x(x)
      !
      ! Since A = A_x e_x:
      !
      ! J_y = dA_x/dz
      !
      ! The flux is evaluated through a section normal to y.
      i_y_1 =  w_y_ys  * i_x * de_zb
      i_y_2 = -e_yt_ys * i_x * dw_z
      i_y_3 = -w_y_ys  * i_x * de_zf
      i_y_4 =  e_yd_ys * i_x * dw_z

      flux_unit = abs(i_y_1 + i_y_2 + i_y_3 + i_y_4)
      flux      = target_current * flux_unit

      if (adjust_amplitude) then

         correction = 1.0_R8P / flux_unit

         print '(A)', mpih%myrankstr//'Current before correction: '//trim(str(flux))
         print '(A)', mpih%myrankstr//trim(str(A))//' Amplitude A('//trim(str(n))//') before correction'
         print '(A)', mpih%myrankstr//'Amplitude scaling factor: '//trim(str(correction))

         amplitude = target_current * correction

         print '(A)', mpih%myrankstr//trim(str(amplitude))//' Amplitude A('//trim(str(n))//') after correction'
         print '(A)', mpih%myrankstr//'Final coil current '//trim(str(n))//': '//trim(str(amplitude*flux_unit))

      else

         print '(A)', mpih%myrankstr//'Amplitude A('//trim(str(n))//') not corrected: '//trim(str(amplitude))
         print '(A)', mpih%myrankstr//'Final coil current '//trim(str(n))//': '//trim(str(amplitude*flux_unit))

      endif

      endsubroutine compute_coil_current_density_flux_analytic_x
   endsubroutine set_rectangular_coil_x

   subroutine set_rectangular_coil_y(self, n, verse)
   !< Set rectangular coil with normal direction parallel to y.
   class(prism_common_object), intent(inout) :: self                    !< Cpu object.
   integer(I4P),               intent(in)    :: n                       !< Coil number.
   real(R8P),                  intent(in)    :: verse                   !< Coil normal direction, +1=+y, -1=-y.
   real(R8P), allocatable                    :: A(:,:,:,:,:)            !< Total coil vector potential field.
   real(R8P), allocatable                    :: J_vec_buffer(:,:,:,:,:) !< Buffer variable for self%coil%J_vec.
   real(R8P), allocatable                    :: A_gc(:,:,:,:)           !< Vector potential contribution from ghost cell reconstruction.
   real(R8P)                                 :: A_1                     !< Vector potential contribution from side 1.
   real(R8P)                                 :: A_2                     !< Vector potential contribution from side 2.
   real(R8P)                                 :: A_3                     !< Vector potential contribution from side 3.
   real(R8P)                                 :: A_4                     !< Vector potential contribution from side 4.
   real(R8P)                                 :: cell_coord(3)           !< Cell-center coordinate vector.
   real(R8P)                                 :: gc_coord(3)             !< Ghost cell center coordinate vector for ghost cell reconstruction.
   real(R8P)                                 :: x_1                     !< Lower x-boundary of the rectangular coil.
   real(R8P)                                 :: x_2                     !< Upper x-boundary of the rectangular coil.
   real(R8P)                                 :: z_1                     !< Lower z-boundary of the rectangular coil.
   real(R8P)                                 :: z_2                     !< Upper z-boundary of the rectangular coil.
   real(R8P)                                 :: F_n                     !< Normal error-function profile.
   real(R8P)                                 :: W_t                     !< Tangential window function.
   real(R8P)                                 :: W_y                     !< Coil thickness window function along y.
   integer(I4P)                              :: i_dir_n, i_dir_a, i_dir_b !< Normal and tangential directions indices for ghost cell reconstruction.
   integer(I4P), parameter                   :: n_faces = 6_I4P         !< Number of ghost-cell face slabs.
   integer(I4P)                              :: i1_f(n_faces), i2_f(n_faces)
   integer(I4P)                              :: j1_f(n_faces), j2_f(n_faces)
   integer(I4P)                              :: k1_f(n_faces), k2_f(n_faces)
   integer(I4P)                              :: b,i,j,k,f               !< Counters.

   ! Schema locale nel piano della spira (z,x):
   !
   !          x
   !          ^
   !          |
   !   x2  3 <-------- 2
   !       |            |
   !       |            |
   !       |            |
   !   x1  1 --------> 4
   !          z1       z2  ---> z
   !
   ! Scelta locale:
   !   a = z
   !   b = x
   !   a x b = z x x = y
   !
   ! Lati, per normale +y:
   !   1: x = x_1, tangente lungo +z
   !   2: z = z_2, tangente lungo +x
   !   3: x = x_2, tangente lungo -z
   !   4: z = z_1, tangente lungo -x

   ! Associate grid, field, and coil data.
   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,         &
            nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                     &
            lx=>self%coil%lx(n), ly=>self%coil%ly(n),                                        &
            y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), normal=>self%coil%normal(n),  &
            nb=>field%nb, x_cell=>field%x_cell, y_cell=>field%y_cell,              &
            z_cell=>field%z_cell, sigma=>self%coil%sigma(n), J_vec=>self%coil%J_vec,         &
            dxyz=>field%dxyz, hs=>self%fdv_half_stencil)

   ! Set rectangular coil boundaries with normal direction parallel to y
   ! and center located at (x_c, y_c, z_c).
   !
   ! Since the local in-plane coordinates are:
   !   a = z
   !   b = x
   !
   ! the first side is x = x_1 and the second side is z = z_2.
   x_1 = -lx/2.0_R8P + x_c
   x_2 = +lx/2.0_R8P + x_c
   z_1 = -ly/2.0_R8P + z_c
   z_2 = +ly/2.0_R8P + z_c

   allocate(A(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   A(:,:,:,:,:) = 0.0_R8P

   allocate(A_gc(1:3, 1:2*hs+1, 1:2*hs+1, 1:2*hs+1))
   A_gc(:,:,:,:) = 0.0_R8P

   allocate(J_vec_buffer(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   J_vec_buffer(:,:,:,:,:) = 0.0_R8P

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]

               ! Finestra nella direzione normale, uguale per tutti i lati.
               W_y = tangential_window(s=cell_coord(2), smin=y_c-sigma, smax=y_c+sigma, sigma=sigma)

               ! Side 1: x = x_1, tangential direction along +z.
               F_n = erf_function(s=cell_coord(1), mu=x_1, sigma=sigma)
               W_t = tangential_window(s=cell_coord(3), smin=z_1, smax=z_2, sigma=sigma)
               A_1 = F_n*W_t*W_y

               ! Side 2: z = z_2, tangential direction along +x.
               F_n = erf_function(s=cell_coord(3), mu=z_2, sigma=sigma)
               W_t = tangential_window(s=cell_coord(1), smin=x_1, smax=x_2, sigma=sigma)
               A_2 = -F_n*W_t*W_y

               ! Side 3: x = x_2, tangential direction along -z.
               F_n = erf_function(s=cell_coord(1), mu=x_2, sigma=sigma)
               W_t = tangential_window(s=cell_coord(3), smin=z_1, smax=z_2, sigma=sigma)
               A_3 = -F_n*W_t*W_y

               ! Side 4: z = z_1, tangential direction along -x.
               F_n = erf_function(s=cell_coord(3), mu=z_1, sigma=sigma)
               W_t = tangential_window(s=cell_coord(1), smin=x_1, smax=x_2, sigma=sigma)
               A_4 = F_n*W_t*W_y

               ! Sum vector potential contributions.
               A(2,i,j,k,b) = A_1 + A_2 + A_3 + A_4
            enddo
         enddo
      enddo
   enddo

   call self%compute_curl(hs=hs, ivar=1_I4P, q=A, curl=J_vec_buffer)

   ! Parte relativa alle ghost cells, visto che compute_curl lavora solo
   ! sui punti interni del dominio.

   i_dir_n = 2_I4P
   i_dir_a = 3_I4P
   i_dir_b = 1_I4P

   ! Faccia -x
   i1_f(1) = 1_I4P-ngc
   i2_f(1) = 0_I4P
   j1_f(1) = 1_I4P
   j2_f(1) = nj
   k1_f(1) = 1_I4P
   k2_f(1) = nk

   ! Faccia +x
   i1_f(2) = ni+1_I4P
   i2_f(2) = ni+ngc
   j1_f(2) = 1_I4P
   j2_f(2) = nj
   k1_f(2) = 1_I4P
   k2_f(2) = nk

   ! Faccia -y
   i1_f(3) = 1_I4P
   i2_f(3) = ni
   j1_f(3) = 1_I4P-ngc
   j2_f(3) = 0_I4P
   k1_f(3) = 1_I4P
   k2_f(3) = nk

   ! Faccia +y
   i1_f(4) = 1_I4P
   i2_f(4) = ni
   j1_f(4) = nj+1_I4P
   j2_f(4) = nj+ngc
   k1_f(4) = 1_I4P
   k2_f(4) = nk

   ! Faccia -z
   i1_f(5) = 1_I4P
   i2_f(5) = ni
   j1_f(5) = 1_I4P
   j2_f(5) = nj
   k1_f(5) = 1_I4P-ngc
   k2_f(5) = 0_I4P

   ! Faccia +z
   i1_f(6) = 1_I4P
   i2_f(6) = ni
   j1_f(6) = 1_I4P
   j2_f(6) = nj
   k1_f(6) = nk+1_I4P
   k2_f(6) = nk+ngc

   do f = 1_I4P, n_faces
      do b = 1, blocks_number
         do k = k1_f(f), k2_f(f)
            do j = j1_f(f), j2_f(f)
               do i = i1_f(f), i2_f(f)

                  gc_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]

                  A_gc = build_A_ghost_stencil(hs=hs, dxyz=dxyz(:,b), gc_coord=gc_coord,            &
                                               x_c=x_c, y_c=y_c, z_c=z_c, sigma=sigma, &
                                               i_dir_n=i_dir_n, i_dir_a=i_dir_a, i_dir_b=i_dir_b,   &
                                               a1=z_1, a2=z_2, b1=x_1, b2=x_2)

                  call compute_curl_fd_centered(s=hs, dxyz=dxyz(:,b), q=A_gc,                       &
                                                curl=J_vec_buffer(:,i,j,k,b))

               enddo
            enddo
         enddo
      enddo
   enddo

   J_vec_buffer = J_vec_buffer * verse

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (maxval(abs(J_vec_buffer(:,i,j,k,b))) < 1.0e-12_R8P) then
                  J_vec_buffer(:,i,j,k,b) = 0.0_R8P
               endif
            enddo
         enddo
      enddo
   enddo

   !J_vec(1:3,:,:,:,:,n) = J_vec(1:3,:,:,:,:,n) + J_vec_buffer
   J_vec(1:3,:,:,:,:,n) = J_vec_buffer

   if (n == 1_I4P) then
      call compute_coil_current_density_flux_analytic_y(x_c=x_c, y_c=y_c, z_c=z_c, x_1=x_1, x_2=x_2, z_1=z_1,   &
                                                        z_2=z_2, A=self%coil%A(n), amplitude=self%coil%coil_amplitude(n), &
                                                        sigma=sigma, n=n, adjust_amplitude=.true.)
   else
      self%coil%coil_amplitude(n) = self%coil%coil_amplitude(1)
      call compute_coil_current_density_flux_analytic_y(x_c=x_c, y_c=y_c, z_c=z_c, x_1=x_1, x_2=x_2, z_1=z_1,   &
                                                        z_2=z_2, A=self%coil%A(n), amplitude=self%coil%coil_amplitude(n), &
                                                        sigma=sigma, n=n, adjust_amplitude=.false.)
   endif

   endassociate

   contains

      subroutine compute_coil_current_density_flux_analytic_y(x_c, y_c, z_c, x_1, x_2, z_1, z_2, sigma, &
                                                               A, amplitude, n, adjust_amplitude)
      !< Analytic amplitude correction for a rectangular coil normal to y.
      !< The current is evaluated as the flux of Jz through a section cutting the side x = x_1.
      real(R8P),    intent(in)    :: x_c              !< Coil center x-coordinate.
      real(R8P),    intent(in)    :: y_c              !< Coil center y-coordinate.
      real(R8P),    intent(in)    :: z_c              !< Coil center z-coordinate.
      real(R8P),    intent(in)    :: x_1              !< Lower x-boundary of the rectangular coil.
      real(R8P),    intent(in)    :: x_2              !< Upper x-boundary of the rectangular coil.
      real(R8P),    intent(in)    :: z_1              !< Lower z-boundary of the rectangular coil.
      real(R8P),    intent(in)    :: z_2              !< Upper z-boundary of the rectangular coil.
      real(R8P),    intent(in)    :: sigma            !< Coil Gaussian smoothing width.
      real(R8P),    intent(in)    :: A                !< Coil amplitude, target.
      real(R8P),    intent(inout) :: amplitude        !< Corrected amplitude.
      integer(I4P), intent(in)    :: n                !< Coil number.
      logical,      intent(in)    :: adjust_amplitude !< If true, correct the coil amplitude.
      real(R8P), parameter        :: alpha = 3.5_R8P  !< Half-width of the integration interval in sigma units.
      real(R8P)                   :: y_l              !< Lower y-boundary of the coil thickness window.
      real(R8P)                   :: y_u              !< Upper y-boundary of the coil thickness window.
      real(R8P)                   :: x_s              !< x-coordinate of the flux-section center.
      real(R8P)                   :: y_s              !< y-coordinate of the flux-section center.
      real(R8P)                   :: z_s              !< z-coordinate of the flux-section plane.
      real(R8P)                   :: x_min            !< Lower x-integration boundary.
      real(R8P)                   :: x_max            !< Upper x-integration boundary.
      real(R8P)                   :: y_min            !< Lower y-integration boundary.
      real(R8P)                   :: y_max            !< Upper y-integration boundary.
      real(R8P)                   :: e_z1_zs          !< E_z1 evaluated at z_s.
      real(R8P)                   :: e_z2_zs          !< E_z2 evaluated at z_s.
      real(R8P)                   :: w_z_zs           !< W_z evaluated at z_s.
      real(R8P)                   :: dw_x             !< Difference W_x(x_max)-W_x(x_min).
      real(R8P)                   :: de_x1            !< Difference E_x1(x_max)-E_x1(x_min).
      real(R8P)                   :: de_x2            !< Difference E_x2(x_max)-E_x2(x_min).
      real(R8P)                   :: i_y              !< Integral of W_y over the y-integration interval.
      real(R8P)                   :: i_z_1            !< Flux contribution associated with A_1.
      real(R8P)                   :: i_z_2            !< Flux contribution associated with A_2.
      real(R8P)                   :: i_z_3            !< Flux contribution associated with A_3.
      real(R8P)                   :: i_z_4            !< Flux contribution associated with A_4.
      real(R8P)                   :: flux_unit        !< Analytic current flux for unit amplitude.
      real(R8P)                   :: flux             !< Analytic current flux before amplitude correction.
      real(R8P)                   :: correction       !< Amplitude correction factor.
      real(R8P)                   :: target_current   !< Target input current.

      target_current = A

      ! Set y-window boundaries:
      !
      ! W_y = 1/2 * [E_y_l - E_y_u].
      y_l = y_c - sigma
      y_u = y_c + sigma

      ! Set flux section for side 1, where the current is directed along +z.
      !
      ! Side 1:
      !   x = x_1
      !   J_z > 0 for verse = +1
      !
      ! The section is normal to z and centered at (x_1, y_c, z_c).
      x_s = x_1
      y_s = y_c
      z_s = z_c

      x_min = x_s - alpha * sigma
      x_max = x_s + alpha * sigma
      y_min = y_s - alpha * sigma
      y_max = y_s + alpha * sigma

      ! Integral of W_y over y:
      !
      ! I_y = 1/2 * [ P_yl(ymax) - P_yl(ymin)
      !              -P_yu(ymax) + P_yu(ymin) ]
      i_y = 0.5_R8P * (                                              &
            erf_primitive_function(s=y_max, mu=y_l, sigma=sigma)     &
          - erf_primitive_function(s=y_min, mu=y_l, sigma=sigma)     &
          - erf_primitive_function(s=y_max, mu=y_u, sigma=sigma)     &
          + erf_primitive_function(s=y_min, mu=y_u, sigma=sigma))

      ! Evaluate z-dependent functions at the flux-section plane.
      e_z1_zs = erf_function(s=z_s, mu=z_1, sigma=sigma)
      e_z2_zs = erf_function(s=z_s, mu=z_2, sigma=sigma)

      w_z_zs = tangential_window(s     = z_s,     &
                                 smin  = z_1,     &
                                 smax  = z_2,     &
                                 sigma = sigma)

      ! Compute Delta W_x = W_x(x_max) - W_x(x_min).
      dw_x = tangential_window(s     = x_max,  &
                               smin  = x_1,    &
                               smax  = x_2,    &
                               sigma = sigma)  &
           - tangential_window(s     = x_min,  &
                               smin  = x_1,    &
                               smax  = x_2,    &
                               sigma = sigma)

      ! Compute Delta E_x1 = E_x1(x_max) - E_x1(x_min).
      de_x1 = erf_function(s=x_max, mu=x_1, sigma=sigma) &
            - erf_function(s=x_min, mu=x_1, sigma=sigma)

      ! Compute Delta E_x2 = E_x2(x_max) - E_x2(x_min).
      de_x2 = erf_function(s=x_max, mu=x_2, sigma=sigma) &
            - erf_function(s=x_min, mu=x_2, sigma=sigma)

      ! Contributions to Iz from the four vector-potential components:
      !
      ! A_1 =  E_x1(x) W_z(z) W_y(y)
      ! A_2 = -E_z2(z) W_x(x) W_y(y)
      ! A_3 = -E_x2(x) W_z(z) W_y(y)
      ! A_4 =  E_z1(z) W_x(x) W_y(y)
      !
      ! Since A = A_y e_y:
      !
      ! J_z = dA_y/dx
      !
      ! The flux is evaluated through a section normal to z.
      i_z_1 =  w_z_zs  * i_y * de_x1
      i_z_2 = -e_z2_zs * i_y * dw_x
      i_z_3 = -w_z_zs  * i_y * de_x2
      i_z_4 =  e_z1_zs * i_y * dw_x

      flux_unit = abs(i_z_1 + i_z_2 + i_z_3 + i_z_4)
      flux      = target_current * flux_unit

      if (adjust_amplitude) then

         correction = 1.0_R8P / flux_unit

         print '(A)', mpih%myrankstr//'Current before correction: '//trim(str(flux))
         print '(A)', mpih%myrankstr//trim(str(A))//' Amplitude A('//trim(str(n))//') before correction'
         print '(A)', mpih%myrankstr//'Amplitude scaling factor: '//trim(str(correction))

         amplitude = target_current * correction

         print '(A)', mpih%myrankstr//trim(str(amplitude))//' Amplitude A('//trim(str(n))//') after correction'
         print '(A)', mpih%myrankstr//'Final coil current '//trim(str(n))//': '//trim(str(amplitude*flux_unit))

      else

         print '(A)', mpih%myrankstr//'Amplitude A('//trim(str(n))//') not corrected: '//trim(str(amplitude))
         print '(A)', mpih%myrankstr//'Final coil current '//trim(str(n))//': '//trim(str(amplitude*flux_unit))

      endif

      endsubroutine compute_coil_current_density_flux_analytic_y
   endsubroutine set_rectangular_coil_y

   subroutine set_rectangular_coil_z(self, n, verse)
   !< Set rectangular coil with normal direction parallel to z.
   class(prism_common_object), intent(inout) :: self                    !< Cpu object.
   integer(I4P),               intent(in)    :: n                       !< Coil number.
   real(R8P),                  intent(in)    :: verse                   !< Coil normal direction, +1=+z, -1=-z.
   real(R8P), allocatable                    :: A(:,:,:,:,:)            !< Total coil vector potential field.
   real(R8P), allocatable                    :: J_vec_buffer(:,:,:,:,:) !< Buffer variable for self%coil%J_vec.
   real(R8P), allocatable                    :: A_gc(:,:,:,:)           !< Vector potential contribution from ghost cell reconstruction.
   real(R8P)                                 :: A_1                     !< Vector potential contribution from side 1.
   real(R8P)                                 :: A_2                     !< Vector potential contribution from side 2.
   real(R8P)                                 :: A_3                     !< Vector potential contribution from side 3.
   real(R8P)                                 :: A_4                     !< Vector potential contribution from side 4.
   real(R8P)                                 :: cell_coord(3)           !< Cell-center coordinate vector.
   real(R8P)                                 :: gc_coord(3)             !< Ghost cell center coordinate vector for ghost cell reconstruction.
   real(R8P)                                 :: x_1                     !< Lower x-boundary of the rectangular coil.
   real(R8P)                                 :: x_2                     !< Upper x-boundary of the rectangular coil.
   real(R8P)                                 :: y_1                     !< Lower y-boundary of the rectangular coil.
   real(R8P)                                 :: y_2                     !< Upper y-boundary of the rectangular coil.
   real(R8P)                                 :: F_n                     !< Normal error-function profile.
   real(R8P)                                 :: W_t                     !< Tangential window function.
   real(R8P)                                 :: W_z                     !< Coil thickness window function along z.
   integer(I4P)                              :: i_dir_n, i_dir_a, i_dir_b !< Normal and tangential directions indices for ghost cell reconstruction.
   integer(I4P), parameter                   :: n_faces = 6_I4P         !< Number of ghost-cell face slabs.
   integer(I4P)                              :: i1_f(n_faces), i2_f(n_faces)
   integer(I4P)                              :: j1_f(n_faces), j2_f(n_faces)
   integer(I4P)                              :: k1_f(n_faces), k2_f(n_faces)
   integer(I4P)                              :: b,i,j,k,f               !< Counters.

   ! Schema locale nel piano della spira (x,y):
   !
   !          y
   !          ^
   !          |
   !   y2  3 <-------- 2
   !       |            |
   !       |            |
   !       |            |
   !   y1  4 --------> 1
   !          x1       x2  ---> x
   !
   ! Scelta locale:
   !   a = x
   !   b = y
   !   a x b = x x y = z
   !
   ! Lati, per normale +z:
   !   1: y = y_1, tangente lungo +x
   !   2: x = x_2, tangente lungo +y
   !   3: y = y_2, tangente lungo -x
   !   4: x = x_1, tangente lungo -y

   ! Associate grid, field, and coil data.
   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,         &
            nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                     &
            lx=>self%coil%lx(n), ly=>self%coil%ly(n),                                        &
            y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), normal=>self%coil%normal(n),  &
            nb=>field%nb, x_cell=>field%x_cell, y_cell=>field%y_cell,              &
            z_cell=>field%z_cell, sigma=>self%coil%sigma(n), J_vec=>self%coil%J_vec,         &
            dxyz=>field%dxyz, hs=>self%fdv_half_stencil)

   ! Set rectangular coil boundaries with normal direction parallel to z
   ! and center located at (x_c, y_c, z_c).
   x_1 = -lx/2.0_R8P + x_c
   x_2 = +lx/2.0_R8P + x_c
   y_1 = -ly/2.0_R8P + y_c
   y_2 = +ly/2.0_R8P + y_c

   allocate(A(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   A(:,:,:,:,:) = 0.0_R8P

   allocate(A_gc(1:3, 1:2*hs+1, 1:2*hs+1, 1:2*hs+1))
   A_gc(:,:,:,:) = 0.0_R8P

   allocate(J_vec_buffer(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   J_vec_buffer(:,:,:,:,:) = 0.0_R8P

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]

               ! Finestra nella direzione normale, uguale per tutti i lati.
               W_z = tangential_window(s=cell_coord(3), smin=z_c-sigma, smax=z_c+sigma, sigma=sigma)

               ! Side 1: y = y_1, tangential direction along +x.
               F_n = erf_function(s=cell_coord(2), mu=y_1, sigma=sigma)
               W_t = tangential_window(s=cell_coord(1), smin=x_1, smax=x_2, sigma=sigma)
               A_1 = F_n*W_t*W_z

               ! Side 2: x = x_2, tangential direction along +y.
               F_n = erf_function(s=cell_coord(1), mu=x_2, sigma=sigma)
               W_t = tangential_window(s=cell_coord(2), smin=y_1, smax=y_2, sigma=sigma)
               A_2 = -F_n*W_t*W_z

               ! Side 3: y = y_2, tangential direction along -x.
               F_n = erf_function(s=cell_coord(2), mu=y_2, sigma=sigma)
               W_t = tangential_window(s=cell_coord(1), smin=x_1, smax=x_2, sigma=sigma)
               A_3 = -F_n*W_t*W_z

               ! Side 4: x = x_1, tangential direction along -y.
               F_n = erf_function(s=cell_coord(1), mu=x_1, sigma=sigma)
               W_t = tangential_window(s=cell_coord(2), smin=y_1, smax=y_2, sigma=sigma)
               A_4 = F_n*W_t*W_z

               ! Sum vector potential contributions.
               A(3,i,j,k,b) = A_1 + A_2 + A_3 + A_4
            enddo
         enddo
      enddo
   enddo

   call self%compute_curl(hs=hs, ivar=1_I4P, q=A, curl=J_vec_buffer)

   ! Parte relativa alle ghost cells, visto che compute_curl lavora solo
   ! sui punti interni del dominio.

   i_dir_n = 3_I4P
   i_dir_a = 1_I4P
   i_dir_b = 2_I4P

   ! Faccia -x
   i1_f(1) = 1_I4P - ngc
   i2_f(1) = 0_I4P
   j1_f(1) = 1_I4P
   j2_f(1) = nj
   k1_f(1) = 1_I4P
   k2_f(1) = nk

   ! Faccia +x
   i1_f(2) = ni + 1_I4P
   i2_f(2) = ni + ngc
   j1_f(2) = 1_I4P
   j2_f(2) = nj
   k1_f(2) = 1_I4P
   k2_f(2) = nk

   ! Faccia -y
   i1_f(3) = 1_I4P
   i2_f(3) = ni
   j1_f(3) = 1_I4P - ngc
   j2_f(3) = 0_I4P
   k1_f(3) = 1_I4P
   k2_f(3) = nk

   ! Faccia +y
   i1_f(4) = 1_I4P
   i2_f(4) = ni
   j1_f(4) = nj + 1_I4P
   j2_f(4) = nj + ngc
   k1_f(4) = 1_I4P
   k2_f(4) = nk

   ! Faccia -z
   i1_f(5) = 1_I4P
   i2_f(5) = ni
   j1_f(5) = 1_I4P
   j2_f(5) = nj
   k1_f(5) = 1_I4P - ngc
   k2_f(5) = 0_I4P

   ! Faccia +z
   i1_f(6) = 1_I4P
   i2_f(6) = ni
   j1_f(6) = 1_I4P
   j2_f(6) = nj
   k1_f(6) = nk + 1_I4P
   k2_f(6) = nk + ngc

   do f = 1_I4P, n_faces
      do b = 1, blocks_number
         do k = k1_f(f), k2_f(f)
            do j = j1_f(f), j2_f(f)
               do i = i1_f(f), i2_f(f)

                  gc_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]

                  A_gc = build_A_ghost_stencil(hs=hs, dxyz=dxyz(:,b), gc_coord=gc_coord,            &
                                               x_c=x_c, y_c=y_c, z_c=z_c, sigma=sigma,              &
                                               i_dir_n=i_dir_n, i_dir_a=i_dir_a, i_dir_b=i_dir_b,   &
                                               a1=x_1, a2=x_2, b1=y_1, b2=y_2)

                  call compute_curl_fd_centered(s=hs, dxyz=dxyz(:,b), q=A_gc,                       &
                                                curl=J_vec_buffer(:,i,j,k,b))

               enddo
            enddo
         enddo
      enddo
   enddo

   J_vec_buffer = J_vec_buffer * verse

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (maxval(abs(J_vec_buffer(:,i,j,k,b))) < 1.0e-12_R8P) then
                  J_vec_buffer(:,i,j,k,b) = 0.0_R8P
               endif
            enddo
         enddo
      enddo
   enddo

   !J_vec(1:3,:,:,:,:,n) = J_vec(1:3,:,:,:,:,n) + J_vec_buffer
   J_vec(1:3,:,:,:,:,n) = J_vec_buffer

   if (n == 1_I4P) then
      call compute_coil_current_density_flux_analytic_z(x_c=x_c, y_c=y_c, z_c=z_c, x_1=x_1, x_2=x_2, y_1=y_1,   &
                                                        y_2=y_2, A=self%coil%A(n), amplitude=self%coil%coil_amplitude(n), &
                                                        sigma=sigma, n=n, adjust_amplitude=.true.)
   else
      self%coil%coil_amplitude(n) = self%coil%coil_amplitude(1)
      call compute_coil_current_density_flux_analytic_z(x_c=x_c, y_c=y_c, z_c=z_c, x_1=x_1, x_2=x_2, y_1=y_1,   &
                                                        y_2=y_2, A=self%coil%A(n), amplitude=self%coil%coil_amplitude(n), &
                                                        sigma=sigma, n=n, adjust_amplitude=.false.)
   endif

   endassociate

   contains

      subroutine compute_coil_current_density_flux_analytic_z(x_c, y_c, z_c, x_1, x_2, y_1, y_2, sigma, &
                                                               A, amplitude, n, adjust_amplitude)
      !< Analytic amplitude correction for a rectangular coil normal to z.
      !< The current is evaluated as the flux of Jx through a section cutting the side y = y_1.
      real(R8P),    intent(in)    :: x_c              !< Coil center x-coordinate.
      real(R8P),    intent(in)    :: y_c              !< Coil center y-coordinate.
      real(R8P),    intent(in)    :: z_c              !< Coil center z-coordinate.
      real(R8P),    intent(in)    :: x_1              !< Lower x-boundary of the rectangular coil.
      real(R8P),    intent(in)    :: x_2              !< Upper x-boundary of the rectangular coil.
      real(R8P),    intent(in)    :: y_1              !< Lower y-boundary of the rectangular coil.
      real(R8P),    intent(in)    :: y_2              !< Upper y-boundary of the rectangular coil.
      real(R8P),    intent(in)    :: sigma            !< Coil Gaussian smoothing width.
      real(R8P),    intent(in)    :: A                !< Coil amplitude, target.
      real(R8P),    intent(inout) :: amplitude        !< Corrected amplitude.
      integer(I4P), intent(in)    :: n                !< Coil number.
      logical,      intent(in)    :: adjust_amplitude !< If true, correct the coil amplitude.
      real(R8P), parameter        :: alpha = 3.5_R8P  !< Half-width of the integration interval in sigma units.
      real(R8P)                   :: z_l              !< Lower z-boundary of the coil thickness window.
      real(R8P)                   :: z_u              !< Upper z-boundary of the coil thickness window.
      real(R8P)                   :: x_s              !< x-coordinate of the flux-section plane.
      real(R8P)                   :: y_s              !< y-coordinate of the flux-section center.
      real(R8P)                   :: z_s              !< z-coordinate of the flux-section center.
      real(R8P)                   :: y_min            !< Lower y-integration boundary.
      real(R8P)                   :: y_max            !< Upper y-integration boundary.
      real(R8P)                   :: z_min            !< Lower z-integration boundary.
      real(R8P)                   :: z_max            !< Upper z-integration boundary.
      real(R8P)                   :: e_x1_xs          !< E_x1 evaluated at x_s.
      real(R8P)                   :: e_x2_xs          !< E_x2 evaluated at x_s.
      real(R8P)                   :: w_x_xs           !< W_x evaluated at x_s.
      real(R8P)                   :: dw_y             !< Difference W_y(y_max)-W_y(y_min).
      real(R8P)                   :: de_y1            !< Difference E_y1(y_max)-E_y1(y_min).
      real(R8P)                   :: de_y2            !< Difference E_y2(y_max)-E_y2(y_min).
      real(R8P)                   :: i_z              !< Integral of W_z over the z-integration interval.
      real(R8P)                   :: i_x_1            !< Flux contribution associated with A_1.
      real(R8P)                   :: i_x_2            !< Flux contribution associated with A_2.
      real(R8P)                   :: i_x_3            !< Flux contribution associated with A_3.
      real(R8P)                   :: i_x_4            !< Flux contribution associated with A_4.
      real(R8P)                   :: flux_unit        !< Analytic current flux for unit amplitude.
      real(R8P)                   :: flux             !< Analytic current flux before amplitude correction.
      real(R8P)                   :: correction       !< Amplitude correction factor.
      real(R8P)                   :: target_current   !< Target input current.

      target_current = A

      ! Set z-window boundaries:
      !
      ! W_z = 1/2 * [E_z_l - E_z_u].
      z_l = z_c - sigma
      z_u = z_c + sigma

      ! Set flux section for side 1, where the current is directed along +x.
      !
      ! Side 1:
      !   y = y_1
      !   J_x > 0 for verse = +1
      !
      ! The section is normal to x and centered at (x_c, y_1, z_c).
      x_s = x_c
      y_s = y_1
      z_s = z_c

      y_min = y_s - alpha * sigma
      y_max = y_s + alpha * sigma
      z_min = z_s - alpha * sigma
      z_max = z_s + alpha * sigma

      ! Integral of W_z over z:
      !
      ! I_z = 1/2 * [ P_zl(zmax) - P_zl(zmin)
      !              -P_zu(zmax) + P_zu(zmin) ]
      i_z = 0.5_R8P * (                                              &
            erf_primitive_function(s=z_max, mu=z_l, sigma=sigma)     &
          - erf_primitive_function(s=z_min, mu=z_l, sigma=sigma)     &
          - erf_primitive_function(s=z_max, mu=z_u, sigma=sigma)     &
          + erf_primitive_function(s=z_min, mu=z_u, sigma=sigma))

      ! Evaluate x-dependent functions at the flux-section plane.
      e_x1_xs = erf_function(s=x_s, mu=x_1, sigma=sigma)
      e_x2_xs = erf_function(s=x_s, mu=x_2, sigma=sigma)

      w_x_xs = tangential_window(s     = x_s,     &
                                 smin  = x_1,     &
                                 smax  = x_2,     &
                                 sigma = sigma)

      ! Compute Delta W_y = W_y(y_max) - W_y(y_min).
      dw_y = tangential_window(s     = y_max,  &
                               smin  = y_1,    &
                               smax  = y_2,    &
                               sigma = sigma)  &
           - tangential_window(s     = y_min,  &
                               smin  = y_1,    &
                               smax  = y_2,    &
                               sigma = sigma)

      ! Compute Delta E_y1 = E_y1(y_max) - E_y1(y_min).
      de_y1 = erf_function(s=y_max, mu=y_1, sigma=sigma) &
            - erf_function(s=y_min, mu=y_1, sigma=sigma)

      ! Compute Delta E_y2 = E_y2(y_max) - E_y2(y_min).
      de_y2 = erf_function(s=y_max, mu=y_2, sigma=sigma) &
            - erf_function(s=y_min, mu=y_2, sigma=sigma)

      ! Contributions to Ix from the four vector-potential components:
      !
      ! A_1 =  E_y1(y) W_x(x) W_z(z)
      ! A_2 = -E_x2(x) W_y(y) W_z(z)
      ! A_3 = -E_y2(y) W_x(x) W_z(z)
      ! A_4 =  E_x1(x) W_y(y) W_z(z)
      !
      ! Since A = A_z e_z:
      !
      ! J_x = dA_z/dy
      !
      ! The flux is evaluated through a section normal to x.
      i_x_1 =  w_x_xs  * i_z * de_y1
      i_x_2 = -e_x2_xs * i_z * dw_y
      i_x_3 = -w_x_xs  * i_z * de_y2
      i_x_4 =  e_x1_xs * i_z * dw_y

      flux_unit = abs(i_x_1 + i_x_2 + i_x_3 + i_x_4)
      flux      = target_current * flux_unit

      if (adjust_amplitude) then

         correction = 1.0_R8P / flux_unit

         print '(A)', mpih%myrankstr//'Current before correction: '//trim(str(flux))
         print '(A)', mpih%myrankstr//trim(str(A))//' Amplitude A('//trim(str(n))//') before correction'
         print '(A)', mpih%myrankstr//'Amplitude scaling factor: '//trim(str(correction))

         amplitude = target_current * correction

         print '(A)', mpih%myrankstr//trim(str(amplitude))//' Amplitude A('//trim(str(n))//') after correction'
         print '(A)', mpih%myrankstr//'Final coil current '//trim(str(n))//': '//trim(str(amplitude*flux_unit))

      else

         print '(A)', mpih%myrankstr//'Amplitude A('//trim(str(n))//') not corrected: '//trim(str(amplitude))
         print '(A)', mpih%myrankstr//'Final coil current '//trim(str(n))//': '//trim(str(amplitude*flux_unit))

      endif

      endsubroutine compute_coil_current_density_flux_analytic_z
   endsubroutine set_rectangular_coil_z

   subroutine set_circular_coil_x(self, n, verse)
   !< Set circular coil with normal direction parallel to x.
   class(prism_common_object), intent(inout) :: self                    !< Cpu object.
   integer(I4P),               intent(in)    :: n                       !< Coil number.
   real(R8P),                  intent(in)    :: verse                   !< Coil normal direction, +1=+x, -1=-x.
   real(R8P), allocatable                    :: A(:,:,:,:,:)            !< Total coil vector potential field.
   real(R8P), allocatable                    :: J_vec_buffer(:,:,:,:,:) !< Buffer variable for self%coil%J_vec.
   real(R8P)                                 :: A_c                     !< Circular vector potential contribution.
   real(R8P)                                 :: cell_coord(3)           !< Cell-center coordinate vector.
   real(R8P)                                 :: r                       !< Radial distance in the y-z plane.
   real(R8P)                                 :: F_r                     !< Radial error-function profile.
   real(R8P)                                 :: W_x                     !< Coil thickness window function along x.
   integer(I4P)                              :: b,i,j,k                 !< Counters.

   ! Associate grid, field, and coil data.
   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,         &
            nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                     &
            y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), r_c=>self%coil%r_coil(n),      &
            normal=>self%coil%normal(n), nb=>field%nb, x_cell=>field%x_cell,             &
            y_cell=>field%y_cell, z_cell=>field%z_cell, sigma=>self%coil%sigma(n),       &
            J_vec=>self%coil%J_vec, hs=>self%fdv_half_stencil)

   allocate(A(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   A(:,:,:,:,:) = 0.0_R8P

   allocate(J_vec_buffer(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   J_vec_buffer(:,:,:,:,:) = 0.0_R8P

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]

               ! Circular coil in the y-z plane.
               r = sqrt((cell_coord(2)-y_c)**2 + (cell_coord(3)-z_c)**2)

               ! Radial profile. The minus sign gives positive circulation for verse = +1.
               F_r = -erf_function(s=r, mu=r_c, sigma=sigma)

               ! Thickness window along the coil normal direction.
               W_x = tangential_window(s=cell_coord(1), smin=x_c-sigma, smax=x_c+sigma, sigma=sigma)

               A_c = F_r*W_x

               ! Store vector potential component normal to the circular coil plane.
               A(1,i,j,k,b) = A_c
            enddo
         enddo
      enddo
   enddo

   call self%compute_curl(hs=hs, ivar=1_I4P, q=A, curl=J_vec_buffer)

   J_vec_buffer = J_vec_buffer * verse

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (maxval(abs(J_vec_buffer(:,i,j,k,b))) < 1e-12_R8P) then
                  J_vec_buffer(:,i,j,k,b) = 0.0_R8P
               endif
            enddo
         enddo
      enddo
   enddo

   J_vec(1:3,:,:,:,:,n) = J_vec(1:3,:,:,:,:,n) + J_vec_buffer

   if (n == 1_I4P) then
      call compute_coil_current_density_flux_analytic_circular_x(x_c=x_c, y_c=y_c, z_c=z_c, r_c=r_c, &
                                                                 A=self%coil%A(n), sigma=sigma, n=n,      &
                                                                 adjust_amplitude=.true.)
   else
      self%coil%A(n) = self%coil%A(1)
      call compute_coil_current_density_flux_analytic_circular_x(x_c=x_c, y_c=y_c, z_c=z_c, r_c=r_c, &
                                                                 A=self%coil%A(n), sigma=sigma, n=n,      &
                                                                 adjust_amplitude=.false.)
   endif

   endassociate

   contains

      subroutine compute_coil_current_density_flux_analytic_circular_x(x_c, y_c, z_c, r_c, sigma, A, n, adjust_amplitude)
      !< Analytic amplitude correction for a circular coil normal to x.
      !< The current is evaluated as the flux of Jy through the lower point z = z_c - r_c.
      real(R8P),    intent(in)    :: x_c              !< Coil center x-coordinate.
      real(R8P),    intent(in)    :: y_c              !< Coil center y-coordinate.
      real(R8P),    intent(in)    :: z_c              !< Coil center z-coordinate.
      real(R8P),    intent(in)    :: r_c              !< Circular coil radius.
      real(R8P),    intent(in)    :: sigma            !< Coil Gaussian smoothing width.
      real(R8P),    intent(inout) :: A                !< Coil amplitude, corrected in place.
      integer(I4P), intent(in)    :: n                !< Coil number.
      logical,      intent(in)    :: adjust_amplitude !< If true, correct the coil amplitude.
      real(R8P), parameter        :: alpha = 3.5_R8P  !< Half-width of the integration interval in sigma units.
      real(R8P)                   :: x_l              !< Lower x-boundary of the coil thickness window.
      real(R8P)                   :: x_u              !< Upper x-boundary of the coil thickness window.
      real(R8P)                   :: x_s              !< x-coordinate of the flux-section center.
      real(R8P)                   :: z_s              !< z-coordinate of the lower flux-section center.
      real(R8P)                   :: x_1              !< Lower x-integration boundary.
      real(R8P)                   :: x_2              !< Upper x-integration boundary.
      real(R8P)                   :: r_1              !< Lower radial integration boundary.
      real(R8P)                   :: r_2              !< Upper radial integration boundary.
      real(R8P)                   :: de_r             !< Difference E_R(r_2)-E_R(r_1).
      real(R8P)                   :: i_x              !< Integral of W_x over the x-integration interval.
      real(R8P)                   :: i_y_c            !< Flux contribution associated with the circular current profile.
      real(R8P)                   :: flux_unit        !< Analytic current flux for unit amplitude.
      real(R8P)                   :: flux             !< Analytic current flux before amplitude correction.
      real(R8P)                   :: correction       !< Amplitude correction factor.
      real(R8P)                   :: target_current   !< Target input current.

      target_current = A

      if (r_c <= alpha*sigma) then
         error stop 'compute_coil_current_density_flux_analytic_circular_x: r_c must be larger than alpha*sigma'
      endif

      ! Set x-window boundaries:
      ! W_x = 1/2 * [E_x_l - E_x_u].
      x_l = x_c - sigma
      x_u = x_c + sigma

      ! Set flux section for the lower point of the circular coil.
      x_s = x_c
      z_s = z_c - r_c

      x_1 = x_s - alpha * sigma
      x_2 = x_s + alpha * sigma

      r_1 = r_c - alpha * sigma
      r_2 = r_c + alpha * sigma

      ! Integral of W_x over x.
      i_x = 0.5_R8P * (                                        &
            erf_primitive_function(s=x_2, mu=x_l, sigma=sigma) &
          - erf_primitive_function(s=x_1, mu=x_l, sigma=sigma) &
          - erf_primitive_function(s=x_2, mu=x_u, sigma=sigma) &
          + erf_primitive_function(s=x_1, mu=x_u, sigma=sigma))

      ! Radial Gaussian flux contribution.
      de_r = erf_function(s=r_2, mu=r_c, sigma=sigma) &
           - erf_function(s=r_1, mu=r_c, sigma=sigma)

      i_y_c = i_x * de_r

      flux_unit = abs(i_y_c)
      flux      = target_current * flux_unit

      if (adjust_amplitude) then

         correction = 1.0_R8P / flux_unit

         print '(A)', mpih%myrankstr//'Current before correction: '//trim(str(flux))
         print '(A)', mpih%myrankstr//trim(str(A))//' Amplitude A('//trim(str(n))//') before correction'
         print '(A)', mpih%myrankstr//'Amplitude scaling factor: '//trim(str(correction))

         A = target_current * correction

         print '(A)', mpih%myrankstr//trim(str(A))//' Amplitude A('//trim(str(n))//') after correction'
         print '(A)', mpih%myrankstr//'Final coil current '//trim(str(n))//': '//trim(str(A*flux_unit))

      else

         print '(A)', mpih%myrankstr//'Amplitude A('//trim(str(n))//') not corrected: '//trim(str(A))
         print '(A)', mpih%myrankstr//'Final coil current '//trim(str(n))//': '//trim(str(flux))

      endif

      endsubroutine compute_coil_current_density_flux_analytic_circular_x

   endsubroutine set_circular_coil_x

   subroutine set_circular_coil_y(self, n, verse)
   !< Set circular coil with normal direction parallel to y.
   class(prism_common_object), intent(inout) :: self                    !< Cpu object.
   integer(I4P),               intent(in)    :: n                       !< Coil number.
   real(R8P),                  intent(in)    :: verse                   !< Coil normal direction, +1=+y, -1=-y.
   real(R8P), allocatable                    :: A(:,:,:,:,:)            !< Total coil vector potential field.
   real(R8P), allocatable                    :: J_vec_buffer(:,:,:,:,:) !< Buffer variable for self%coil%J_vec.
   real(R8P)                                 :: A_c                     !< Circular vector potential contribution.
   real(R8P)                                 :: cell_coord(3)           !< Cell-center coordinate vector.
   real(R8P)                                 :: r                       !< Radial distance in the x-z plane.
   real(R8P)                                 :: F_r                     !< Radial error-function profile.
   real(R8P)                                 :: W_y                     !< Coil thickness window function along y.
   integer(I4P)                              :: b,i,j,k                 !< Counters.

   ! Associate grid, field, and coil data.
   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,         &
            nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                     &
            y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), r_c=>self%coil%r_coil(n),      &
            normal=>self%coil%normal(n), nb=>field%nb, x_cell=>field%x_cell,             &
            y_cell=>field%y_cell, z_cell=>field%z_cell, sigma=>self%coil%sigma(n),       &
            J_vec=>self%coil%J_vec, hs=>self%fdv_half_stencil)

   allocate(A(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   A(:,:,:,:,:) = 0.0_R8P

   allocate(J_vec_buffer(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   J_vec_buffer(:,:,:,:,:) = 0.0_R8P

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]

               ! Circular coil in the x-z plane.
               r = sqrt((cell_coord(1)-x_c)**2 + (cell_coord(3)-z_c)**2)

               ! Radial profile. The minus sign gives positive circulation for verse = +1.
               F_r = -erf_function(s=r, mu=r_c, sigma=sigma)

               ! Thickness window along the coil normal direction.
               W_y = tangential_window(s=cell_coord(2), smin=y_c-sigma, smax=y_c+sigma, sigma=sigma)

               A_c = F_r*W_y

               ! Store vector potential component normal to the circular coil plane.
               A(2,i,j,k,b) = A_c
            enddo
         enddo
      enddo
   enddo

   call self%compute_curl(hs=hs, ivar=1_I4P, q=A, curl=J_vec_buffer)

   J_vec_buffer = J_vec_buffer * verse

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (maxval(abs(J_vec_buffer(:,i,j,k,b))) < 1e-12_R8P) then
                  J_vec_buffer(:,i,j,k,b) = 0.0_R8P
               endif
            enddo
         enddo
      enddo
   enddo

   J_vec(1:3,:,:,:,:,n) = J_vec(1:3,:,:,:,:,n) + J_vec_buffer

   if (n == 1_I4P) then
      call compute_coil_current_density_flux_analytic_circular_y(x_c=x_c, y_c=y_c, z_c=z_c, r_c=r_c, &
                                                                 A=self%coil%A(n), sigma=sigma, n=n,      &
                                                                 adjust_amplitude=.true.)
   else
      self%coil%A(n) = self%coil%A(1)
      call compute_coil_current_density_flux_analytic_circular_y(x_c=x_c, y_c=y_c, z_c=z_c, r_c=r_c, &
                                                                 A=self%coil%A(n), sigma=sigma, n=n,      &
                                                                 adjust_amplitude=.false.)
   endif

   endassociate

   contains

      subroutine compute_coil_current_density_flux_analytic_circular_y(x_c, y_c, z_c, r_c, sigma, A, n, adjust_amplitude)
      !< Analytic amplitude correction for a circular coil normal to y.
      !< The current is evaluated as the flux of Jz through the left point x = x_c - r_c.
      real(R8P),    intent(in)    :: x_c              !< Coil center x-coordinate.
      real(R8P),    intent(in)    :: y_c              !< Coil center y-coordinate.
      real(R8P),    intent(in)    :: z_c              !< Coil center z-coordinate.
      real(R8P),    intent(in)    :: r_c              !< Circular coil radius.
      real(R8P),    intent(in)    :: sigma            !< Coil Gaussian smoothing width.
      real(R8P),    intent(inout) :: A                !< Coil amplitude, corrected in place.
      integer(I4P), intent(in)    :: n                !< Coil number.
      logical,      intent(in)    :: adjust_amplitude !< If true, correct the coil amplitude.
      real(R8P), parameter        :: alpha = 3.5_R8P  !< Half-width of the integration interval in sigma units.
      real(R8P)                   :: y_l              !< Lower y-boundary of the coil thickness window.
      real(R8P)                   :: y_u              !< Upper y-boundary of the coil thickness window.
      real(R8P)                   :: x_s              !< x-coordinate of the left flux-section center.
      real(R8P)                   :: y_s              !< y-coordinate of the flux-section center.
      real(R8P)                   :: y_1              !< Lower y-integration boundary.
      real(R8P)                   :: y_2              !< Upper y-integration boundary.
      real(R8P)                   :: r_1              !< Lower radial integration boundary.
      real(R8P)                   :: r_2              !< Upper radial integration boundary.
      real(R8P)                   :: de_r             !< Difference E_R(r_2)-E_R(r_1).
      real(R8P)                   :: i_y              !< Integral of W_y over the y-integration interval.
      real(R8P)                   :: i_z_c            !< Flux contribution associated with the circular current profile.
      real(R8P)                   :: flux_unit        !< Analytic current flux for unit amplitude.
      real(R8P)                   :: flux             !< Analytic current flux before amplitude correction.
      real(R8P)                   :: correction       !< Amplitude correction factor.
      real(R8P)                   :: target_current   !< Target input current.

      target_current = A

      if (r_c <= alpha*sigma) then
         error stop 'compute_coil_current_density_flux_analytic_circular_y: r_c must be larger than alpha*sigma'
      endif

      ! Set y-window boundaries:
      ! W_y = 1/2 * [E_y_l - E_y_u].
      y_l = y_c - sigma
      y_u = y_c + sigma

      ! Set flux section for the left point of the circular coil.
      x_s = x_c - r_c
      y_s = y_c

      y_1 = y_s - alpha * sigma
      y_2 = y_s + alpha * sigma

      r_1 = r_c - alpha * sigma
      r_2 = r_c + alpha * sigma

      ! Integral of W_y over y.
      i_y = 0.5_R8P * (                                        &
            erf_primitive_function(s=y_2, mu=y_l, sigma=sigma) &
          - erf_primitive_function(s=y_1, mu=y_l, sigma=sigma) &
          - erf_primitive_function(s=y_2, mu=y_u, sigma=sigma) &
          + erf_primitive_function(s=y_1, mu=y_u, sigma=sigma))

      ! Radial Gaussian flux contribution.
      de_r = erf_function(s=r_2, mu=r_c, sigma=sigma) &
           - erf_function(s=r_1, mu=r_c, sigma=sigma)

      i_z_c = i_y * de_r

      flux_unit = abs(i_z_c)
      flux      = target_current * flux_unit

      if (adjust_amplitude) then

         correction = 1.0_R8P / flux_unit

         print '(A)', mpih%myrankstr//'Current before correction: '//trim(str(flux))
         print '(A)', mpih%myrankstr//trim(str(A))//' Amplitude A('//trim(str(n))//') before correction'
         print '(A)', mpih%myrankstr//'Amplitude scaling factor: '//trim(str(correction))

         A = target_current * correction

         print '(A)', mpih%myrankstr//trim(str(A))//' Amplitude A('//trim(str(n))//') after correction'
         print '(A)', mpih%myrankstr//'Final coil current '//trim(str(n))//': '//trim(str(A*flux_unit))

      else

         print '(A)', mpih%myrankstr//'Amplitude A('//trim(str(n))//') not corrected: '//trim(str(A))
         print '(A)', mpih%myrankstr//'Final coil current '//trim(str(n))//': '//trim(str(flux))

      endif

      endsubroutine compute_coil_current_density_flux_analytic_circular_y

   endsubroutine set_circular_coil_y

   subroutine set_circular_coil_z(self, n, verse)
   !< Set circular coil with normal direction parallel to z.
   class(prism_common_object), intent(inout) :: self                    !< Cpu object.
   integer(I4P),               intent(in)    :: n                       !< Coil number.
   real(R8P),                  intent(in)    :: verse                   !< Coil normal direction, +1=+z, -1=-z.
   real(R8P), allocatable                    :: A(:,:,:,:,:)            !< Total coil vector potential field.
   real(R8P), allocatable                    :: J_vec_buffer(:,:,:,:,:) !< Buffer variable for self%coil%J_vec.
   real(R8P)                                 :: A_c                     !< Circular vector potential contribution.
   real(R8P)                                 :: cell_coord(3)           !< Cell-center coordinate vector.
   real(R8P)                                 :: r                       !< Radial distance in the x-y plane.
   real(R8P)                                 :: F_r                     !< Radial error-function profile.
   real(R8P)                                 :: W_z                     !< Coil thickness window function along z.
   integer(I4P)                              :: b,i,j,k                 !< Counters.

   ! Associate grid, field, and coil data.
   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,          &
            nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                      &
            y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), r_c=>self%coil%r_coil(n),      &
            normal=>self%coil%normal(n), nb=>field%nb, x_cell=>field%x_cell,             &
            y_cell=>field%y_cell, z_cell=>field%z_cell, sigma=>self%coil%sigma(n),       &
            J_vec=>self%coil%J_vec, hs=>self%fdv_half_stencil)

   allocate(A(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   A(:,:,:,:,:) = 0.0_R8P

   allocate(J_vec_buffer(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   J_vec_buffer(:,:,:,:,:) = 0.0_R8P

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]

               ! Circular coil in the x-y plane.
               r = sqrt((cell_coord(1)-x_c)**2 + (cell_coord(2)-y_c)**2)

               ! Radial profile. The minus sign gives positive circulation for verse = +1.
               F_r = -erf_function(s=r, mu=r_c, sigma=sigma)

               ! Thickness window along the coil normal direction.
               W_z = tangential_window(s=cell_coord(3), smin=z_c-sigma, smax=z_c+sigma, sigma=sigma)

               A_c = F_r*W_z

               ! Store vector potential component normal to the circular coil plane.
               A(3,i,j,k,b) = A_c
            enddo
         enddo
      enddo
   enddo

   call self%compute_curl(hs=hs, ivar=1_I4P, q=A, curl=J_vec_buffer)

   J_vec_buffer = J_vec_buffer * verse

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (maxval(abs(J_vec_buffer(:,i,j,k,b))) < 1e-12_R8P) then
                  J_vec_buffer(:,i,j,k,b) = 0.0_R8P
               endif
            enddo
         enddo
      enddo
   enddo

   J_vec(1:3,:,:,:,:,n) = J_vec(1:3,:,:,:,:,n) + J_vec_buffer

   if (n == 1_I4P) then
      call compute_coil_current_density_flux_analytic_circular_z(x_c=x_c, y_c=y_c, z_c=z_c, r_c=r_c, &
                                                                 A=self%coil%A(n), sigma=sigma, n=n,      &
                                                                 adjust_amplitude=.true.)
   else
      self%coil%A(n) = self%coil%A(1)
      call compute_coil_current_density_flux_analytic_circular_z(x_c=x_c, y_c=y_c, z_c=z_c, r_c=r_c, &
                                                                 A=self%coil%A(n), sigma=sigma, n=n,      &
                                                                 adjust_amplitude=.false.)
   endif

   endassociate

   contains

      subroutine compute_coil_current_density_flux_analytic_circular_z(x_c, y_c, z_c, r_c, sigma, A, n, adjust_amplitude)
      !< Analytic amplitude correction for a circular coil normal to z.
      !< The current is evaluated as the flux of Jy through the left point x = x_c - r_c.
      real(R8P),    intent(in)    :: x_c              !< Coil center x-coordinate.
      real(R8P),    intent(in)    :: y_c              !< Coil center y-coordinate.
      real(R8P),    intent(in)    :: z_c              !< Coil center z-coordinate.
      real(R8P),    intent(in)    :: r_c              !< Circular coil radius.
      real(R8P),    intent(in)    :: sigma            !< Coil Gaussian smoothing width.
      real(R8P),    intent(inout) :: A                !< Coil amplitude, corrected in place.
      integer(I4P), intent(in)    :: n                !< Coil number.
      logical,      intent(in)    :: adjust_amplitude !< If true, correct the coil amplitude.
      real(R8P), parameter        :: alpha = 3.5_R8P  !< Half-width of the integration interval in sigma units.
      real(R8P)                   :: z_l              !< Lower z-boundary of the coil thickness window.
      real(R8P)                   :: z_u              !< Upper z-boundary of the coil thickness window.
      real(R8P)                   :: x_s              !< x-coordinate of the left flux-section center.
      real(R8P)                   :: z_s              !< z-coordinate of the flux-section center.
      real(R8P)                   :: z_1              !< Lower z-integration boundary.
      real(R8P)                   :: z_2              !< Upper z-integration boundary.
      real(R8P)                   :: r_1              !< Lower radial integration boundary.
      real(R8P)                   :: r_2              !< Upper radial integration boundary.
      real(R8P)                   :: de_r             !< Difference E_R(r_2)-E_R(r_1).
      real(R8P)                   :: i_z              !< Integral of W_z over the z-integration interval.
      real(R8P)                   :: i_y_c            !< Flux contribution associated with the circular current profile.
      real(R8P)                   :: flux_unit        !< Analytic current flux for unit amplitude.
      real(R8P)                   :: flux             !< Analytic current flux before amplitude correction.
      real(R8P)                   :: correction       !< Amplitude correction factor.
      real(R8P)                   :: target_current   !< Target input current.

      target_current = A

      if (r_c <= alpha*sigma) then
         error stop 'compute_coil_current_density_flux_analytic_circular_z: r_c must be larger than alpha*sigma'
      endif

      ! Set z-window boundaries:
      ! W_z = 1/2 * [E_z_l - E_z_u].
      z_l = z_c - sigma
      z_u = z_c + sigma

      ! Set flux section for the left point of the circular coil.
      x_s = x_c - r_c
      z_s = z_c

      z_1 = z_s - alpha * sigma
      z_2 = z_s + alpha * sigma

      r_1 = r_c - alpha * sigma
      r_2 = r_c + alpha * sigma

      ! Integral of W_z over z.
      i_z = 0.5_R8P * (                                        &
            erf_primitive_function(s=z_2, mu=z_l, sigma=sigma) &
          - erf_primitive_function(s=z_1, mu=z_l, sigma=sigma) &
          - erf_primitive_function(s=z_2, mu=z_u, sigma=sigma) &
          + erf_primitive_function(s=z_1, mu=z_u, sigma=sigma))

      ! Radial Gaussian flux contribution.
      de_r = erf_function(s=r_2, mu=r_c, sigma=sigma) &
           - erf_function(s=r_1, mu=r_c, sigma=sigma)

      ! At x = x_c-r_c, the positive-normal circulation has Jy < 0.
      i_y_c = -i_z * de_r

      flux_unit = abs(i_y_c)
      flux      = target_current * flux_unit

      if (adjust_amplitude) then

         correction = 1.0_R8P / flux_unit

         print '(A)', mpih%myrankstr//'Current before correction: '//trim(str(flux))
         print '(A)', mpih%myrankstr//trim(str(A))//' Amplitude A('//trim(str(n))//') before correction'
         print '(A)', mpih%myrankstr//'Amplitude scaling factor: '//trim(str(correction))

         A = target_current * correction

         print '(A)', mpih%myrankstr//trim(str(A))//' Amplitude A('//trim(str(n))//') after correction'
         print '(A)', mpih%myrankstr//'Final coil current '//trim(str(n))//': '//trim(str(A*flux_unit))

      else

         print '(A)', mpih%myrankstr//'Amplitude A('//trim(str(n))//') not corrected: '//trim(str(A))
         print '(A)', mpih%myrankstr//'Final coil current '//trim(str(n))//': '//trim(str(flux))

      endif

      endsubroutine compute_coil_current_density_flux_analytic_circular_z

   endsubroutine set_circular_coil_z

   subroutine set_solenoid_x(self, n, verse)
   !< Set solenoid with axis direction parallel to x.
   class(prism_common_object), intent(inout) :: self                    !< Cpu object.
   integer(I4P),               intent(in)    :: n                       !< Coil number.
   real(R8P),                  intent(in)    :: verse                   !< Solenoid axis direction, +1=+x, -1=-x.
   real(R8P), allocatable                    :: A(:,:,:,:,:)            !< Total solenoid vector potential field.
   real(R8P), allocatable                    :: J_vec_buffer(:,:,:,:,:) !< Buffer variable for self%coil%J_vec.
   real(R8P)                                 :: cell_coord(3)           !< Cell-center coordinate vector.
   real(R8P)                                 :: rho                     !< Radial distance in the y-z plane.
   real(R8P)                                 :: F_n                     !< Radial error-function profile.
   real(R8P)                                 :: W_x                     !< Axial window function along x.
   integer(I4P)                              :: b,i,j,k                 !< Counters.

   ! Associate grid, field, and coil data.
   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,         &
            nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                     &
            y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), r_coil=>self%coil%r_coil(n),  &
            l_sol=>self%coil%l_solenoid(n), windings=>self%coil%windings(n),                 &
            nb=>field%nb, x_cell=>field%x_cell, y_cell=>field%y_cell,              &
            z_cell=>field%z_cell, sigma=>self%coil%sigma(n), J_vec=>self%coil%J_vec,         &
            hs=>self%fdv_half_stencil)

   allocate(A(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   A(:,:,:,:,:) = 0.0_R8P

   allocate(J_vec_buffer(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   J_vec_buffer(:,:,:,:,:) = 0.0_R8P

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]

               rho = sqrt((cell_coord(2)-y_c)**2 + (cell_coord(3)-z_c)**2)

               F_n = erf_function(s=rho, mu=r_coil, sigma=sigma)

               W_x = tangential_window(s     = cell_coord(1),     &
                                        smin  = x_c-0.5_R8P*l_sol, &
                                        smax  = x_c+0.5_R8P*l_sol, &
                                        sigma = sigma)

               A(1,i,j,k,b) = -F_n*W_x
            enddo
         enddo
      enddo
   enddo

   call self%compute_curl(hs=hs, ivar=1_I4P, q=A, curl=J_vec_buffer)

   J_vec_buffer = J_vec_buffer * verse

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (maxval(abs(J_vec_buffer(:,i,j,k,b))) < 1e-12_R8P) then
                  J_vec_buffer(:,i,j,k,b) = 0.0_R8P
               endif
            enddo
         enddo
      enddo
   enddo

   J_vec(1:3,:,:,:,:,n) = J_vec(1:3,:,:,:,:,n) + J_vec_buffer

   if (n == 1_I4P) then
      call compute_solenoid_current_density_flux_analytic_x(x_c=x_c, y_c=y_c, z_c=z_c, r_coil=r_coil, &
                                                            l_sol=l_sol, windings=windings,           &
                                                            A=self%coil%A(n), sigma=sigma, n=n,            &
                                                            adjust_amplitude=.true.)
   else
      self%coil%A(n) = self%coil%A(1)
      call compute_solenoid_current_density_flux_analytic_x(x_c=x_c, y_c=y_c, z_c=z_c, r_coil=r_coil, &
                                                            l_sol=l_sol, windings=windings,           &
                                                            A=self%coil%A(n), sigma=sigma, n=n,            &
                                                            adjust_amplitude=.false.)
   endif

   endassociate

   contains

      subroutine compute_solenoid_current_density_flux_analytic_x(x_c, y_c, z_c, r_coil, l_sol, windings, &
                                                                   sigma, A, n, adjust_amplitude)
      !< Analytic amplitude correction for a solenoid with axis parallel to x.
      !< The current is evaluated as the flux of Jz through the lower radial section y = y_c-r_coil.
      real(R8P),    intent(in)    :: x_c              !< Solenoid center x-coordinate.
      real(R8P),    intent(in)    :: y_c              !< Solenoid center y-coordinate.
      real(R8P),    intent(in)    :: z_c              !< Solenoid center z-coordinate.
      real(R8P),    intent(in)    :: r_coil           !< Solenoid radius.
      real(R8P),    intent(in)    :: l_sol            !< Solenoid length.
      real(R8P),    intent(in)    :: windings         !< Number of solenoid windings.
      real(R8P),    intent(in)    :: sigma            !< Solenoid Gaussian smoothing width.
      real(R8P),    intent(inout) :: A                !< Solenoid amplitude, corrected in place.
      integer(I4P), intent(in)    :: n                !< Coil number.
      logical,      intent(in)    :: adjust_amplitude !< If true, correct the solenoid amplitude.
      real(R8P),    parameter     :: alpha = 3.5_R8P  !< Half-width of the radial integration interval in sigma units.
      real(R8P)                   :: x_l              !< Lower x-boundary of the solenoid axial window.
      real(R8P)                   :: x_u              !< Upper x-boundary of the solenoid axial window.
      real(R8P)                   :: x_1              !< Lower x-integration boundary.
      real(R8P)                   :: x_2              !< Upper x-integration boundary.
      real(R8P)                   :: rho_1            !< Lower radial integration boundary.
      real(R8P)                   :: rho_2            !< Upper radial integration boundary.
      real(R8P)                   :: i_x              !< Integral of W_x over the x-integration interval.
      real(R8P)                   :: de_rho           !< Difference E_R(rho_2)-E_R(rho_1).
      real(R8P)                   :: i_z              !< Analytic flux contribution of Jz.
      real(R8P)                   :: flux_unit        !< Analytic current flux for unit amplitude.
      real(R8P)                   :: flux             !< Analytic current flux before amplitude correction.
      real(R8P)                   :: correction       !< Amplitude correction factor.
      real(R8P)                   :: target_current   !< Target total current, equal to input A*windings.

      x_l = x_c - 0.5_R8P*l_sol
      x_u = x_c + 0.5_R8P*l_sol

      x_1   = x_l - alpha*sigma
      x_2   = x_u + alpha*sigma
      rho_1 = r_coil - alpha*sigma
      rho_2 = r_coil + alpha*sigma

      i_x = 0.5_R8P * (                                        &
            erf_primitive_function(s=x_2, mu=x_l, sigma=sigma) &
          - erf_primitive_function(s=x_1, mu=x_l, sigma=sigma) &
          - erf_primitive_function(s=x_2, mu=x_u, sigma=sigma) &
          + erf_primitive_function(s=x_1, mu=x_u, sigma=sigma))

      de_rho = erf_function(s=rho_2, mu=r_coil, sigma=sigma) &
             - erf_function(s=rho_1, mu=r_coil, sigma=sigma)

      ! For A_x = -E_R(rho) W_x(x), at y = y_c-r_coil, Jz is negative.
      i_z = -i_x*de_rho

      flux_unit = abs(i_z)

      if (flux_unit <= tiny(1.0_R8P)) then
         error stop 'compute_solenoid_current_density_flux_analytic_x: null analytic unit flux'
      endif

      if (adjust_amplitude) then

         target_current = A * windings
         flux           = target_current * flux_unit
         correction     = 1.0_R8P / flux_unit

         print '(A)', mpih%myrankstr//'Current before correction: '//trim(str(flux))
         print '(A)', mpih%myrankstr//trim(str(target_current))//' Target current A('//trim(str(n))//')*N before correction'
         print '(A)', mpih%myrankstr//'Amplitude scaling factor: '//trim(str(correction))

         A = target_current * correction

         print '(A)', mpih%myrankstr//trim(str(A))//' Amplitude A('//trim(str(n))//') after correction'
         print '(A)', mpih%myrankstr//'Final solenoid current '//trim(str(n))//': '//trim(str(A*flux_unit))

      else

         flux = A * flux_unit

         print '(A)', mpih%myrankstr//'Amplitude A('//trim(str(n))//') not corrected: '//trim(str(A))
         print '(A)', mpih%myrankstr//'Final solenoid current '//trim(str(n))//': '//trim(str(flux))

      endif

      endsubroutine compute_solenoid_current_density_flux_analytic_x

endsubroutine set_solenoid_x

subroutine set_solenoid_y(self, n, verse)
   !< Set solenoid with axis direction parallel to y.
   class(prism_common_object), intent(inout) :: self                    !< Cpu object.
   integer(I4P),               intent(in)    :: n                       !< Coil number.
   real(R8P),                  intent(in)    :: verse                   !< Solenoid axis direction, +1=+y, -1=-y.
   real(R8P), allocatable                    :: A(:,:,:,:,:)            !< Total solenoid vector potential field.
   real(R8P), allocatable                    :: J_vec_buffer(:,:,:,:,:) !< Buffer variable for self%coil%J_vec.
   real(R8P)                                 :: cell_coord(3)           !< Cell-center coordinate vector.
   real(R8P)                                 :: rho                     !< Radial distance in the x-z plane.
   real(R8P)                                 :: F_n                     !< Radial error-function profile.
   real(R8P)                                 :: W_y                     !< Axial window function along y.
   integer(I4P)                              :: b,i,j,k                 !< Counters.

   ! Associate grid, field, and coil data.
   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,         &
            nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                     &
            y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), r_coil=>self%coil%r_coil(n),  &
            l_sol=>self%coil%l_solenoid(n), windings=>self%coil%windings(n),                 &
            nb=>field%nb, x_cell=>field%x_cell, y_cell=>field%y_cell,              &
            z_cell=>field%z_cell, sigma=>self%coil%sigma(n), J_vec=>self%coil%J_vec,         &
            hs=>self%fdv_half_stencil)

   allocate(A(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   A(:,:,:,:,:) = 0.0_R8P

   allocate(J_vec_buffer(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   J_vec_buffer(:,:,:,:,:) = 0.0_R8P

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]

               rho = sqrt((cell_coord(1)-x_c)**2 + (cell_coord(3)-z_c)**2)

               F_n = erf_function(s=rho, mu=r_coil, sigma=sigma)

               W_y = tangential_window(s     = cell_coord(2),     &
                                        smin  = y_c-0.5_R8P*l_sol, &
                                        smax  = y_c+0.5_R8P*l_sol, &
                                        sigma = sigma)

               A(2,i,j,k,b) = -F_n*W_y
            enddo
         enddo
      enddo
   enddo

   call self%compute_curl(hs=hs, ivar=1_I4P, q=A, curl=J_vec_buffer)

   J_vec_buffer = J_vec_buffer * verse

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (maxval(abs(J_vec_buffer(:,i,j,k,b))) < 1e-12_R8P) then
                  J_vec_buffer(:,i,j,k,b) = 0.0_R8P
               endif
            enddo
         enddo
      enddo
   enddo

   J_vec(1:3,:,:,:,:,n) = J_vec(1:3,:,:,:,:,n) + J_vec_buffer

   if (n == 1_I4P) then
      call compute_solenoid_current_density_flux_analytic_y(x_c=x_c, y_c=y_c, z_c=z_c, r_coil=r_coil, &
                                                            l_sol=l_sol, windings=windings,           &
                                                            A=self%coil%A(n), sigma=sigma, n=n,            &
                                                            adjust_amplitude=.true.)
   else
      self%coil%A(n) = self%coil%A(1)
      call compute_solenoid_current_density_flux_analytic_y(x_c=x_c, y_c=y_c, z_c=z_c, r_coil=r_coil, &
                                                            l_sol=l_sol, windings=windings,           &
                                                            A=self%coil%A(n), sigma=sigma, n=n,            &
                                                            adjust_amplitude=.false.)
   endif

   endassociate

   contains

      subroutine compute_solenoid_current_density_flux_analytic_y(x_c, y_c, z_c, r_coil, l_sol, windings, &
                                                                  sigma, A, n, adjust_amplitude)
      !< Analytic amplitude correction for a solenoid with axis parallel to y.
      !< The current is evaluated as the flux of Jz through the section x = x_c-r_coil.
      real(R8P),    intent(in)    :: x_c              !< Solenoid center x-coordinate.
      real(R8P),    intent(in)    :: y_c              !< Solenoid center y-coordinate.
      real(R8P),    intent(in)    :: z_c              !< Solenoid center z-coordinate.
      real(R8P),    intent(in)    :: r_coil           !< Solenoid radius.
      real(R8P),    intent(in)    :: l_sol            !< Solenoid length.
      real(R8P),    intent(in)    :: windings         !< Number of solenoid windings.
      real(R8P),    intent(in)    :: sigma            !< Solenoid Gaussian smoothing width.
      real(R8P),    intent(inout) :: A                !< Solenoid amplitude, corrected in place.
      integer(I4P), intent(in)    :: n                !< Coil number.
      logical,      intent(in)    :: adjust_amplitude !< If true, correct the solenoid amplitude.
      real(R8P),    parameter     :: alpha = 3.5_R8P  !< Half-width of the radial integration interval in sigma units.
      real(R8P)                   :: y_l              !< Lower y-boundary of the solenoid axial window.
      real(R8P)                   :: y_u              !< Upper y-boundary of the solenoid axial window.
      real(R8P)                   :: y_1              !< Lower y-integration boundary.
      real(R8P)                   :: y_2              !< Upper y-integration boundary.
      real(R8P)                   :: rho_1            !< Lower radial integration boundary.
      real(R8P)                   :: rho_2            !< Upper radial integration boundary.
      real(R8P)                   :: i_y              !< Integral of W_y over the y-integration interval.
      real(R8P)                   :: de_rho           !< Difference E_R(rho_2)-E_R(rho_1).
      real(R8P)                   :: i_z              !< Analytic flux contribution of Jz.
      real(R8P)                   :: flux_unit        !< Analytic current flux for unit amplitude.
      real(R8P)                   :: flux             !< Analytic current flux before amplitude correction.
      real(R8P)                   :: correction       !< Amplitude correction factor.
      real(R8P)                   :: target_current   !< Target total current, equal to input A*windings.

      y_l = y_c - 0.5_R8P*l_sol
      y_u = y_c + 0.5_R8P*l_sol

      y_1   = y_l - alpha*sigma
      y_2   = y_u + alpha*sigma
      rho_1 = r_coil - alpha*sigma
      rho_2 = r_coil + alpha*sigma

      i_y = 0.5_R8P * (                                        &
            erf_primitive_function(s=y_2, mu=y_l, sigma=sigma) &
          - erf_primitive_function(s=y_1, mu=y_l, sigma=sigma) &
          - erf_primitive_function(s=y_2, mu=y_u, sigma=sigma) &
          + erf_primitive_function(s=y_1, mu=y_u, sigma=sigma))

      de_rho = erf_function(s=rho_2, mu=r_coil, sigma=sigma) &
             - erf_function(s=rho_1, mu=r_coil, sigma=sigma)

      ! For A_y = -E_R(rho) W_y(y), at x = x_c-r_coil, Jz is positive.
      i_z = i_y*de_rho

      flux_unit = abs(i_z)

      if (flux_unit <= tiny(1.0_R8P)) then
         error stop 'compute_solenoid_current_density_flux_analytic_y: null analytic unit flux'
      endif

      if (adjust_amplitude) then

         target_current = A * windings
         flux           = target_current * flux_unit
         correction     = 1.0_R8P / flux_unit

         print '(A)', mpih%myrankstr//'Current before correction: '//trim(str(flux))
         print '(A)', mpih%myrankstr//trim(str(target_current))//' Target current A('//trim(str(n))//')*N before correction'
         print '(A)', mpih%myrankstr//'Amplitude scaling factor: '//trim(str(correction))

         A = target_current * correction

         print '(A)', mpih%myrankstr//trim(str(A))//' Amplitude A('//trim(str(n))//') after correction'
         print '(A)', mpih%myrankstr//'Final solenoid current '//trim(str(n))//': '//trim(str(A*flux_unit))

      else

         flux = A * flux_unit

         print '(A)', mpih%myrankstr//'Amplitude A('//trim(str(n))//') not corrected: '//trim(str(A))
         print '(A)', mpih%myrankstr//'Final solenoid current '//trim(str(n))//': '//trim(str(flux))

      endif

      endsubroutine compute_solenoid_current_density_flux_analytic_y

   endsubroutine set_solenoid_y

   subroutine set_solenoid_z(self, n, verse)
   !< Set solenoid with axis direction parallel to z.
   class(prism_common_object), intent(inout) :: self                    !< Cpu object.
   integer(I4P),               intent(in)    :: n                       !< Coil number.
   real(R8P),                  intent(in)    :: verse                   !< Solenoid axis direction, +1=+z, -1=-z.
   real(R8P), allocatable                    :: A(:,:,:,:,:)            !< Total solenoid vector potential field.
   real(R8P), allocatable                    :: J_vec_buffer(:,:,:,:,:) !< Buffer variable for self%coil%J_vec.
   real(R8P)                                 :: cell_coord(3)           !< Cell-center coordinate vector.
   real(R8P)                                 :: rho                     !< Radial distance in the x-y plane.
   real(R8P)                                 :: F_n                     !< Radial error-function profile.
   real(R8P)                                 :: W_z                     !< Axial window function along z.
   integer(I4P)                              :: b,i,j,k                 !< Counters.

   ! Associate grid, field, and coil data.
   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,         &
            nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                     &
            y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), r_coil=>self%coil%r_coil(n),  &
            l_sol=>self%coil%l_solenoid(n), windings=>self%coil%windings(n),                 &
            nb=>field%nb, x_cell=>field%x_cell, y_cell=>field%y_cell,              &
            z_cell=>field%z_cell, sigma=>self%coil%sigma(n), J_vec=>self%coil%J_vec,         &
            hs=>self%fdv_half_stencil)

   allocate(A(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   A(:,:,:,:,:) = 0.0_R8P

   allocate(J_vec_buffer(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   J_vec_buffer(:,:,:,:,:) = 0.0_R8P

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]

               rho = sqrt((cell_coord(1)-x_c)**2 + (cell_coord(2)-y_c)**2)

               F_n = erf_function(s=rho, mu=r_coil, sigma=sigma)

               W_z = tangential_window(s     = cell_coord(3),     &
                                        smin  = z_c-0.5_R8P*l_sol, &
                                        smax  = z_c+0.5_R8P*l_sol, &
                                        sigma = sigma)

               A(3,i,j,k,b) = -F_n*W_z
            enddo
         enddo
      enddo
   enddo

   call self%compute_curl(hs=hs, ivar=1_I4P, q=A, curl=J_vec_buffer)

   J_vec_buffer = J_vec_buffer * verse

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (maxval(abs(J_vec_buffer(:,i,j,k,b))) < 1e-12_R8P) then
                  J_vec_buffer(:,i,j,k,b) = 0.0_R8P
               endif
            enddo
         enddo
      enddo
   enddo

   J_vec(1:3,:,:,:,:,n) = J_vec(1:3,:,:,:,:,n) + J_vec_buffer

   if (n == 1_I4P) then
      call compute_solenoid_current_density_flux_analytic_z(x_c=x_c, y_c=y_c, z_c=z_c, r_coil=r_coil, &
                                                            l_sol=l_sol, windings=windings,           &
                                                            A=self%coil%A(n), sigma=sigma, n=n,            &
                                                            adjust_amplitude=.true.)
   else
      self%coil%A(n) = self%coil%A(1)
      call compute_solenoid_current_density_flux_analytic_z(x_c=x_c, y_c=y_c, z_c=z_c, r_coil=r_coil, &
                                                            l_sol=l_sol, windings=windings,           &
                                                            A=self%coil%A(n), sigma=sigma, n=n,            &
                                                            adjust_amplitude=.false.)
   endif

   endassociate

   contains

      subroutine compute_solenoid_current_density_flux_analytic_z(x_c, y_c, z_c, r_coil, l_sol, windings, &
                                                                  sigma, A, n, adjust_amplitude)
      !< Analytic amplitude correction for a solenoid with axis parallel to z.
      !< The current is evaluated as the flux of Jy through the section x = x_c-r_coil.
      real(R8P),    intent(in)    :: x_c              !< Solenoid center x-coordinate.
      real(R8P),    intent(in)    :: y_c              !< Solenoid center y-coordinate.
      real(R8P),    intent(in)    :: z_c              !< Solenoid center z-coordinate.
      real(R8P),    intent(in)    :: r_coil           !< Solenoid radius.
      real(R8P),    intent(in)    :: l_sol            !< Solenoid length.
      real(R8P),    intent(in)    :: windings         !< Number of solenoid windings.
      real(R8P),    intent(in)    :: sigma            !< Solenoid Gaussian smoothing width.
      real(R8P),    intent(inout) :: A                !< Solenoid amplitude, corrected in place.
      integer(I4P), intent(in)    :: n                !< Coil number.
      logical,      intent(in)    :: adjust_amplitude !< If true, correct the solenoid amplitude.
      real(R8P), parameter        :: alpha = 3.5_R8P  !< Half-width of the radial integration interval in sigma units.
      real(R8P)                   :: z_l              !< Lower z-boundary of the solenoid axial window.
      real(R8P)                   :: z_u              !< Upper z-boundary of the solenoid axial window.
      real(R8P)                   :: z_1              !< Lower z-integration boundary.
      real(R8P)                   :: z_2              !< Upper z-integration boundary.
      real(R8P)                   :: rho_1            !< Lower radial integration boundary.
      real(R8P)                   :: rho_2            !< Upper radial integration boundary.
      real(R8P)                   :: i_z              !< Integral of W_z over the z-integration interval.
      real(R8P)                   :: de_rho           !< Difference E_R(rho_2)-E_R(rho_1).
      real(R8P)                   :: i_y              !< Analytic flux contribution of Jy.
      real(R8P)                   :: flux_unit        !< Analytic current flux for unit amplitude.
      real(R8P)                   :: flux             !< Analytic current flux before amplitude correction.
      real(R8P)                   :: correction       !< Amplitude correction factor.
      real(R8P)                   :: target_current   !< Target total current, equal to input A*windings.

      z_l = z_c - 0.5_R8P*l_sol
      z_u = z_c + 0.5_R8P*l_sol

      z_1   = z_l - alpha*sigma
      z_2   = z_u + alpha*sigma
      rho_1 = r_coil - alpha*sigma
      rho_2 = r_coil + alpha*sigma

      i_z = 0.5_R8P * (                                        &
            erf_primitive_function(s=z_2, mu=z_l, sigma=sigma) &
          - erf_primitive_function(s=z_1, mu=z_l, sigma=sigma) &
          - erf_primitive_function(s=z_2, mu=z_u, sigma=sigma) &
          + erf_primitive_function(s=z_1, mu=z_u, sigma=sigma))

      de_rho = erf_function(s=rho_2, mu=r_coil, sigma=sigma) &
             - erf_function(s=rho_1, mu=r_coil, sigma=sigma)

      ! For A_z = -E_R(rho) W_z(z), at x = x_c-r_coil, Jy is negative.
      i_y = -i_z*de_rho

      flux_unit = abs(i_y)

      if (flux_unit <= tiny(1.0_R8P)) then
         error stop 'compute_solenoid_current_density_flux_analytic_z: null analytic unit flux'
      endif

      if (adjust_amplitude) then

         target_current = A * windings
         flux           = target_current * flux_unit
         correction     = 1.0_R8P / flux_unit

         print '(A)', mpih%myrankstr//'Current before correction: '//trim(str(flux))
         print '(A)', mpih%myrankstr//trim(str(target_current))//' Target current A('//trim(str(n))//')*N before correction'
         print '(A)', mpih%myrankstr//'Amplitude scaling factor: '//trim(str(correction))

         A = target_current * correction

         print '(A)', mpih%myrankstr//trim(str(A))//' Amplitude A('//trim(str(n))//') after correction'
         print '(A)', mpih%myrankstr//'Final solenoid current '//trim(str(n))//': '//trim(str(A*flux_unit))

      else

         flux = A * flux_unit

         print '(A)', mpih%myrankstr//'Amplitude A('//trim(str(n))//') not corrected: '//trim(str(A))
         print '(A)', mpih%myrankstr//'Final solenoid current '//trim(str(n))//': '//trim(str(flux))

      endif

      endsubroutine compute_solenoid_current_density_flux_analytic_z

   endsubroutine set_solenoid_z

   function erf_function(s, mu, sigma) result(res)
   real(R8P), intent(in) :: s, mu, sigma
   real(R8P)             :: res

   res = erf((s - mu)/(sigma*sqrt(2.0_R8P)))
   endfunction erf_function

   function tangential_window(s, smin, smax, sigma) result(res)
   real(R8P), intent(in) :: s, smin, smax, sigma
   real(R8P)             :: res

   res = 0.5_R8P * (erf_function(s, smin, sigma) - erf_function(s, smax, sigma))
   endfunction tangential_window

   function erf_primitive_function(s, mu, sigma) result(res)
   !< Primitive of erf((s-mu)/(sqrt(2)*sigma)).
   real(R8P), intent(in) :: s
   real(R8P), intent(in) :: mu
   real(R8P), intent(in) :: sigma
   real(R8P)             :: res

   res = (s - mu) * erf_function(s=s, mu=mu, sigma=sigma) &
       + sigma * sqrt(2.0_R8P / acos(-1.0_R8P))           &
       * exp(-((s - mu)**2) / (2.0_R8P * sigma**2))

   endfunction erf_primitive_function

   function build_A_ghost_stencil(hs, dxyz, gc_coord, x_c, y_c, z_c, sigma, &
                                 i_dir_n, i_dir_a, i_dir_b, a1, a2, b1, b2) result(A_gc)
   !< Build ghost cell A array to compute J_vec value through the curl, useful in multiblock application
   integer(I4P), intent(in)  :: hs                                    !< Half-stencil size for the curl computation.
   real(R8P),    intent(in)  :: dxyz(3)                               !< Space steps.
   real(R8P),    intent(in)  :: gc_coord(3)                           !< Ghost cell coordinates.
   real(R8P),    intent(in)  :: x_c, y_c, z_c                         !< Center coordinates.
   real(R8P),    intent(in)  :: sigma                                 !< Width of the Gaussian profile.
   real(R8P),    intent(in)  :: a1, a2                                !< Local coordinates of the rectangular coil along the first tangential direction.
   real(R8P),    intent(in)  :: b1, b2                                !< Local coordinates of the rectangular coil along the second tangential direction.
   integer(I4P), intent(in)  :: i_dir_n, i_dir_a, i_dir_b             !< Direction indices for the normal and tangential directions.
   real(R8P)                 :: A_gc(3, 1:2*hs+1, 1:2*hs+1, 1:2*hs+1) !< Output ghost cell A array for the curl stencil.
   real(R8P),    allocatable :: coordinates_matrix(:,:,:,:)           !< Coordinates matrix for the ghost cell stencil.

   integer(I4P)              :: i,j,k                         !< Counters.

   allocate(coordinates_matrix(1:3, 1:2*hs+1, 1:2*hs+1, 1:2*hs+1))
   do i = 1, 2*hs+1
      do j = 1, 2*hs+1
         do k = 1, 2*hs+1
            coordinates_matrix(:,i,j,k) = gc_coord + [(i-1-hs)*dxyz(1), (j-1-hs)*dxyz(2), (k-1-hs)*dxyz(3)]
         enddo
      enddo
   enddo

   do i = 1, 2*hs+1
      do j = 1, 2*hs+1
         do k = 1, 2*hs+1
            A_gc(:,i,j,k) = rectangular_coil_A_at_point(xp=coordinates_matrix(:,i,j,k), &
                                                        xc=[x_c, y_c, z_c], &
                                                        idir_n=i_dir_n, idir_a=i_dir_a, idir_b=i_dir_b, &
                                                        a1=a1, a2=a2, b1=b1, b2=b2, &
                                                        sigma=sigma)
         enddo
      enddo
   enddo
   endfunction build_A_ghost_stencil

   function rectangular_coil_A_at_point(xp, xc, idir_n, idir_a, idir_b, a1, a2, b1, b2, sigma) result(Ap)
   real(R8P),    intent(in) :: xp(3)   !< Punto di valutazione.
   real(R8P),    intent(in) :: xc(3)   !< Centro spira.
   integer(I4P), intent(in) :: idir_n  !< Componente normale: 1=x, 2=y, 3=z.
   integer(I4P), intent(in) :: idir_a  !< Prima coordinata locale nel piano spira.
   integer(I4P), intent(in) :: idir_b  !< Seconda coordinata locale nel piano spira.
   real(R8P),    intent(in) :: a1, a2  !< Estremi lato locale a.
   real(R8P),    intent(in) :: b1, b2  !< Estremi lato locale b.
   real(R8P),    intent(in) :: sigma   !< Smearing unico della spira.
   real(R8P)                :: Ap(3)

   real(R8P) :: x_a, x_b, x_n
   real(R8P) :: A_1, A_2, A_3, A_4
   real(R8P) :: F_l, W_t, W_n

   Ap(:) = 0.0_R8P

   ! Schema locale nel piano della spira (a,b):
   !
   !          b
   !          ^
   !          |
   !   b2  3 <-------- 2
   !       |            |
   !       |            |
   !       |            |
   !   b1  4 --------> 1
   !          a1       a2  ---> a
   !
   ! Lati:
   !   1: b = b1, tangente lungo +a
   !   2: a = a2, tangente lungo +b
   !   3: b = b2, tangente lungo -a
   !   4: a = a1, tangente lungo -b

   x_a = xp(idir_a)
   x_b = xp(idir_b)
   x_n = xp(idir_n)

   W_n = tangential_window(s=x_n, smin=xc(idir_n)-sigma, &
                                smax=xc(idir_n)+sigma, &
                                sigma=sigma)

   F_l = erf_function(s=x_b, mu=b1, sigma=sigma)
   W_t = tangential_window(s=x_a, smin=a1, smax=a2, sigma=sigma)
   A_1 = F_l*W_t*W_n

   F_l = erf_function(s=x_a, mu=a2, sigma=sigma)
   W_t = tangential_window(s=x_b, smin=b1, smax=b2, sigma=sigma)
   A_2 = -F_l*W_t*W_n

   F_l = erf_function(s=x_b, mu=b2, sigma=sigma)
   W_t = tangential_window(s=x_a, smin=a1, smax=a2, sigma=sigma)
   A_3 = -F_l*W_t*W_n

   F_l = erf_function(s=x_a, mu=a1, sigma=sigma)
   W_t = tangential_window(s=x_b, smin=b1, smax=b2, sigma=sigma)
   A_4 = F_l*W_t*W_n

   Ap(idir_n) = (A_1 + A_2 + A_3 + A_4)

   endfunction rectangular_coil_A_at_point

endmodule adam_prism_common_object
