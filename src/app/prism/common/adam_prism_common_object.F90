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
! PRISM singleton objects
use :: adam_prism_bc_global,       only : bc
use :: adam_prism_ic_global,       only : ic
use :: adam_prism_numerics_global, only : numerics
use :: adam_prism_physics_global,  only : physics
use :: adam_prism_rk_bc_global,    only : rk_bc
use :: adam_prism_time_global,     only : time
! third party modules
use :: motion
use :: penf
use :: stringifor

implicit none
private
public :: prism_common_object

type, extends(equation_object) :: prism_common_object
   !< Maxwell equations system class definition, common data to all backends.
   type(prism_fWLayer_object)            :: fWLayer            !< fWLayer handler.
   type(prism_coil_object)               :: coil               !< Coils handler.
   type(prism_external_fields_object)    :: external_fields    !< External fields handler.
   type(prism_pic_object)                :: pic                !< Particle-in-Cell (PIC) handler.
   type(prism_particle_injection_object) :: particle_injection !< Particle injection handler.
   type(prism_leapfrog_pic_object)       :: leapfrog_pic       !< Leapfrog PIC integrator.
   type(prism_rk_pic_object)             :: rk_pic             !< RK PIC integrator.
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
   contains
      procedure, pass(self) :: allocate_common          !< Allocate common data.
      procedure, pass(self) :: compute_auxiliary_fields !< Compute auxiliary fields.
      procedure, pass(self) :: initialize               !< Initialize the equation common data.
      ! IO methods
      procedure, pass(self) :: load_restart_files      !< Load restart files.
      procedure, pass(self) :: save_energy_error       !< Save energy error history.
      procedure, pass(self) :: save_energy_history     !< Save energy history.
      procedure, pass(self) :: save_divergence_history !< Save divergence history.
      procedure, pass(self) :: save_restart_files      !< Save restart files.
      procedure, pass(self) :: save_xh5f               !< Save simulation data in XH5F format.
      ! coils initialization methods
      procedure, pass(self) :: compute_coil_current_density_flux     !< Compute coil current density fluxes for Maxwell equations.
      procedure, pass(self) :: compute_solenoid_current_density_flux !< Compute solenoid current density fluxes for Maxwell equations.
      procedure, pass(self) :: initialize_coils                      !< Initialize coils.
      procedure, pass(self) :: set_rectangular_coil_x                !< Subroutine to set a rectangular coil source with +-x normal
      procedure, pass(self) :: set_rectangular_coil_y                !< Subroutine to set a rectangular coil source with +-y normal
      procedure, pass(self) :: set_rectangular_coil_z                !< Subroutine to set a rectangular coil source with +-z normal
      procedure, pass(self) :: set_circular_coil_x                   !< Subroutine to set a circular coil source with +-x normal
      procedure, pass(self) :: set_circular_coil_y                   !< Subroutine to set a circular coil source with +-y normal
      procedure, pass(self) :: set_circular_coil_z                   !< Subroutine to set a circular coil source with +-z normal
      procedure, pass(self) :: set_solenoid_x                        !< Subroutine to set a solenoid source with +-x normal
      procedure, pass(self) :: set_solenoid_y                        !< Subroutine to set a solenoid source with +-y normal
      procedure, pass(self) :: set_solenoid_z                        !< Subroutine to set a solenoid source with +-z normal
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

   if (physics%physical_model == PIC_PHYSICAL_MODEL) then
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

   associate(hs=>self%fdv_half_stencil,var_jx=>physics%var_jx)
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
   call numerics%initialize(file_parameters=file_parameters)
   call physics%initialize(file_parameters=file_parameters,                              &
                                reconstruction_vars=numerics%reconstruction_vars,        &
                                div_corr_var=numerics%div_corr_var,                      &
                                constrained_transport_D=numerics%constrained_transport_D,&
                                constrained_transport_B=numerics%constrained_transport_B)
   self%nv_c   => physics%nv_c
   self%nv_s   => physics%nv_s
   self%nv_cl  => physics%nv_cl
   !self%nv_pic => physics%nv_pic
   call self%equation_object%initialize(filename=filename, memory_avail=memory_avail, nv=physics%nv, verbose=verbose_)
   call bc%initialize(file_parameters=file_parameters)
   call grid%set_bc_type(bc_type=bc%bc_type)
   if (physics%physical_model == PIC_PHYSICAL_MODEL) &
      call self%pic%initialize(file_parameters=file_parameters)
   if (physics%physical_model == PIC_PHYSICAL_MODEL) &
      call self%particle_injection%initialize(file_parameters=file_parameters, pic=self%pic)
   call time%initialize(file_parameters=file_parameters)
   call ic%initialize(file_parameters=file_parameters)
   call self%fWLayer%initialize(file_parameters=file_parameters, physics=physics)
   call self%coil%initialize(file_parameters=file_parameters)
   call self%external_fields%initialize(file_parameters=file_parameters)
   if (numerics%scheme_time==NUM_SCHEME_TIME_RUNGE_KUTTA) &
      call rk_bc%initialize(file_parameters=file_parameters, rk=rk, physics=physics)
   if (physics%physical_model == PIC_PHYSICAL_MODEL) then
      if (self%pic%scheme_time==NUM_SCHEME_TIME_PIC_LEAPFROG) &
         call self%leapfrog_pic%initialize(file_parameters=file_parameters, pic=self%pic)
      if (self%pic%scheme_time==NUM_SCHEME_TIME_PIC_RUNGE_KUTTA) &
         call self%rk_pic%initialize(file_parameters=file_parameters, rk=rk, pic=self%pic)
   endif
   call check_ngc_number
   call self%allocate_common
   call io_initialize
   endassociate
   if (verbose_) call mpih%print_message('prism_common_object%initialize finish')
   contains
      subroutine check_ngc_number
      !< Check if the number of ghost cells is consistent with the numerical schemes used, if not an error is echoed and
      !< the simulation is stop.

      if (numerics%scheme_space==NUM_SCHEME_SPACE_WENO) then
         if (weno%S > grid%ngc) &
            call mpih%error_stop(msg=': ghost cells number (ngc) must be >= of weno stencil number (weno%S):'//&
                                      ' ngc='//trim(str(grid%ngc))//' weno%S='//trim(str(weno%S)))
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

      add_rho = physics%physical_model == PIC_PHYSICAL_MODEL
      associate(add_phi=>numerics%constrained_transport_D, &
                add_psi=>numerics%constrained_transport_B)
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

                   q_name = ['div_D ','div_B ','div_J ','div_04','div_05','div_06']
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

   call self%adam%load_restart_files(basename=self%io%restart_basename, t=t, time=time, q=self%q)
   call self%adam%make_comm_local_maps_ghost_bc
   endsubroutine load_restart_files

   subroutine save_divergence_history(self, is_to_open, is_to_close)
   !< Save divergence history.
   class(prism_common_object), intent(inout)        :: self        !< The equation.
   logical,                    intent(in), optional :: is_to_open  !< Flag to open  file before first saving.
   logical,                    intent(in), optional :: is_to_close !< Flag to close file after last saving.
   real(R8P)                                        :: max_div_D   !< Maximum of divergence of D field.
   real(R8P)                                        :: max_div_B   !< Maximum of divergence of B
   real(R8P)                                        :: max_div_J   !< Maximum of divergence of J field.
   real(R8P)                                        :: r           !< Auxiliary variable to identify fWL presence

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, C=>self%fWLayer%C, hs=>self%fdv_half_stencils(1))
   r = nint(real(C)/(real(C)+1_I4P))
   if (time%is_to_save(it_save=self%io%divergence_history_save)) then
      max_div_D = maxval(abs(self%divergence(1,1+r*(C+hs):ni-r*(C+hs-1_I4P),1+r*(C+hs):nj-r*(C+hs-1_I4P), &
                                             1+r*(C+hs):nk-r*(C+hs-1_I4P),:)))
      max_div_B = maxval(abs(self%divergence(2,1+r*(C+hs):ni-r*(C+hs-1_I4P),1+r*(C+hs):nj-r*(C+hs-1_I4P), &
                                             1+r*(C+hs):nk-r*(C+hs-1_I4P),:)))
      max_div_J = maxval(abs(self%divergence(3,1+r*(C+hs):ni-r*(C+hs-1_I4P),1+r*(C+hs):nj-r*(C+hs-1_I4P), &
                                             1+r*(C+hs):nk-r*(C+hs-1_I4P),:)))
      call self%io%save_divergence_history(it=time%it,time=time%time,blocks_number=self%blocks_number, &
                                           div_D=max_div_D,div_B=max_div_B,div_J=max_div_J, &
                                           is_to_open=is_to_open,is_to_close=is_to_close)
   endif
   endassociate
   endsubroutine save_divergence_history

   subroutine save_energy_error(self, is_to_open, is_to_close)
   !< Save energy error history.
   class(prism_common_object), intent(inout)        :: self        !< The equation.
   logical,                    intent(in), optional :: is_to_open  !< Flag to open  file before first saving.
   logical,                    intent(in), optional :: is_to_close !< Flag to close file after last saving.

   if (time%is_to_save(it_save=self%io%energy_error_save)) then
      call self%io%save_energy_error(it=time%it,time=time%time,blocks_number=self%blocks_number,                  &
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

   if (time%is_to_save(it_save=self%io%energy_history_save)) then
      call self%io%save_energy_history(it=time%it,time=time%time,blocks_number=self%blocks_number, &
                                       energy_D=self%energy_D,energy_B=self%energy_B,                        &
                                       coil_power=self%coil_power,Poynting_flux=self%Poynting_flux,          &
                                       is_to_open=is_to_open,is_to_close=is_to_close)
   endif
   endsubroutine save_energy_history

   subroutine save_restart_files(self)
   !< Save restart files.
   class(prism_common_object), intent(inout) :: self !< The equation.

   call mpih%barrier(tictoc=.true.)
   call mpih%print_message('save restart files t: '//trim(str(time%it,.true.))//', time: '//&
                                    trim(str(time%time,.true.)))
   call self%adam%save_restart_files(basename=self%io%restart_basename, t=time%it, time=time%time, q=self%q)
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
   call mpih%print_message('save HDF5 files t: '//trim(str(time%it,.true.))//', time: '//&
                                trim(str(time%time,.true.)))

   output_basename_ = trim(self%io%output_basename)//'-'//trim(strz(time%it,9))
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
      call self%open_block_xh5f(xh5f=xh5f, b=b, nijk=nijk, t=time%it, time=time%time)

      call self%io%save_field(xh5f=xh5f, block_name=bn, ijk=ijk, nijk=nijk, q=self%q(:,:,:,:,b), q_name=self%q_name)

      if (self%coil%total_coils_number>0) then
         call self%io%save_field(xh5f=xh5f, block_name=bn, ijk=ijk, nijk=nijk, &
                                 q=self%coil%coil_flag(:,:,:,b), q_name=self%coil%coil_flag_name)
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
   subroutine compute_coil_current_density_flux(self, n, adjust_amplitude)
   !< Subroutine to adjust current amplitude in order to match the input one
   class(prism_common_object), intent(inout) :: self             !< Cpu object.
   integer(I4P),               intent(in)    :: n                !< Coil number.
   logical,                    intent(in)    :: adjust_amplitude !< If true, adjust amplitude
   real(R8P)                                 :: x_s, y_s, z_s    !< Flux center coordinates
   real(R8P)                                 :: flux, correction !< Computed flux
   real(R8P)                                 :: l                !< Auxiliary variable to compute flux at a certain distance
   integer(I4P)                              :: i_s, j_s, k_s    !< Flux center cell coordinates
   integer(I4P)                              :: i, j, k          !< Counter.

   associate(blocks_number=>self%blocks_number, ni=>self%ni, nj=>self%nj,                 &
             nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                      &
             lx=>self%coil%lx(n), ly=>self%coil%ly(n), coil_flag =>self%coil%coil_flag,   &
             r_coil=>self%coil%r_coil(n), coil_type=>self%coil%coil_type(n),              &
             y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), dx=>field%dxyz(1,:), &
             dy=>field%dxyz(2,:), dz=>field%dxyz(3,:), normal=>self%coil%normal(n),       &
             nb=>field%nb, x_cell=>field%x_cell, y_cell=>field%y_cell,                    &
             z_cell=>field%z_cell, sigma=>self%coil%sigma(n), J_vec=>self%coil%J_vec,     &
             e_min => grid%domain_emin, e_max => grid%domain_emax,                        &
             q=>self%q)

   !Per ora la imposto per griglia uniforme monoblocco. Vedremo come estendere il problema
   flux = 0.0_R8P
   correction = 1.0_R8P

   if (coil_type == COIL_TYPE_RECTANGULAR) then
      l = lx
   elseif (coil_type == COIL_TYPE_CIRCULAR) then
      l = r_coil*2._R8P
   endif

   if (normal == NORMAL_P_X .or. normal == NORMAL_M_X) then !Valuto corrente lungo z

      x_s = x_c
      y_s = y_c - l/2
      z_s = z_c
      i_s = ceiling((x_s - e_min(1)) / dx(1))
      j_s = ceiling((y_s - e_min(2)) / dy(1))
      k_s = ceiling((z_s - e_min(3)) / dz(1))
      !i_s = floor( (x_s - e_min(1)) / dx(1) ) + 1_I4P
      !j_s = floor( (y_s - e_min(2)) / dy(1) ) + 1_I4P
      !k_s = floor( (z_s - e_min(3)) / dz(1) ) + 1_I4P

      ! Clamp su celle fisiche
      i_s = max(1_I4P, min(ni, i_s))
      j_s = max(1_I4P, min(nj, j_s))
      k_s = max(1_I4P, min(nk, k_s))

      do j = j_s - nint(3.5_R8P*self%coil%sigma(n)), j_s + nint(3.5_R8P*self%coil%sigma(n))
         do i = i_s - nint(3.5_R8P*self%coil%sigma(n)), i_s + nint(3.5_R8P*self%coil%sigma(n))
               flux = flux + J_vec(3,i,j,k_s,1,n)*dx(1)*dy(1)
         enddo
      enddo

   elseif (normal == NORMAL_P_Y .or. normal == NORMAL_M_Y) then !Valuto corrente lungo z

      x_s = x_c - l/2
      y_s = y_c
      z_s = z_c
      i_s = ceiling((x_s - e_min(1)) / dx(1))
      j_s = ceiling((y_s - e_min(2)) / dy(1))
      k_s = ceiling((z_s - e_min(3)) / dz(1))
      !i_s = floor( (x_s - e_min(1)) / dx(1) ) + 1_I4P
      !j_s = floor( (y_s - e_min(2)) / dy(1) ) + 1_I4P
      !k_s = floor( (z_s - e_min(3)) / dz(1) ) + 1_I4P

      ! Clamp su celle fisiche
      i_s = max(1_I4P, min(ni, i_s))
      j_s = max(1_I4P, min(nj, j_s))
      k_s = max(1_I4P, min(nk, k_s))

      do j = j_s - nint(3.5_R8P*self%coil%sigma(n)), j_s + nint(3.5_R8P*self%coil%sigma(n))
         do i = i_s - nint(3.5_R8P*self%coil%sigma(n)), i_s + nint(3.5_R8P*self%coil%sigma(n))
               flux = flux + J_vec(3,i,j,k_s,1,n)*dx(1)*dy(1)
         enddo
      enddo

   elseif (normal == NORMAL_P_Z .or. normal == NORMAL_M_Z) then !Valuto corrente lungo y

      x_s = x_c - l/2
      y_s = y_c
      z_s = z_c
      i_s = ceiling((x_s - e_min(1)) / dx(1))
      j_s = ceiling((y_s - e_min(2)) / dy(1))
      k_s = ceiling((z_s - e_min(3)) / dz(1))
      !i_s = floor( (x_s - e_min(1)) / dx(1) ) + 1_I4P
      !j_s = floor( (y_s - e_min(2)) / dy(1) ) + 1_I4P
      !k_s = floor( (z_s - e_min(3)) / dz(1) ) + 1_I4P

      ! Clamp su celle fisiche
      i_s = max(1_I4P, min(ni, i_s))
      j_s = max(1_I4P, min(nj, j_s))
      k_s = max(1_I4P, min(nk, k_s))

      do k = k_s - nint(3.5_R8P*self%coil%sigma(n)), k_s + nint(3.5_R8P*self%coil%sigma(n))
         do i = i_s - nint(3.5_R8P*self%coil%sigma(n)), i_s + nint(3.5_R8P*self%coil%sigma(n))
               flux = flux + J_vec(2,i,j_s,k,1,n)*dx(1)*dz(1)
         enddo
      enddo

   endif
   flux = abs(flux)*self%coil%A(n)

   if (adjust_amplitude) then
      print '(A)', mpih%myrankstr//'Valore corrente pre correzione: '//trim(str(flux))
      print '(A)', mpih%myrankstr//trim(str(self%coil%A(n)))//'Ampiezza A('//trim(str(n))//') pre correzione'
      correction = (self%coil%A(n)/flux)
      print '(A)', mpih%myrankstr//'Scaling factor ampiezza: '//trim(str(correction))
      self%coil%A(n) = self%coil%A(n)*correction
      print '(A)', mpih%myrankstr//trim(str(self%coil%A(n)))//'Ampiezza A('//trim(str(n))//') post correzione'
      print '(A)', mpih%myrankstr//'Valore finale corrente spira'//trim(str(n))//': '//trim(str(flux*correction))
   else
      print '(A)', mpih%myrankstr//'Ampiezza A('//trim(str(n))//') non corretta: '//trim(str(self%coil%A(n)))
      print '(A)', mpih%myrankstr//'Valore finale corrente spira'//trim(str(n))//': '//trim(str(flux))
   endif

   endassociate
   endsubroutine compute_coil_current_density_flux

   subroutine compute_solenoid_current_density_flux(self, n, adjust_amplitude)
   !< Subroutine to adjust current amplitude in order to match the input one
   class(prism_common_object), intent(inout) :: self             !< Cpu object.
   integer(I4P),               intent(in)    :: n                !< Coil number.
   logical,                    intent(in)    :: adjust_amplitude !< If true, adjust amplitude
   real(R8P)                                 :: x_s, y_s, z_s    !< Flux center coordinates
   real(R8P)                                 :: flux, correction !< Computed flux
   integer(I4P)                              :: delta            !< Auxiliary variable to compute flux
   integer(I4P)                              :: i_s, j_s, k_s    !< Flux center cell coordinates
   integer(I4P)                              :: i, j, k          !< Counter.

   associate(blocks_number=>self%blocks_number, ni=>self%ni, nj=>self%nj,                 &
             nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                      &
             lx=>self%coil%lx(n), ly=>self%coil%ly(n), coil_flag =>self%coil%coil_flag,   &
             l_sol=>self%coil%l_solenoid(n), windings=>self%coil%windings(n),             &
             r_coil=>self%coil%r_coil(n), coil_type=>self%coil%coil_type(n),              &
             y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), dx=>field%dxyz(1,:), &
             dy=>field%dxyz(2,:), dz=>field%dxyz(3,:), normal=>self%coil%normal(n),       &
             nb=>field%nb, x_cell=>field%x_cell, y_cell=>field%y_cell,                    &
             z_cell=>field%z_cell, sigma=>self%coil%sigma(n), J_vec=>self%coil%J_vec,     &
             e_min => grid%domain_emin, e_max => grid%domain_emax,                        &
             q=>self%q)

   !Per ora la imposto per griglia uniforme monoblocco. Vedremo come estendere il problema
   flux = 0.0_R8P
   correction = 1.0_R8P

   if (normal == NORMAL_P_X .or. normal == NORMAL_M_X) then !Valuto corrente lungo z

      x_s = x_c
      y_s = y_c - r_coil
      z_s = z_c
      i_s = ceiling((x_s - e_min(1)) / dx(1))
      j_s = ceiling((y_s - e_min(2)) / dy(1))
      k_s = ceiling((z_s - e_min(3)) / dz(1))
      !i_s = floor( (x_s - e_min(1)) / dx(1) ) + 1_I4P
      !j_s = floor( (y_s - e_min(2)) / dy(1) ) + 1_I4P
      !k_s = floor( (z_s - e_min(3)) / dz(1) ) + 1_I4P
      delta = ceiling((l_sol/2)/dx(1))
      ! Clamp su celle fisiche
      i_s = max(1_I4P, min(ni, i_s))
      j_s = max(1_I4P, min(nj, j_s))
      k_s = max(1_I4P, min(nk, k_s))

      do j = j_s - nint(3.5_R8P*sigma), j_s + nint(3.5_R8P*sigma)
         do i = (i_s - delta) - nint(3.5_R8P*sigma), (i_s + delta) + nint(3.5_R8P*sigma)
               flux = flux + J_vec(3,i,j,k_s,1,n)*dx(1)*dy(1)
         enddo
      enddo

   elseif (normal == NORMAL_P_Y .or. normal == NORMAL_M_Y) then !Valuto corrente lungo z

      x_s = x_c - r_coil
      y_s = y_c
      z_s = z_c
      i_s = ceiling((x_s - e_min(1)) / dx(1))
      j_s = ceiling((y_s - e_min(2)) / dy(1))
      k_s = ceiling((z_s - e_min(3)) / dz(1))
      !i_s = floor( (x_s - e_min(1)) / dx(1) ) + 1_I4P
      !j_s = floor( (y_s - e_min(2)) / dy(1) ) + 1_I4P
      !k_s = floor( (z_s - e_min(3)) / dz(1) ) + 1_I4P
      delta = ceiling((l_sol/2)/dy(1))

      ! Clamp su celle fisiche
      i_s = max(1_I4P, min(ni, i_s))
      j_s = max(1_I4P, min(nj, j_s))
      k_s = max(1_I4P, min(nk, k_s))

      do j = (j_s - delta) - nint(3.5_R8P*sigma), (j_s + delta) + nint(3.5_R8P*sigma)
         do i = i_s - nint(3.5_R8P*sigma), i_s + nint(3.5_R8P*sigma)
               flux = flux + J_vec(3,i,j,k_s,1,n)*dx(1)*dy(1)
         enddo
      enddo

   elseif (normal == NORMAL_P_Z .or. normal == NORMAL_M_Z) then !Valuto corrente lungo y

      x_s = x_c - r_coil
      y_s = y_c
      z_s = z_c
      i_s = ceiling((x_s - e_min(1)) / dx(1))
      j_s = ceiling((y_s - e_min(2)) / dy(1))
      k_s = ceiling((z_s - e_min(3)) / dz(1))
      !i_s = floor( (x_s - e_min(1)) / dx(1) ) + 1_I4P
      !j_s = floor( (y_s - e_min(2)) / dy(1) ) + 1_I4P
      !k_s = floor( (z_s - e_min(3)) / dz(1) ) + 1_I4P
      delta = ceiling((l_sol/2)/dz(1))

      ! Clamp su celle fisiche
      i_s = max(1_I4P, min(ni, i_s))
      j_s = max(1_I4P, min(nj, j_s))
      k_s = max(1_I4P, min(nk, k_s))

      do k = (k_s - delta) - nint(3.5_R8P*sigma), (k_s + delta) + nint(3.5_R8P*sigma)
         do i = i_s - nint(3.5_R8P*sigma), i_s + nint(3.5_R8P*sigma)
               flux = flux + J_vec(2,i,j_s,k,1,n)*dx(1)*dz(1)
         enddo
      enddo

   endif
   flux = abs(flux)*self%coil%A(n)*windings !Valore calcolato (corrente complessiva, deve essere N*A)
                                            !A*windings è il valore target
   if (adjust_amplitude) then
      print '(A)', mpih%myrankstr//'Valore corrente pre correzione: '//trim(str(flux))
      print '(A)', mpih%myrankstr//trim(str(self%coil%A(n)*windings))//'Ampiezza A('//trim(str(n))//')*N pre correzione'
      correction = (self%coil%A(n)*windings/flux)
      print '(A)', mpih%myrankstr//'Scaling factor ampiezza: '//trim(str(correction))
      self%coil%A(n) = self%coil%A(n)*correction*windings
      print '(A)', mpih%myrankstr//trim(str(self%coil%A(n)))//'Ampiezza A('//trim(str(n))//')*N post correzione'
      print '(A)', mpih%myrankstr//'Valore finale corrente solenoide'//trim(str(n))//': '//trim(str(flux*correction))
   else
      print '(A)', mpih%myrankstr//'Ampiezza A('//trim(str(n))//') non corretta: '//trim(str(self%coil%A(n)))
      print '(A)', mpih%myrankstr//'Valore finale corrente solenoide'//trim(str(n))//': '//trim(str(flux))
   endif

   endassociate
   endsubroutine compute_solenoid_current_density_flux

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
      call self%compute_divergence(hs=self%fdv_half_stencils(1),ivar=1_I4P,q=self%coil%J_vec(1:3,:,:,:,:,n),&
                                   divergence=self%divergence(3,:,:,:,:))
      print '(A)', mpih%myrankstr//'Divergenza J vec della spira: ' &
                  //trim(str(n))//' pari a: '//trim(str(maxval(abs(self%divergence(3,:,:,:,:)))))
   enddo
   endsubroutine initialize_coils

   subroutine set_rectangular_coil_x(self, n, verse)
   class(prism_common_object),  intent(inout) :: self                    !< Cpu object.
   integer(I4P),                intent(in)    :: n                       !< Coil number.
   real(R8P),                   intent(in)    :: verse                   !< Coil normal direction, +1=+x, -1=-x.
   real(R8P),                   allocatable   :: A(:,:,:,:,:)            !< Campo vettoriale totale della spira
   real(R8P),                   allocatable   :: J_vec_buffer(:,:,:,:,:) !< Variabile buffer per coil%J_vec
   real(R8P)                                  :: A_1, A_2, A_3, A_4      !< Campo vettoriale lati spira
   real(R8P)                                  :: cell_coord(3)           !< Vettore posizione centro cella
   real(R8P)                                  :: y_d, y_t, z_b, z_f
   real(R8P)                                  :: F_n, W_t, W_x
   integer(I4P)                               :: b,i,j,k                 !< Counter.

   !associo per dati su posizioni delle celle e contatori
   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,               &
            nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                      &
            lx=>self%coil%lx(n), ly=>self%coil%ly(n), coil_flag =>self%coil%coil_flag,   &
            y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), dx=>field%dxyz(1,:), &
            dy=>field%dxyz(2,:), dz=>field%dxyz(3,:), normal=>self%coil%normal(n),       &
            nb=>field%nb, x_cell=>field%x_cell, y_cell=>field%y_cell,                    &
            z_cell=>field%z_cell, sigma=>self%coil%sigma(n), J_vec=>self%coil%J_vec,     &
            hs=>self%fdv_half_stencil)

   !Fisso estremi della spira rettangolare con normale parallela a x, e centro in (x_c, y_c, z_c)
   y_d = -lx/2 + y_c
   y_t = +lx/2 + y_c
   z_b = -ly/2 + z_c
   z_f = +ly/2 + z_c

   allocate(A(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   A(:,:,:,:,:) = 0.0_R8P
   allocate(J_vec_buffer(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   J_vec_buffer(:,:,:,:,:) = 0.0_R8P

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
               !Primo filo (y = y_d, tangenziale lungo z)
               F_n = erf_function(s=cell_coord(2), mu = y_d, sigma = sigma*dy(b))
               W_t = tangential_window(s=cell_coord(3), smin = z_b, smax = z_f, sigma = sigma*dz(b))
               W_x = tangential_window(s=cell_coord(1), smin = x_c-sigma*dx(b), smax = x_c+sigma*dx(b), sigma=sigma*dx(b))
               A_1 = F_n*W_t*W_x
               !Secondo filo (z = z_f, tangenziale lungo y)
               F_n = erf_function(s=cell_coord(3), mu = z_f, sigma = sigma*dz(b))
               W_t = tangential_window(s=cell_coord(2), smin = y_d, smax = y_t, sigma = sigma*dy(b))
               W_x = tangential_window(s=cell_coord(1), smin = x_c-sigma*dx(b), smax = x_c+sigma*dx(b), sigma=sigma*dx(b))
               A_2 = -F_n*W_t*W_x
               !Terzo filo (y = y_t, tangenziale lungo z)
               F_n = erf_function(s=cell_coord(2), mu = y_t, sigma = sigma*dy(b))
               W_t = tangential_window(s=cell_coord(3), smin = z_b, smax = z_f, sigma = sigma*dz(b))
               W_x = tangential_window(s=cell_coord(1), smin = x_c-sigma*dx(b), smax = x_c+sigma*dx(b), sigma=sigma*dx(b))
               A_3 = -F_n*W_t*W_x
               !Quarto filo (z = z_b, tangenziale lungo y)
               F_n = erf_function(s=cell_coord(3), mu = z_b, sigma = sigma*dz(b))
               W_t = tangential_window(s=cell_coord(2), smin = y_d, smax = y_t, sigma = sigma*dy(b))
               W_x = tangential_window(s=cell_coord(1), smin = x_c-sigma*dx(b), smax = x_c+sigma*dx(b), sigma=sigma*dx(b))
               A_4 = F_n*W_t*W_x
               !Somma dei campi vettoriali
               A(1,i,j,k,b) = A_1 + A_2 + A_3 + A_4
            enddo
         enddo
      enddo
   enddo

   call self%compute_curl(hs=hs, ivar=1_I4P,q=A,curl=J_vec_buffer)
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

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (J_vec_buffer(1,i,j,k,b) /= 0.0_R8P .or. J_vec_buffer(2,i,j,k,b) /= 0.0_R8P &
                  .or. J_vec_buffer(3,i,j,k,b) /= 0.0_R8P) then
                  coil_flag(i,j,k,b) = n
               endif
            enddo
         enddo
      enddo
   enddo

   J_vec(1:3,:,:,:,:,n) = J_vec(1:3,:,:,:,:,n) + J_vec_buffer
   endassociate

   if (n == 1_I4P) then
      call self%compute_coil_current_density_flux(n=n, adjust_amplitude=.true.)
   else
      self%coil%A(n) = self%coil%A(1)
      call self%compute_coil_current_density_flux(n=n, adjust_amplitude=.false.)
   endif
   endsubroutine set_rectangular_coil_x

   subroutine set_rectangular_coil_y(self, n, verse)
   class(prism_common_object), intent(inout) :: self                    !< Cpu object.
   integer(I4P),               intent(in)    :: n                       !< Coil number.
   real(R8P),                  intent(in)    :: verse                   !< Coil normal direction, +1=+y, -1=-y.
   real(R8P),                  allocatable   :: A(:,:,:,:,:)            !< Campo vettoriale totale della spira
   real(R8P),                  allocatable   :: J_vec_buffer(:,:,:,:,:) !< Variabile buffer per coil%J_vec
   real(R8P)                                 :: A_1, A_2, A_3, A_4      !< Campo vettoriale lati spira
   real(R8P)                                 :: cell_coord(3)           !< Vettore posizione centro cella
   real(R8P)                                 :: x_l, x_r, z_b, z_f
   real(R8P)                                 :: F_n, W_t, W_y
   integer(I4P)                              :: b,i,j,k                 !< Counter.

   !associo per dati su posizioni delle celle e contatori
   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,               &
            nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                      &
            lx=>self%coil%lx(n), ly=>self%coil%ly(n), coil_flag =>self%coil%coil_flag,   &
            y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), dx=>field%dxyz(1,:), &
            dy=>field%dxyz(2,:), dz=>field%dxyz(3,:), normal=>self%coil%normal(n),       &
            nb=>field%nb, x_cell=>field%x_cell, y_cell=>field%y_cell,                    &
            z_cell=>field%z_cell, sigma=>self%coil%sigma(n), J_vec=>self%coil%J_vec,     &
            hs=>self%fdv_half_stencil)

   !Fisso estremi della spira rettangolare con normale parallela a y, e centro in (x_c, y_c, z_c)
   x_l = -lx/2 + x_c
   x_r = +lx/2 + x_c
   z_b = -ly/2 + z_c
   z_f = +ly/2 + z_c

   allocate(A(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   A(:,:,:,:,:) = 0.0_R8P
   allocate(J_vec_buffer(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   J_vec_buffer(:,:,:,:,:) = 0.0_R8P

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
               !Primo filo (z = z_b, tangenziale lungo x)
               F_n = erf_function(s=cell_coord(3), mu = z_b, sigma = sigma*dz(b))
               W_t = tangential_window(s=cell_coord(1), smin = x_l, smax = x_r, sigma = sigma*dx(b))
               W_y = tangential_window(s=cell_coord(2), smin = y_c-sigma*dy(b), smax = y_c+sigma*dy(b), sigma=sigma*dy(b))
               A_1 = F_n*W_t*W_y
               !Secondo filo (x = x_r, tangenziale lungo z)
               F_n = erf_function(s=cell_coord(1), mu = x_r, sigma = sigma*dx(b))
               W_t = tangential_window(s=cell_coord(3), smin = z_b, smax = z_f, sigma = sigma*dz(b))
               W_y = tangential_window(s=cell_coord(2), smin = y_c-sigma*dy(b), smax = y_c+sigma*dy(b), sigma=sigma*dy(b))
               A_2 = -F_n*W_t*W_y
               !Terzo filo (z = z_f, tangenziale lungo x)
               F_n = erf_function(s=cell_coord(3), mu = z_f, sigma = sigma*dz(b))
               W_t = tangential_window(s=cell_coord(1), smin = x_l, smax = x_r, sigma = sigma*dx(b))
               W_y = tangential_window(s=cell_coord(2), smin = y_c-sigma*dy(b), smax = y_c+sigma*dy(b), sigma=sigma*dy(b))
               A_3 = -F_n*W_t*W_y
               !Quarto filo (x = x_l, tangenziale lungo z)
               F_n = erf_function(s=cell_coord(1), mu = x_l, sigma = sigma*dx(b))
               W_t = tangential_window(s=cell_coord(3), smin = z_b, smax = z_f, sigma = sigma*dz(b))
               W_y = tangential_window(s=cell_coord(2), smin = y_c-sigma*dy(b), smax = y_c+sigma*dy(b), sigma=sigma*dy(b))
               A_4 = F_n*W_t*W_y
               !Somma dei campi vettoriali
               A(2,i,j,k,b) = A_1 + A_2 + A_3 + A_4
            enddo
         enddo
      enddo
   enddo

   call self%compute_curl(hs=hs,ivar=1_I4P,q=A,curl=J_vec_buffer)
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

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (J_vec_buffer(1,i,j,k,b) /= 0.0_R8P .or. J_vec_buffer(2,i,j,k,b) /= 0.0_R8P &
                  .or. J_vec_buffer(3,i,j,k,b) /= 0.0_R8P) then
                  coil_flag(i,j,k,b) = n
               endif
            enddo
         enddo
      enddo
   enddo

   J_vec(1:3,:,:,:,:,n) = J_vec(1:3,:,:,:,:,n) + J_vec_buffer
   endassociate
   if (n == 1_I4P) then
      call self%compute_coil_current_density_flux(n=n, adjust_amplitude=.true.)
   else
      self%coil%A(n) = self%coil%A(1)
      call self%compute_coil_current_density_flux(n=n, adjust_amplitude=.false.)
   endif
   endsubroutine set_rectangular_coil_y

   subroutine set_rectangular_coil_z(self, n, verse)
   class(prism_common_object),      intent(inout) :: self                    !< Cpu object.
   integer(I4P),                    intent(in)    :: n                       !< Coil number.
   real(R8P),                       intent(in)    :: verse                   !< Coil normal direction, +1=+z, -1=-z.
   real(R8P),                       allocatable   :: A(:,:,:,:,:)            !< Campo vettoriale totale della spira, somma dei campi A_1, A_2, A_3 e A_4
   real(R8P),                       allocatable   :: J_vec_buffer(:,:,:,:,:) !< Variabile buffer per il campo di corrente da assegnare alla variabile coil%J_vec
   real(R8P)                                      :: A_1, A_2, A_3, A_4      !< Campo vettoriale lati spira
   real(R8P)                                      :: c_c(3)                  !< Vettore posizione centro spira
   real(R8P)                                      :: cell_coord(3)           !< Vettore posizione centro cella
   real(R8P)                                      :: x_l, x_r, y_d, y_t
   real(R8P)                                      :: F_n, W_t, W_z
   integer(I4P)                                   :: b,i,j,k                 !< Counter.

   !associo per dati su posizioni delle celle e contatori
   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,               &
            nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                      &
            lx=>self%coil%lx(n), ly=>self%coil%ly(n), coil_flag =>self%coil%coil_flag,   &
            y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), dx=>field%dxyz(1,:), &
            dy=>field%dxyz(2,:), dz=>field%dxyz(3,:), normal=>self%coil%normal(n),       &
            nb=>field%nb, x_cell=>field%x_cell, y_cell=>field%y_cell,                    &
            z_cell=>field%z_cell, sigma=>self%coil%sigma(n), J_vec=>self%coil%J_vec,     &
            hs=>self%fdv_half_stencil)

   !Fisso estremi della spira rettangolare con normale parallela a z, e centro in (x_c, y_c, z_c)
   x_l = -lx/2 +x_c
   x_r = +lx/2 +x_c
   y_d = -ly/2 +y_c
   y_t = +ly/2 +y_c

   allocate(A(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   A(:,:,:,:,:) = 0.0_R8P
   allocate(J_vec_buffer(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   J_vec_buffer(:,:,:,:,:) = 0.0_R8P

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
               !Primo filo
               F_n = erf_function(s=cell_coord(2), mu = y_d, sigma = sigma*dy(b))
               W_t = tangential_window(s=cell_coord(1), smin = x_l, smax = x_r, sigma = sigma*dx(b))
               W_z = tangential_window(s=cell_coord(3), smin =z_c-sigma*dz(b), smax=z_c+sigma*dz(b), sigma=sigma*dz(b))
               A_1 = F_n*W_t*W_z
               !Secondo filo
               F_n = erf_function(s=cell_coord(1), mu = x_r, sigma = sigma*dx(b))
               W_t = tangential_window(s=cell_coord(2), smin = y_d, smax = y_t, sigma = sigma*dy(b))
               W_z = tangential_window(s=cell_coord(3), smin =z_c-sigma*dz(b), smax=z_c+sigma*dz(b), sigma=sigma*dz(b))
               A_2 = -F_n*W_t*W_z
               !Terzo filo
               F_n = erf_function(s=cell_coord(2), mu = y_t, sigma = sigma*dy(b))
               W_t = tangential_window(s=cell_coord(1), smin = x_l, smax = x_r, sigma = sigma*dx(b))
               W_z = tangential_window(s=cell_coord(3), smin =z_c-sigma*dz(b), smax=z_c+sigma*dz(b), sigma=sigma*dz(b))
               A_3 = -F_n*W_t*W_z
               !Quarto filo
               F_n = erf_function(s=cell_coord(1), mu = x_l, sigma = sigma*dx(b))
               W_t = tangential_window(s=cell_coord(2), smin = y_d, smax = y_t, sigma = sigma*dy(b))
               W_z = tangential_window(s=cell_coord(3), smin =z_c-sigma*dz(b), smax=z_c+sigma*dz(b), sigma=sigma*dz(b))
               A_4 = F_n*W_t*W_z
               !Somma dei campi vettoriali
               A(3,i,j,k,b) = A_1 + A_2+ A_3 + A_4
            enddo
         enddo
      enddo
   enddo
   call self%compute_curl(hs=hs, ivar=1_I4P,q=A,curl=J_vec_buffer)
   J_vec_buffer = J_vec_buffer * verse !* 10.0_R8P**4._R8P

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

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (J_vec_buffer(1,i,j,k,b) /= 0.0_R8P .or. J_vec_buffer(2,i,j,k,b) /= 0.0_R8P &
                  .or. J_vec_buffer(3,i,j,k,b) /= 0.0_R8P) then
                  coil_flag(i,j,k,b) = n
               endif
            enddo
         enddo
      enddo
   enddo
   J_vec(1:3,:,:,:,:,n) = J_vec(1:3,:,:,:,:,n) + J_vec_buffer
   endassociate

   if (n == 1_I4P) then
      call self%compute_coil_current_density_flux(n=n, adjust_amplitude=.true.)
   else
      self%coil%A(n) = self%coil%A(1)
      call self%compute_coil_current_density_flux(n=n, adjust_amplitude=.false.)
   endif

   !!Calcolo divergenza di J (ampiezza massima) prima della correzione alla Poisson
   !self%coil%J_vec(n,1:3,:,:,:,:) = self%coil%J_vec(n,1:3,:,:,:,:) * self%coil%A(n)
   !self%divergence(3,:,:,:,:) = 0.0_R8P
   !call self%compute_divergence(ivar=1_I4P,q=self%coil%J_vec(n,1:3,:,:,:,:),divergence=self%divergence(3,:,:,:,:))
   !print *, 'Divergenza J max pre Poisson: ', maxval(abs(self%divergence(3,:,:,:,:)))
!
   !!Applico correzione alla Poisson per ridurre divergenza residua
   !call self%impose_div_coil_correction(ivar=1_I4P, q=self%coil%J_vec(n,1:3,:,:,:,:))
   !self%divergence(3,:,:,:,:) = 0.0_R8P
   !call self%compute_divergence(ivar=1_I4P,q=self%coil%J_vec(n,1:3,:,:,:,:),divergence=self%divergence(3,:,:,:,:))
   !print *, 'Divergenza J max post Poisson: ', maxval(abs(self%divergence(3,:,:,:,:)))
!
   !!Riscalo J_vec a versore
   !self%coil%J_vec(n,1:3,:,:,:,:)=self%coil%J_vec(n,1:3,:,:,:,:)/self%coil%A(n)
!
!
   !!!Prove per capire se il problema è il round-off
   !!print *, 'max|J| pre  = ', maxval(abs(self%coil%J_vec(n,1:3,:,:,:,:)))
   !!J_vec_buffer = self%coil%J_vec(n,1:3,:,:,:,:)
   !!self%coil%J_vec(n,1:3,:,:,:,:)=self%coil%J_vec(n,1:3,:,:,:,:)/self%coil%A(n)*self%coil%A(n)
   !!print *, 'max|J| post = ', maxval(abs(self%coil%J_vec(n,1:3,:,:,:,:)))
   !!print *, 'max|dJ|     = ', maxval(abs( self%coil%J_vec(n,1:3,:,:,:,:) - J_vec_buffer(1:3,:,:,:,:) ))
   !!print *, 'rel err max = ', maxval(abs(self%coil%J_vec(n,1:3,:,:,:,:) - J_vec_buffer)) / &
   !!                           maxval(abs(J_vec_buffer))
   !!self%divergence(3,:,:,:,:) = 0.0_R8P
   !!call self%compute_divergence(ivar=1_I4P,q=self%coil%J_vec,divergence=self%divergence(3,:,:,:,:))
   !!print *, 'Divergenza J max post Poisson: ', maxval(abs(self%divergence(3,:,:,:,:)))
!
!
!
   !!Verifico che la corrente complessiva sia effettivamente quella richiesta (printo solo, non modifico)
   !call self%compute_coil_current_density_flux(n=n, adjust_amplitude=.false.)
!
   !!Calcolo le divergenze di J_vec e J come J = A*J_vec. Non dovrebbe cambiare nulla riseptto alla divergenza di Jmax
   !!precedentemente calcolata e stampata (a meno di errori legati al roundoff numerico)
   !self%divergence(3,:,:,:,:) = 0.0_R8P
   !call self%compute_divergence(ivar=1_I4P,q=self%coil%J_vec(n,1:3,:,:,:,:), &
   !                             divergence=self%divergence(3,:,:,:,:))
   !print *, 'Divergenza finale J_vec: ', maxval(abs(self%divergence(3,:,:,:,:)))
   !self%divergence(3,:,:,:,:) = 0.0_R8P
   !call self%compute_divergence(ivar=1_I4P,q=self%coil%J_vec(n,1:3,:,:,:,:)*self%coil%A(n), &
   !                             divergence=self%divergence(3,:,:,:,:))
   !print *, 'Divergenza finale corrente: ', maxval(abs(self%divergence(3,:,:,:,:)))
   endsubroutine set_rectangular_coil_z

   subroutine set_circular_coil_x(self, n, verse)
   class(prism_common_object),      intent(inout) :: self                    !< Cpu object.
   integer(I4P),                    intent(in)    :: n                       !< Coil number.
   real(R8P),                       intent(in)    :: verse                   !< Coil normal direction, +1=+x, -1=-x.
   real(R8P),                       allocatable   :: A(:,:,:,:,:)            !< Campo vettoriale della spira
   real(R8P),                       allocatable   :: J_vec_buffer(:,:,:,:,:) !< Buffer per il campo di corrente
   real(R8P)                                      :: cell_coord(3)           !< Vettore posizione centro cella
   real(R8P)                                      :: rho                     !< Distanza radiale dal centro spira nel piano yz
   real(R8P)                                      :: F_n, W_x                !< Funzioni di shaping
   real(R8P)                                      :: sigma_rho               !< Spessore radiale effettivo
   integer(I4P)                                   :: b,i,j,k                 !< Counter.

   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,                    &
             nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                          &
             y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), rc=>self%coil%r_coil(n), &
             coil_flag=>self%coil%coil_flag, dx=>field%dxyz(1,:), dy=>field%dxyz(2,:),        &
             dz=>field%dxyz(3,:), nb=>field%nb, x_cell=>field%x_cell,                         &
             y_cell=>field%y_cell, z_cell=>field%z_cell, sigma=>self%coil%sigma(n),           &
             J_vec=>self%coil%J_vec, hs=>self%fdv_half_stencil)

   allocate(A(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   A(:,:,:,:,:) = 0.0_R8P

   allocate(J_vec_buffer(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   J_vec_buffer(:,:,:,:,:) = 0.0_R8P

   do b=1, blocks_number
      sigma_rho = sigma * max(dy(b),dz(b))
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]

               rho = sqrt((cell_coord(2)-y_c)**2 + (cell_coord(3)-z_c)**2)

               F_n = erf_function(s=rho, mu=rc, sigma=sigma_rho)
               W_x = tangential_window(s=cell_coord(1), smin=x_c-sigma*dx(b), smax=x_c+sigma*dx(b), &
                                       sigma=sigma*dx(b))

               A(1,i,j,k,b) = -F_n * W_x
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

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (J_vec_buffer(1,i,j,k,b) /= 0.0_R8P .or. J_vec_buffer(2,i,j,k,b) /= 0.0_R8P .or. &
                   J_vec_buffer(3,i,j,k,b) /= 0.0_R8P) then
                  coil_flag(i,j,k,b) = n
               endif
            enddo
         enddo
      enddo
   enddo

   J_vec(1:3,:,:,:,:,n) = J_vec(1:3,:,:,:,:,n) + J_vec_buffer
   endassociate

   !if (n == 1_I4P) then
      call self%compute_coil_current_density_flux(n=n, adjust_amplitude=.true.)
   !else
   !   self%coil%A(n) = self%coil%A(1)
   !   call self%compute_coil_current_density_flux(n=n, adjust_amplitude=.false.)
   !endif

   endsubroutine set_circular_coil_x

   subroutine set_circular_coil_y(self, n, verse)
   class(prism_common_object),      intent(inout) :: self                    !< Cpu object.
   integer(I4P),                    intent(in)    :: n                       !< Coil number.
   real(R8P),                       intent(in)    :: verse                   !< Coil normal direction, +1=+y, -1=-y.
   real(R8P),                       allocatable   :: A(:,:,:,:,:)            !< Campo vettoriale della spira
   real(R8P),                       allocatable   :: J_vec_buffer(:,:,:,:,:) !< Buffer per il campo di corrente
   real(R8P)                                      :: cell_coord(3)           !< Vettore posizione centro cella
   real(R8P)                                      :: rho                     !< Distanza radiale dal centro spira nel piano xz
   real(R8P)                                      :: F_n, W_y                !< Funzioni di shaping
   real(R8P)                                      :: sigma_rho               !< Spessore radiale effettivo
   integer(I4P)                                   :: b,i,j,k                 !< Counter.

   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,                    &
             nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                          &
             y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), rc=>self%coil%r_coil(n), &
             coil_flag=>self%coil%coil_flag, dx=>field%dxyz(1,:), dy=>field%dxyz(2,:),        &
             dz=>field%dxyz(3,:), nb=>field%nb, x_cell=>field%x_cell,                         &
             y_cell=>field%y_cell, z_cell=>field%z_cell, sigma=>self%coil%sigma(n),           &
             J_vec=>self%coil%J_vec, hs=>self%fdv_half_stencil)

   allocate(A(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   A(:,:,:,:,:) = 0.0_R8P

   allocate(J_vec_buffer(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   J_vec_buffer(:,:,:,:,:) = 0.0_R8P

   do b=1, blocks_number
      sigma_rho = sigma * max(dx(b),dz(b))
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]

               rho = sqrt((cell_coord(1)-x_c)**2 + (cell_coord(3)-z_c)**2)

               F_n = erf_function(s=rho, mu=rc, sigma=sigma_rho)
               W_y = tangential_window(s=cell_coord(2), smin=y_c-sigma*dy(b), smax=y_c+sigma*dy(b), &
                                       sigma=sigma*dy(b))

               A(2,i,j,k,b) = -F_n * W_y
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

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (J_vec_buffer(1,i,j,k,b) /= 0.0_R8P .or. J_vec_buffer(2,i,j,k,b) /= 0.0_R8P .or. &
                   J_vec_buffer(3,i,j,k,b) /= 0.0_R8P) then
                  coil_flag(i,j,k,b) = n
               endif
            enddo
         enddo
      enddo
   enddo

   J_vec(1:3,:,:,:,:,n) = J_vec(1:3,:,:,:,:,n) + J_vec_buffer
   endassociate

   !if (n == 1_I4P) then
      call self%compute_coil_current_density_flux(n=n, adjust_amplitude=.true.)
   !else
   !   self%coil%A(n) = self%coil%A(1)
   !   call self%compute_coil_current_density_flux(n=n, adjust_amplitude=.false.)
   !endif

   endsubroutine set_circular_coil_y

   subroutine set_circular_coil_z(self, n, verse)
   class(prism_common_object),      intent(inout) :: self                    !< Cpu object.
   integer(I4P),                    intent(in)    :: n                       !< Coil number.
   real(R8P),                       intent(in)    :: verse                   !< Coil normal direction, +1=+z, -1=-z.
   real(R8P),                       allocatable   :: A(:,:,:,:,:)            !< Campo vettoriale della spira
   real(R8P),                       allocatable   :: J_vec_buffer(:,:,:,:,:) !< Buffer per il campo di corrente
   real(R8P)                                      :: c_c(3)                  !< Vettore posizione centro spira
   real(R8P)                                      :: cell_coord(3)           !< Vettore posizione centro cella
   real(R8P)                                      :: rho                     !< Distanza radiale dal centro spira nel piano xy
   real(R8P)                                      :: F_n, W_z                !< Funzioni di shaping
   real(R8P)                                      :: sigma_rho               !< Spessore radiale effettivo
   integer(I4P)                                   :: b,i,j,k                 !< Counter.

   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,                    &
             nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                          &
             y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), rc=>self%coil%r_coil(n), &
             coil_flag=>self%coil%coil_flag, dx=>field%dxyz(1,:), dy=>field%dxyz(2,:),        &
             dz=>field%dxyz(3,:), nb=>field%nb, x_cell=>field%x_cell,                         &
             y_cell=>field%y_cell, z_cell=>field%z_cell, sigma=>self%coil%sigma(n),           &
             J_vec=>self%coil%J_vec, hs=>self%fdv_half_stencil)

   allocate(A(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   A(:,:,:,:,:) = 0.0_R8P

   allocate(J_vec_buffer(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   J_vec_buffer(:,:,:,:,:) = 0.0_R8P

   do b=1, blocks_number
      sigma_rho = sigma * max(dx(b),dy(b))
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]

               rho = sqrt((cell_coord(1)-x_c)**2 + (cell_coord(2)-y_c)**2)

               F_n = erf_function(s=rho, mu=rc, sigma=sigma_rho)
               W_z = tangential_window(s=cell_coord(3), smin=z_c-sigma*dz(b), smax=z_c+sigma*dz(b), &
                                       sigma=sigma*dz(b))

               A(3,i,j,k,b) = -F_n * W_z
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

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (J_vec_buffer(1,i,j,k,b) /= 0.0_R8P .or. J_vec_buffer(2,i,j,k,b) /= 0.0_R8P .or. &
                   J_vec_buffer(3,i,j,k,b) /= 0.0_R8P) then
                  coil_flag(i,j,k,b) = n
               endif
            enddo
         enddo
      enddo
   enddo

   J_vec(1:3,:,:,:,:,n) = J_vec(1:3,:,:,:,:,n) + J_vec_buffer

   endassociate

   !if (n == 1_I4P) then
      call self%compute_coil_current_density_flux(n=n, adjust_amplitude=.true.)
   !else
   !   self%coil%A(n) = self%coil%A(1)
   !   call self%compute_coil_current_density_flux(n=n, adjust_amplitude=.false.)
   !endif

   endsubroutine set_circular_coil_z

   subroutine set_solenoid_x(self, n, verse)
   class(prism_common_object),      intent(inout) :: self                    !< Cpu object.
   integer(I4P),                    intent(in)    :: n                       !< Coil number.
   real(R8P),                       intent(in)    :: verse                   !< Solenoid normal direction, +1=+x, -1=-x.
   real(R8P),                       allocatable   :: A(:,:,:,:,:)            !< Campo vettoriale del solenoide
   real(R8P),                       allocatable   :: J_vec_buffer(:,:,:,:,:) !< Variabile buffer per il campo di corrente
   real(R8P)                                      :: cell_coord(3)           !< Vettore posizione centro cella
   real(R8P)                                      :: rho                     !< Distanza radiale nel piano yz
   real(R8P)                                      :: F_n, W_x                !< Shaping functions
   real(R8P)                                      :: sigma_rho               !< Spessore radiale effettivo
   integer(I4P)                                   :: b,i,j,k                 !< Counter.

   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,                       &
            nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                              &
            y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), r_coil=>self%coil%r_coil(n), &
            l_sol=>self%coil%l_solenoid(n), coil_flag=>self%coil%coil_flag,                      &
            dx=>field%dxyz(1,:), dy=>field%dxyz(2,:), dz=>field%dxyz(3,:),                       &
            nb=>field%nb, x_cell=>field%x_cell, y_cell=>field%y_cell,                            &
            z_cell=>field%z_cell, sigma=>self%coil%sigma(n), J_vec=>self%coil%J_vec,             &
            hs=>self%fdv_half_stencil)

   allocate(A(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   A(:,:,:,:,:) = 0.0_R8P

   allocate(J_vec_buffer(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   J_vec_buffer(:,:,:,:,:) = 0.0_R8P

   do b=1, blocks_number
      sigma_rho = sigma * max(dy(b),dz(b))
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]

               rho = sqrt((cell_coord(2)-y_c)**2 + (cell_coord(3)-z_c)**2)

               F_n = erf_function(s=rho, mu=r_coil, sigma=sigma_rho)
               W_x = tangential_window(s=cell_coord(1), smin=x_c-0.5_R8P*l_sol, smax=x_c+0.5_R8P*l_sol, &
                                       sigma=sigma*dx(b))

               A(1,i,j,k,b) = - F_n * W_x
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

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (J_vec_buffer(1,i,j,k,b) /= 0.0_R8P .or. J_vec_buffer(2,i,j,k,b) /= 0.0_R8P .or. &
                   J_vec_buffer(3,i,j,k,b) /= 0.0_R8P) then
                  coil_flag(i,j,k,b) = n
               endif
            enddo
         enddo
      enddo
   enddo

   J_vec(1:3,:,:,:,:,n) = J_vec(1:3,:,:,:,:,n) + J_vec_buffer
   endassociate

   !if (n == 1_I4P) then
      call self%compute_solenoid_current_density_flux(n=n, adjust_amplitude=.true.)
   !else
   !   self%coil%A(n) = self%coil%A(1)
   !   call self%compute_coil_current_density_flux(n=n, adjust_amplitude=.false.)
   !endif
   endsubroutine set_solenoid_x

   subroutine set_solenoid_y(self, n, verse)
   class(prism_common_object),      intent(inout) :: self                    !< Cpu object.
   integer(I4P),                    intent(in)    :: n                       !< Coil number.
   real(R8P),                       intent(in)    :: verse                   !< Solenoid normal direction, +1=+y, -1=-y.
   real(R8P),                       allocatable   :: A(:,:,:,:,:)            !< Campo vettoriale del solenoide
   real(R8P),                       allocatable   :: J_vec_buffer(:,:,:,:,:) !< Variabile buffer per il campo di corrente
   real(R8P)                                      :: cell_coord(3)           !< Vettore posizione centro cella
   real(R8P)                                      :: rho                     !< Distanza radiale nel piano xz
   real(R8P)                                      :: F_n, W_y                !< Shaping functions
   real(R8P)                                      :: sigma_rho               !< Spessore radiale effettivo
   integer(I4P)                                   :: b,i,j,k                 !< Counter.

   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,                       &
            nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                              &
            y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), r_coil=>self%coil%r_coil(n), &
            l_sol=>self%coil%l_solenoid(n), coil_flag=>self%coil%coil_flag,                      &
            dx=>field%dxyz(1,:), dy=>field%dxyz(2,:), dz=>field%dxyz(3,:),                       &
            nb=>field%nb, x_cell=>field%x_cell, y_cell=>field%y_cell,                            &
            z_cell=>field%z_cell, sigma=>self%coil%sigma(n), J_vec=>self%coil%J_vec,             &
            hs=>self%fdv_half_stencil)

   allocate(A(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   A(:,:,:,:,:) = 0.0_R8P

   allocate(J_vec_buffer(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   J_vec_buffer(:,:,:,:,:) = 0.0_R8P

   do b=1, blocks_number
      sigma_rho = sigma * max(dx(b),dz(b))
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]

               rho = sqrt((cell_coord(1)-x_c)**2 + (cell_coord(3)-z_c)**2)

               F_n = erf_function(s=rho, mu=r_coil, sigma=sigma_rho)
               W_y = tangential_window(s=cell_coord(2), smin=y_c-0.5_R8P*l_sol, smax=y_c+0.5_R8P*l_sol, &
                                       sigma=sigma*dy(b))

               A(2,i,j,k,b) = - F_n * W_y
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

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (J_vec_buffer(1,i,j,k,b) /= 0.0_R8P .or. J_vec_buffer(2,i,j,k,b) /= 0.0_R8P .or. &
                   J_vec_buffer(3,i,j,k,b) /= 0.0_R8P) then
                  coil_flag(i,j,k,b) = n
               endif
            enddo
         enddo
      enddo
   enddo

   J_vec(1:3,:,:,:,:,n) = J_vec(1:3,:,:,:,:,n) + J_vec_buffer
   endassociate

   !if (n == 1_I4P) then
      call self%compute_solenoid_current_density_flux(n=n, adjust_amplitude=.true.)
   !else
   !   self%coil%A(n) = self%coil%A(1)
   !   call self%compute_coil_current_density_flux(n=n, adjust_amplitude=.false.)
   !endif
   endsubroutine set_solenoid_y

   subroutine set_solenoid_z(self, n, verse)
   class(prism_common_object),      intent(inout) :: self                    !< Cpu object.
   integer(I4P),                    intent(in)    :: n                       !< Coil number.
   real(R8P),                       intent(in)    :: verse                   !< Solenoid normal direction, +1=+z, -1=-z.
   real(R8P),                       allocatable   :: A(:,:,:,:,:)            !< Campo vettoriale del solenoide
   real(R8P),                       allocatable   :: J_vec_buffer(:,:,:,:,:) !< Variabile buffer per il campo di corrente
   real(R8P)                                      :: cell_coord(3)           !< Vettore posizione centro cella
   real(R8P)                                      :: rho                     !< Distanza radiale nel piano xy
   real(R8P)                                      :: F_n, W_z                !< Shaping functions
   real(R8P)                                      :: sigma_rho               !< Spessore radiale effettivo
   integer(I4P)                                   :: b,i,j,k                 !< Counter.

   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,                       &
            nk=>grid%nk, ngc=>grid%ngc, x_c=>self%coil%x_center(n),                              &
            y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), r_coil=>self%coil%r_coil(n), &
            l_sol=>self%coil%l_solenoid(n), coil_flag=>self%coil%coil_flag,                      &
            dx=>field%dxyz(1,:), dy=>field%dxyz(2,:), dz=>field%dxyz(3,:),                       &
            nb=>field%nb, x_cell=>field%x_cell, y_cell=>field%y_cell,                            &
            z_cell=>field%z_cell, sigma=>self%coil%sigma(n), J_vec=>self%coil%J_vec,             &
            hs=>self%fdv_half_stencil)

   allocate(A(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   A(:,:,:,:,:) = 0.0_R8P

   allocate(J_vec_buffer(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   J_vec_buffer(:,:,:,:,:) = 0.0_R8P

   do b=1, blocks_number
      sigma_rho = sigma * max(dx(b),dy(b))
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]

               rho = sqrt((cell_coord(1)-x_c)**2 + (cell_coord(2)-y_c)**2)

               F_n = erf_function(s=rho, mu=r_coil, sigma=sigma_rho)
               W_z = tangential_window(s=cell_coord(3), smin=z_c-0.5_R8P*l_sol, smax=z_c+0.5_R8P*l_sol, &
                                       sigma=sigma*dz(b))

               A(3,i,j,k,b) = - F_n * W_z
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

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (J_vec_buffer(1,i,j,k,b) /= 0.0_R8P .or. J_vec_buffer(2,i,j,k,b) /= 0.0_R8P .or. &
                   J_vec_buffer(3,i,j,k,b) /= 0.0_R8P) then
                  coil_flag(i,j,k,b) = n
               endif
            enddo
         enddo
      enddo
   enddo

   J_vec(1:3,:,:,:,:,n) = J_vec(1:3,:,:,:,:,n) + J_vec_buffer
   endassociate

   !if (n == 1_I4P) then
      call self%compute_solenoid_current_density_flux(n=n, adjust_amplitude=.true.)
   !else
   !   self%coil%A(n) = self%coil%A(1)
   !   call self%compute_coil_current_density_flux(n=n, adjust_amplitude=.false.)
   !endif
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
endmodule adam_prism_common_object
