!< ADAM, Maxwell equations system class definition, common data to all backends.
module adam_prism_common_object

! ADAM modules
use adam_common_library
! PRISM modules
use adam_prism_bc_object
use adam_prism_coil_object
use adam_prism_external_fields_object
use adam_prism_fWLayer_object
use adam_prism_ic_object
use adam_prism_io_object
use adam_prism_leapfrog_pic_object
use adam_prism_numerics_object
use adam_prism_physics_object
use adam_prism_pic_object
use adam_prism_particle_injection_object
use adam_prism_rk_pic_object
use adam_prism_rk_bc_object
use adam_prism_time_object
! third party modules
use penf
! use ISO_C_BINDING

implicit none
private
public :: prism_common_object

type :: prism_common_object
   !< Maxwell equations system class definition, common data to all backends.
   ! ADAM library objects
   type(mpih_object)           :: mpih          !< MPI handler.
   type(adam_object)           :: adam          !< ADAM.
   type(field_object), pointer :: field=>null() !< The field.
   type(grid_object),  pointer :: grid=>null()  !< The grid.
   type(amr_object)            :: amr           !< AMR marker handler.
   type(ib_object)             :: ib            !< Immersed Boundary (IB) handler.
   type(slices_object)         :: slices        !< Slices handler.
   type(blanesmoan_object)     :: blanesmoan    !< Blanes-Moan integrator.
   type(cfm_object)            :: cfm           !< Commutator-Free Magnus integrator.
   type(leapfrog_object)       :: leapfrog      !< Leapfrog integrator.
   type(rk_object)             :: rk            !< RK integrator.
   type(weno_object)           :: weno          !< WENO reconstructor.
   type(flail_object)          :: flail         !< Linear algebra methods handler.
   ! PRISM library objects
   type(prism_io_object)                 :: io                 !< IO handler.
   type(prism_numerics_object)           :: numerics           !< Numerics handler.
   type(prism_physics_object)            :: physics            !< Fluids physiscs handler.
   type(prism_ic_object)                 :: ic                 !< Initial Conditions (IC) handler.
   type(prism_bc_object)                 :: bc                 !< Boundary Conditions (BC) handler.
   type(prism_rk_bc_object)              :: rk_bc              !< RK integrator for BC.
   type(prism_time_object)               :: time               !< Time handler.
   type(prism_fWLayer_object)            :: fWLayer            !< fWLayer handler.
   type(prism_coil_object)               :: coil               !< Coils handler.
   type(prism_external_fields_object)    :: external_fields    !< External fields handler.
   type(prism_pic_object)                :: pic                !< Particle-in-Cell (PIC) handler.
   type(prism_particle_injection_object) :: particle_injection !< Particle injection handler.
   type(prism_leapfrog_pic_object)       :: leapfrog_pic       !< Leapfrog PIC integrator.
   type(prism_rk_pic_object)             :: rk_pic             !< RK PIC integrator.
   ! grid/field data replica for easy handling
   integer(I4P), pointer :: ngc=>null()           !< Number of ghost cells.
   integer(I4P), pointer :: ni=>null()            !< Number of cells in i direction.
   integer(I4P), pointer :: nj=>null()            !< Number of cells in j direction.
   integer(I4P), pointer :: nk=>null()            !< Number of cells in k direction.
   integer(I4P), pointer :: nb=>null()            !< Total blocks number for MPI.
   integer(I4P), pointer :: blocks_number=>null() !< Actual blocks number.
   integer(I4P), pointer :: nv=>null()            !< Number of variables in q vector.
   integer(I4P), pointer :: nv_c=>null()          !< Number of conservative variables in q vector.
   integer(I4P), pointer :: nv_s=>null()          !< Number of source variables in q vector.
   integer(I4P), pointer :: nv_cl=>null()         !< Number of divergence cleaning variables in q vector.
   ! fields data [1:nv,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb].
   real(R8P),    allocatable :: q(         :,:,:,:,:) !< Conservative cell centered variables.
   real(R8P),    allocatable :: dq(        :,:,:,:,:) !< Residuals right hand side.
   real(R8P),    allocatable :: q_pic(           :,:) !< PIC variables.
   real(R8P),    allocatable :: pic_fields(      :,:) !< Fields value at particle locations.
   real(R8P),    allocatable :: curl(      :,:,:,:,:) !< Curl fields.
   real(R8P),    allocatable :: divergence(:,:,:,:,:) !< Divergence fields.
   character(3), allocatable :: q_name(:)             !< Fields names [1:nv].
   ! auxiliary data
   real(R8P), allocatable :: energy_D(:)                !< Energy of field D, time history.
   real(R8P), allocatable :: energy_B(:)                !< Energy of field B, time history.
   real(R8P), allocatable :: coil_power(:)              !< Power of coils, time history.
   real(R8P), allocatable :: Poynting_flux(:)           !< Total Poynting flux from boundary, time history.
   real(R8P)              :: rms_energy_error_D=0.0_R8P !< RMS energy error of D field.
   real(R8P)              :: rms_energy_error_B=0.0_R8P !< RMS energy error of B field.
   !< Pointer (abstract) TBP.
   procedure(compute_curl_interface),       pass(self),pointer :: compute_curl       =>null()!< Compute curl of vector field.
   procedure(compute_derivative1_interface),pass(self),pointer :: compute_derivative1=>null()!< Compute derivative1 of scalar field.
   procedure(compute_derivative2_interface),pass(self),pointer :: compute_derivative2=>null()!< Compute derivative2 of scalar field.
   procedure(compute_derivative4_interface),pass(self),pointer :: compute_derivative4=>null()!< Compute derivative4 of scalar field.
   procedure(compute_divergence_interface), pass(self),pointer :: compute_divergence =>null()!< Compute divergence of vector field.
   procedure(compute_gradient_interface),   pass(self),pointer :: compute_gradient   =>null()!< Compute gradient of scalar field.
   procedure(compute_laplacian_interface),  pass(self),pointer :: compute_laplacian  =>null()!< Compute laplacian of scalar field.
   contains
      procedure, pass(self) :: allocate_common         !< Allocate common data.
      procedure, pass(self) :: initialize_common       !< Initialize the equation common data.
      ! IO methods
      procedure, pass(self) :: load_restart_files      !< Load restart files.
      procedure, pass(self) :: save_energy_error       !< Save energy error history.
      procedure, pass(self) :: save_energy_history     !< Save energy history.
      procedure, pass(self) :: save_divergence_history !< Save divergence history.
      procedure, pass(self) :: save_restart_files      !< Save restart files.
      procedure, pass(self) :: save_xh5f               !< Save simulation data in XH5F format.
      ! FDV operators numerical methods
      procedure, pass(self) :: compute_curl_fd        !< Compute curl of vector field by finite difference.
      procedure, pass(self) :: compute_curl_fv        !< Compute curl of vector field by finite volume.
      procedure, pass(self) :: compute_derivative1_fd !< Compute derivative1 of scalar fields, finite difference schemes.
      procedure, pass(self) :: compute_derivative1_fv !< Compute derivative1 of scalar fields, finite volume schemes.
      procedure, pass(self) :: compute_derivative2_fd !< Compute derivative2 of scalar fields, finite difference schemes.
      procedure, pass(self) :: compute_derivative2_fv !< Compute derivative2 of scalar fields, finite volume schemes.
      procedure, pass(self) :: compute_derivative4_fd !< Compute derivative4 of scalar fields, finite difference schemes.
      procedure, pass(self) :: compute_divergence_fd  !< Compute divergence of vector field by finite difference.
      procedure, pass(self) :: compute_divergence_fv  !< Compute divergence of vector field by finite volume.
      procedure, pass(self) :: compute_gradient_fd    !< Compute gradient of scalar field, finite difference schemes.
      procedure, pass(self) :: compute_gradient_fv    !< Compute gradient of scalar field, finite volume schemes.
      procedure, pass(self) :: compute_laplacian_fd   !< Compute laplacian of scalar field, finite difference schemes.
      procedure, pass(self) :: compute_laplacian_fv   !< Compute laplacian of scalar field, finite volume schemes.
      ! coils initialization methods
      procedure, pass(self) :: compute_coil_current_density_flux !< Compute coil current density fluxes for Maxwell equations.
      procedure, pass(self) :: initialize_coils                  !< Initialize coils.
      procedure, pass(self) :: set_rectangular_coil_x            !< Subroutine to set a rectangular coil source with +-x normal
      procedure, pass(self) :: set_rectangular_coil_y            !< Subroutine to set a rectangular coil source with +-y normal
      procedure, pass(self) :: set_rectangular_coil_z            !< Subroutine to set a rectangular coil source with +-z normal
endtype prism_common_object

interface
   subroutine compute_curl_interface(self, ivar, q, curl)
   !< Compute curl of vector fields, div(q(ivar:ivar+2).
   import :: prism_common_object, I4P, R8P
   class(prism_common_object), intent(in)    :: self                                            !< The equation.
   integer(I4P),               intent(in)    :: ivar                                            !< Start index of variable of q.
   real(R8P),                  intent(in)    :: q(   1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),                  intent(inout) :: curl(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Curl.
   endsubroutine compute_curl_interface

   subroutine compute_derivative1_interface(self, dir, ivar, q, dq_ds)
   !< Compute derivative1 of scalar fields, dq(ivar)/ds.
   import :: prism_common_object, I4P, R8P
   class(prism_common_object), intent(in)    :: self                                         !< The equation.
   integer(I4P),               intent(in)    :: dir                                          !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),               intent(in)    :: ivar                                         !< Start index of (vec.) variable of q.
   real(R8P),                  intent(in)    :: q(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Field variables.
   real(R8P),                  intent(inout) :: dq_ds(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Derivative1, dq/ds.
   endsubroutine compute_derivative1_interface

   subroutine compute_derivative2_interface(self, dir, ivar, q, d2q_ds2)
   !< Compute derivative2 of scalar fields, d2q(ivar)/ds2.
   import :: prism_common_object, I4P, R8P
   class(prism_common_object), intent(in)    :: self                                            !< The equation.
   integer(I4P),               intent(in)    :: dir                                             !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),               intent(in)    :: ivar                                            !< Start index of variable of q.
   real(R8P),                  intent(in)    :: q(1:,   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),                  intent(inout) :: d2q_ds2(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative2, d2q/ds2.
   endsubroutine compute_derivative2_interface

   subroutine compute_derivative4_interface(self, dir, ivar, q, d4q_ds4)
   !< Compute derivative4 of scalar fields, d4q(ivar)/ds4.
   import :: prism_common_object, I4P, R8P
   class(prism_common_object), intent(in)    :: self                                            !< The equation.
   integer(I4P),               intent(in)    :: dir                                             !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),               intent(in)    :: ivar                                            !< Start index of variable of q.
   real(R8P),                  intent(in)    :: q(1:,   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),                  intent(inout) :: d4q_ds4(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative4, d4q/ds4.
   endsubroutine compute_derivative4_interface

   subroutine compute_divergence_interface(self, ivar, q, divergence)
   !< Compute divergence of vector fields, div(q(ivar:ivar+2).
   import :: prism_common_object, I4P, R8P
   class(prism_common_object), intent(in)    :: self                                               !< The equation.
   integer(I4P),               intent(in)    :: ivar                                               !< Start index of field of q.
   real(R8P),                  intent(in)    :: q(1:,      1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),                  intent(inout) :: divergence(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
   endsubroutine compute_divergence_interface

   subroutine compute_gradient_interface(self, ivar, q, gradient)
   !< Compute gradient of scalar variable q(ivar).
   import :: prism_common_object, I4P, R8P
   class(prism_common_object), intent(in)    :: self                                               !< The equation.
   integer(I4P),               intent(in)    :: ivar                                               !< Index of scalar variable of q.
   real(R8P),                  intent(in)    :: q(       1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Field variables.
   real(R8P),                  intent(inout) :: gradient(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Gradient.
   endsubroutine compute_gradient_interface

   subroutine compute_laplacian_interface(self, ivar, q, laplacian)
   !< Compute laplacian of scalar variable q(ivar).
   import :: prism_common_object, I4P, R8P
   class(prism_common_object), intent(in)    :: self                                              !< The equation.
   integer(I4P),               intent(in)    :: ivar                                              !< Index of scalar variable of q.
   real(R8P),                  intent(in)    :: q(     1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),                  intent(inout) :: laplacian(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Laplacian.
   endsubroutine compute_laplacian_interface
endinterface
contains
   subroutine allocate_common(self)
   !< Allocate common data.
   class(prism_common_object), intent(inout) :: self !< The equation.

   associate(nv=>self%nv, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, nb=>self%nb, &
             particle_number=>self%pic%particle_number)
   call allocate_variable(var=self%dq,               &
                          ulb=reshape([1,nv,         &
                                       1-ngc,ni+ngc, &
                                       1-ngc,nj+ngc, &
                                       1-ngc,nk+ngc, &
                                       1,nb],[2,5]), &
                          msg=self%mpih%myrankstr//'prism_common_object%allocate_common(dq) ', verbose=.true.)
   self%dq = 0._R8P
   call allocate_variable(var=self%divergence,        &
                          ulb=reshape([1,self%nv,    &
                                       1-ngc,ni+ngc, &
                                       1-ngc,nj+ngc, &
                                       1-ngc,nk+ngc, &
                                       1,nb],[2,5]), &
                          msg=self%mpih%myrankstr//'prism_common_object%allocate_common(divergence) ', verbose=.true.)
   self%divergence = 0._R8P
   self%divergence(4,:,:,:,:) = self%fWLayer%f(1,:,:,:,:)
   self%divergence(5,:,:,:,:) = self%fWLayer%f(2,:,:,:,:)
   self%divergence(6,:,:,:,:) = self%fWLayer%f(3,:,:,:,:)
   call allocate_variable(var=self%curl,             &
                          ulb=reshape([1,self%nv,    &
                                       1-ngc,ni+ngc, &
                                       1-ngc,nj+ngc, &
                                       1-ngc,nk+ngc, &
                                       1,nb],[2,5]), &
                          msg=self%mpih%myrankstr//'prism_common_object%allocate_common(curl) ', verbose=.true.)
   self%curl = 0._R8P

   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
      call allocate_variable(var=self%q_pic,                          &
                             ulb=reshape([1,8,                        &
                                          1,particle_number],[2,2]),  &
                             msg=self%mpih%myrankstr//'prism_common_object%allocate_common(q_pic) ', verbose=.true.)
      self%q_pic = 0._R8P
      call allocate_variable(var=self%pic_fields,                     &
                             ulb=reshape([1,6,                        &
                                          1,particle_number],[2,2]),  &
                             msg=self%mpih%myrankstr//'prism_common_object%allocate_common(pic_fields) ', verbose=.true.)
      self%pic_fields = 0._R8P
   endif
   endassociate
   endsubroutine allocate_common

   subroutine initialize_common(self, field, filename, memory_avail, do_mpi_init, verbose)
   !< Initialize the equation common data.
   class(prism_common_object), intent(inout), target :: self         !< The equation.
   type(field_object),         intent(inout)         :: field        !< The field.
   character(*),               intent(in)            :: filename     !< Input file name.
   real(R8P),                  intent(in),value      :: memory_avail !< Memory available for single MPI process.
   logical,                    intent(in), optional  :: do_mpi_init  !< Flag to activate MPI init call.
   logical,                    intent(in), optional  :: verbose      !< Trigger verbose output.
   logical                                           :: verbose_     !< Trigger verbose output, local variable.
   integer(I8P)                                      :: nodes_number !< Allocated nodes on tree.
   integer(I4P)                                      :: nb           !< Number of allocated blocks.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   call self%mpih%initialize(do_mpi_init=do_mpi_init, verbose=verbose_)
   if (verbose_) call self%mpih%print_message('prism_common_object%initialize start')
   call self%io%initialize(filename=trim(filename))
   associate(file_parameters=>self%io%file_parameters)
   call self%bc%initialize(file_parameters=file_parameters)
   call self%numerics%initialize(file_parameters=file_parameters)
   call self%physics%initialize(file_parameters=file_parameters, reconstruction_vars=self%numerics%reconstruction_vars, &
                                 div_corr_var=self%numerics%div_corr_var,                                               &
                                 constrained_transport_D=self%numerics%constrained_transport_D,                         &
                                 constrained_transport_B=self%numerics%constrained_transport_B)
   call self%adam%grid%initialize(file_parameters=file_parameters,bc_type=self%bc%bc_type, verbose=.true.)
   call self%adam%compute_blocks_number(memory_avail=memory_avail, fields_number=80, nb=nb, nodes_number=nodes_number)
   call self%adam%initialize(file_parameters=file_parameters, &
                             do_tree_init=.true.,             &
                             do_maps_init=.true.,             &
                             do_field_init=.true.,            &
                             nv=self%physics%nv, nb=1, nodes_number=11_I8P, q=self%q) !nb = nb !nodes_number = nodes_number
   call associate_adam_data(grid=self%adam%grid, field=self%adam%field, physics=self%physics)
   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) call self%pic%initialize(file_parameters=file_parameters, &
                                                                                   field=self%field)
   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) &
      call self%particle_injection%initialize(file_parameters=file_parameters, pic=self%pic)
   call self%adam%refine_uniform(refinement_levels=self%adam%tree%iu_ref_levels, do_blocks_reorder=.false.,q=self%q)
   call self%adam%prune(ijkl_prune=self%adam%tree%ijkl_prune, do_blocks_reorder=.false.,q=self%q)
   call self%amr%initialize(file_parameters=file_parameters)
   call self%field%compute_metrics
   call self%time%initialize(file_parameters=file_parameters)
   call self%ic%initialize(file_parameters=file_parameters)
   call self%fWLayer%initialize(file_parameters=file_parameters, physics=self%physics, field=self%field)
   call self%coil%initialize(file_parameters=file_parameters, field=self%field)
   call self%ib%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field)
   call self%slices%initialize(file_parameters=file_parameters)
   call self%external_fields%initialize(file_parameters=file_parameters)
   if (self%numerics%scheme_time==NUM_SCHEME_TIME_RUNGE_KUTTA) &
      call self%rk%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field)
   if (self%numerics%scheme_time==NUM_SCHEME_TIME_RUNGE_KUTTA) &
      call self%rk_bc%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field, rk=self%rk, &
                              physics=self%physics)
   if (self%numerics%scheme_time==NUM_SCHEME_TIME_LEAPFROG) &
      call self%leapfrog%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field)
   if (self%numerics%scheme_time==NUM_SCHEME_TIME_BLANES_MOAN) &
      call self%blanesmoan%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field)
   if (self%numerics%scheme_time==NUM_SCHEME_TIME_CFM) &
      call self%cfm%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field)
   if (self%numerics%scheme_space==NUM_SCHEME_SPACE_WENO) &
      call self%weno%initialize(file_parameters=file_parameters, nb=self%nb, ngc=self%ngc, ni=self%ni, nj=self%nj, nk=self%nk)
   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
      if (self%pic%scheme_time==NUM_SCHEME_TIME_PIC_LEAPFROG) &
         call self%leapfrog_pic%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field, pic=self%pic)
      if (self%pic%scheme_time==NUM_SCHEME_TIME_PIC_RUNGE_KUTTA) &
         call self%rk_pic%initialize(file_parameters=file_parameters, rk=self%rk, pic=self%pic)
   endif
   call self%flail%initialize(file_parameters=file_parameters)
   call check_ngc_number
   call self%allocate_common
   call io_initialize

   select case(self%numerics%scheme_space)
   case(NUM_SCHEME_SPACE_WENO)
      self%compute_curl        => compute_curl_fv
      self%compute_derivative1 => compute_derivative1_fv
      self%compute_derivative2 => compute_derivative2_fv
      self%compute_divergence  => compute_divergence_fv
      self%compute_gradient    => compute_gradient_fv
      self%compute_laplacian   => compute_laplacian_fv
   case(NUM_SCHEME_SPACE_FD_CENTERED)
      self%compute_curl        => compute_curl_fd
      self%compute_derivative1 => compute_derivative1_fd
      self%compute_derivative2 => compute_derivative2_fd
      self%compute_divergence  => compute_divergence_fd
      self%compute_gradient    => compute_gradient_fd
      self%compute_laplacian   => compute_laplacian_fd
   case(NUM_SCHEME_SPACE_FV_CENTERED)
      self%compute_curl        => compute_curl_fv
      self%compute_derivative1 => compute_derivative1_fv
      self%compute_derivative2 => compute_derivative2_fv
      self%compute_divergence  => compute_divergence_fv
      self%compute_gradient    => compute_gradient_fv
      self%compute_laplacian   => compute_laplacian_fv
   endselect

   endassociate
   if (verbose_) call self%mpih%print_message('prism_common_object%initialize finish')
   contains
      subroutine associate_adam_data(grid, field, physics)
      !< Associate objects data to equation for easy handling.
      type(grid_object),          intent(in), target :: grid    !< The grid.
      type(field_object),         intent(in), target :: field   !< The field.
      type(prism_physics_object), intent(in), target :: physics !< The physics.

      self%grid          => grid
      self%field         => field
      self%blocks_number => field%blocks_number
      self%ni            => field%grid%ni
      self%nj            => field%grid%nj
      self%nk            => field%grid%nk
      self%ngc           => field%grid%ngc
      self%nb            => field%nb
      self%nv            => physics%nv
      self%nv_c          => physics%nv_c
      self%nv_s          => physics%nv_s
      self%nv_cl         => physics%nv_cl
      !self%nv_pic        => physics%nv_pic
      endsubroutine associate_adam_data

      subroutine check_ngc_number
      !< Check if the number of ghost cells is consistent with the numerical schemes used, if not an error is echoed and
      !< the simulation is stop.

      if (self%numerics%scheme_space==NUM_SCHEME_SPACE_WENO) then
         if (self%weno%S > self%grid%ngc) &
            call self%mpih%error_stop(msg=': ghost cells number (ngc) must be >= of weno stencil number (weno%S):'//&
                                      ' ngc='//trim(str(self%grid%ngc))//' weno%S='//trim(str(self%weno%S)))
      endif
      if (self%numerics%fdv_half_stencil > self%grid%ngc) &
         call self%mpih%error_stop(msg=': ghost cells number (ngc) must be >= of FDV half stencil number (numerics%fdv_hs):'//&
                                  ' ngc='//trim(str(self%grid%ngc))//' numerics%fdv_hs='//trim(str(self%numerics%fdv_half_stencil)))
      endsubroutine check_ngc_number

      subroutine io_initialize
      !< Initialize IO data.
      character(6),  allocatable :: q1_R8P_name(:) !< Variables names buffer.
      character(5),  allocatable :: q2_R8P_name(:) !< Variables names buffer.
      character(7),  allocatable :: q3_R8P_name(:) !< Variables names buffer.
      character(7),  allocatable :: q4_R8P_name(:) !< Variables names buffer.
      integer(I4P)               :: c              !< Counter.

      call self%adam%io%initialize(grid=self%adam%grid, field=self%adam%field)
      if (self%physics%physical_model == EM_PHYSICAL_MODEL) then
         select case(self%numerics%div_corr_var)
         case(DIV_CORR_VAR_POISS)
            self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ']
            q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_Jx','res_Jy','res_Jz']
            q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09']
         case(DIV_CORR_VAR_HYPER)
            if (self%numerics%constrained_transport_D .and. .not.self%numerics%constrained_transport_B) then
               self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','phi','Jx ','Jy ','Jz ']
               q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_ph','res_Jx','res_Jy','res_Jz']
               q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09','div10']
            elseif (.not.self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
               self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','psi','Jx ','Jy ','Jz ']
               q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_ps','res_Jx','res_Jy','res_Jz']
               q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09','div10']
            elseif (self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
               self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','phi','psi','Jx ','Jy ','Jz ']
               q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_ph','res_ps','res_Jx','res_Jy','res_Jz']
               q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09','div10','div11']
            endif
         case default
            self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ']
            q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_Jx','res_Jy','res_Jz']
            q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09']
         endselect
         q4_R8P_name = ['curlD_x','curlD_y','curlD_z','curlB_x','curlB_y','curlB_z','curlJ_x','curlJ_y','curlJ_z']
         if (self%io%save_residual_fields) call self%adam%io%register_aux_field(q1_R8P=self%dq,q1_R8P_name=q1_R8P_name)
         if (self%io%save_divergence_fields) call self%adam%io%register_aux_field(q2_R8P=self%divergence,q2_R8P_name=q2_R8P_name)
         if (self%coil%total_coils_number>0) then
            q3_R8P_name = ['j_vec_1','j_vec_2','j_vec_3','f_Gauss']
            call self%adam%io%register_aux_field(q3_R8P=self%coil%j_vec(:,:,:,:,:,1),q3_R8P_name=q3_R8P_name)
            call self%adam%io%register_aux_field(s1_I4P=self%coil%coil_flag,s1_I4P_name='coil_flag')
         endif
         if (self%io%save_curl_fields) call self%adam%io%register_aux_field(q4_R8P=self%curl,q4_R8P_name=q4_R8P_name)
      elseif(self%physics%physical_model == PIC_PHYSICAL_MODEL) then
         select case(self%numerics%div_corr_var)
         case(DIV_CORR_VAR_POISS)
            self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ','rho']
            q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_Jx','res_Jy','res_Jz','res_rh']
            q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09','div10']
         case(DIV_CORR_VAR_HYPER)
            if (self%numerics%constrained_transport_D .and. .not.self%numerics%constrained_transport_B) then
               self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','phi','Jx ','Jy ','Jz ','rho']
               q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_ph','res_Jx','res_Jy','res_Jz','res_rh']
               q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09','div10','div11']
            elseif (.not.self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
               self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','psi','Jx ','Jy ','Jz ','rho']
               q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_ps','res_Jx','res_Jy','res_Jz','res_rh']
               q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09','div10','div11']
            elseif (self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
               self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','phi','psi','Jx ','Jy ','Jz ','rho']
               q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_ph','res_ps','res_Jx','res_Jy', &
                              'res_Jz','res_rh']
               q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09','div10','div11','div12']
            endif
         case default
            self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ','rho']
            q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_Jx','res_Jy','res_Jz','res_rh']
            q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09','div10']
         endselect
         q4_R8P_name = ['curlD_x','curlD_y','curlD_z','curlB_x','curlB_y','curlB_z','curlJ_x','curlJ_y','curlJ_z']
         if (self%io%save_residual_fields) call self%adam%io%register_aux_field(q1_R8P=self%dq,q1_R8P_name=q1_R8P_name)
         if (self%io%save_divergence_fields) call self%adam%io%register_aux_field(q2_R8P=self%divergence,q2_R8P_name=q2_R8P_name)
         if (self%coil%total_coils_number>0) then
            q3_R8P_name = ['j_vec_1','j_vec_2','j_vec_3','f_Gauss']
            call self%adam%io%register_aux_field(q3_R8P=self%coil%j_vec(:,:,:,:,:,1),q3_R8P_name=q3_R8P_name)
            call self%adam%io%register_aux_field(s1_I4P=self%coil%coil_flag,s1_I4P_name='coil_flag')
         endif
         if (self%io%save_curl_fields) call self%adam%io%register_aux_field(q4_R8P=self%curl,q4_R8P_name=q4_R8P_name)
      endif
      endsubroutine io_initialize
   endsubroutine initialize_common

   ! IO methods
   subroutine load_restart_files(self, t, time)
   !< Save restart files.
   class(prism_common_object), intent(inout) :: self !< The equation.
   integer(I4P),               intent(out)   :: t    !< Time iteration.
   real(R8P),                  intent(out)   :: time !< Time.

   call self%adam%load_restart_files(basename=self%io%restart_basename, t=t, time=time, q=self%q)
   call self%adam%make_comm_local_maps_ghost_bc
   endsubroutine load_restart_files

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

   subroutine save_divergence_history(self, is_to_open, is_to_close)
   !< Save divergence history.
   class(prism_common_object), intent(inout)        :: self        !< The equation.
   logical,                    intent(in), optional :: is_to_open  !< Flag to open  file before first saving.
   logical,                    intent(in), optional :: is_to_close !< Flag to close file after last saving.
   real(R8P)                                        :: max_div_D    !< Maximum of divergence of D field.
   real(R8P)                                        :: max_div_B    !< Maximum of divergence of B
   real(R8P)                                        :: max_div_J    !< Maximum of divergence of J field.

   if (self%time%is_to_save(it_save=self%io%divergence_history_save)) then
      max_div_D = maxval(abs(self%divergence(1,:,:,:,:)))
      max_div_B = maxval(abs(self%divergence(2,:,:,:,:)))
      max_div_J = maxval(abs(self%divergence(3,:,:,:,:)))
      call self%io%save_divergence_history(it=self%time%it,time=self%time%time,blocks_number=self%blocks_number, &
                                           div_D=max_div_D,div_B=max_div_B,div_J=max_div_J, &
                                           is_to_open=is_to_open,is_to_close=is_to_close)
   endif
   endsubroutine save_divergence_history

   subroutine save_restart_files(self)
   !< Save restart files.
   class(prism_common_object), intent(inout) :: self !< The equation.

   call self%mpih%barrier(tictoc=.true.)
   call self%mpih%print_message('save restart files t: '//trim(str(self%time%it,.true.))//', time: '//&
                                    trim(str(self%time%time,.true.)))
   call self%adam%save_restart_files(basename=self%io%restart_basename, t=self%time%it, time=self%time%time, q=self%q)
   call self%save_xh5f(output_basename=self%io%restart_basename)
   call self%mpih%barrier(tictoc=.true.)
   endsubroutine save_restart_files

   subroutine save_xh5f(self, output_basename, with_ghost)
   !< Save simulation data in HDF5 format.
   class(prism_common_object), intent(inout)        :: self             !< The equation.
   character(*),               intent(in), optional :: output_basename  !< Output basename.
   logical,                    intent(in), optional :: with_ghost       !< Flag to save ghost cells.
   character(:), allocatable                        :: output_basename_ !< Output basename, local var.

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

   ! FDV operators numerical methods
   subroutine compute_curl_fd(self, ivar, q, curl)
   !< Compute curl of vector fields, div(q(ivar:ivar+2), using finite difference schemes.
   class(prism_common_object), intent(in)    :: self                                            !< The equation.
   integer(I4P),               intent(in)    :: ivar                                            !< Start index of variable of q.
   real(R8P),                  intent(in)    :: q(   1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),                  intent(inout) :: curl(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Curl.
   integer(I4P)                              :: i,j,k,b                                         !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz=>self%field%dxyz, &
             hs=>self%numerics%fdv_half_stencils(1))
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      call compute_curl_fd_centered(s=hs,dxyz=dxyz(:,b),q=q(ivar:ivar+2,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,b),curl=curl(ivar:,i,j,k,b))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_curl_fd

   subroutine compute_curl_fv(self, ivar, q, curl)
   !< Compute curl of vector fields, div(q(ivar:ivar+2), using finite volume schemes.
   class(prism_common_object), intent(in)    :: self                                            !< The equation.
   integer(I4P),               intent(in)    :: ivar                                            !< Start index of variable of q.
   real(R8P),                  intent(in)    :: q(   1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),                  intent(inout) :: curl(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Curl.
   integer(I4P)                              :: i,j,k,b                                         !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz=>self%field%dxyz, &
             hs=>self%numerics%fdv_half_stencils(1))
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      call compute_curl_fv_centered(s=hs,dxyz=dxyz(:,b),q=q(ivar:ivar+2,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,b),curl=curl(ivar:,i,j,k,b))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_curl_fv

   subroutine compute_derivative1_fd(self, dir, ivar, q, dq_ds)
   !< Compute derivative1 of scalar fields, dq(ivar)/ds, using finite difference schemes.
   class(prism_common_object), intent(in)    :: self                                         !< The equation.
   integer(I4P),               intent(in)    :: dir                                          !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),               intent(in)    :: ivar                                         !< Start index of (vec.) variable of q.
   real(R8P),                  intent(in)    :: q(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Field variables.
   real(R8P),                  intent(inout) :: dq_ds(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Derivative1, dq/ds.
   integer(I4P)                              :: i,j,k,b                                      !< Counter.
   integer(I4P)                              :: is,js,ks                                     !< Stencils.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz=>self%field%dxyz, &
             hs=>self%numerics%fdv_half_stencils(1))
   select case(dir)
   case(1)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative1_fd_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i-hs:i+hs,j,k,b),dq_ds=dq_ds(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(2)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative1_fd_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j-hs:j+hs,k,b),dq_ds=dq_ds(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(3)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative1_fd_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j,k-hs:k+hs,b),dq_ds=dq_ds(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   endselect
   endassociate
   endsubroutine compute_derivative1_fd

   subroutine compute_derivative1_fv(self, dir, ivar, q, dq_ds)
   !< Compute derivative1 of scalar fields, dq(ivar)/ds, using finite volume schemes.
   class(prism_common_object), intent(in)    :: self                                         !< The equation.
   integer(I4P),               intent(in)    :: dir                                          !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),               intent(in)    :: ivar                                         !< Start index of (vec.) variable of q.
   real(R8P),                  intent(in)    :: q(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Field variables.
   real(R8P),                  intent(inout) :: dq_ds(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Derivative1, dq/ds.
   integer(I4P)                              :: i,j,k,b                                      !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz=>self%field%dxyz, &
             hs=>self%numerics%fdv_half_stencils(1))
   select case(dir)
   case(1)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative1_fv_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i-hs:i+hs,j,k,b),dq_ds=dq_ds(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(2)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative1_fv_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j-hs:j+hs,k,b),dq_ds=dq_ds(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(3)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative1_fv_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j,k-hs:k+hs,b),dq_ds=dq_ds(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   endselect
   endassociate
   endsubroutine compute_derivative1_fv

   subroutine compute_derivative2_fd(self, dir, ivar, q, d2q_ds2)
   !< Compute derivative2 of scalar fields, d2q(ivar)/ds2, using finite difference schemes.
   class(prism_common_object), intent(in)    :: self                                           !< The equation.
   integer(I4P),               intent(in)    :: dir                                            !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),               intent(in)    :: ivar                                           !< Start index of (vec.) variable of q.
   real(R8P),                  intent(in)    :: q(1:,   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Field variables.
   real(R8P),                  intent(inout) :: d2q_ds2(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Derivative2, d2q/ds2.
   integer(I4P)                              :: i,j,k,b                                        !< Counter.
   integer(I4P)                              :: is,js,ks                                       !< Stencils.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz=>self%field%dxyz, &
             hs=>self%numerics%fdv_half_stencils(2))
   select case(dir)
   case(1)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative2_fd_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i-hs:i+hs,j,k,b),d2q_ds2=d2q_ds2(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(2)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative2_fd_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j-hs:j+hs,k,b),d2q_ds2=d2q_ds2(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(3)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative2_fd_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j,k-hs:k+hs,b),d2q_ds2=d2q_ds2(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   endselect
   endassociate
   endsubroutine compute_derivative2_fd

   subroutine compute_derivative2_fv(self, dir, ivar, q, d2q_ds2)
   !< Compute derivative2 of scalar fields, d2q(ivar)/ds2, using finite volume schemes.
   class(prism_common_object), intent(in)    :: self                                            !< The equation.
   integer(I4P),               intent(in)    :: dir                                             !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),               intent(in)    :: ivar                                            !< Start index of variable of q.
   real(R8P),                  intent(in)    :: q(1:,   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),                  intent(inout) :: d2q_ds2(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative2, d2q/ds2.
   integer(I4P)                              :: i,j,k,b                                         !< Counter.
   integer(I4P)                              :: is,js,ks                                        !< Stencils.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz=>self%field%dxyz, &
             hs=>self%numerics%fdv_half_stencils(2))
   select case(dir)
   case(1)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         ! call compute_derivative2_fv_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i-hs:i+hs,j,k,b),d2q_ds2=d2q_ds2(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(2)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         ! call compute_derivative2_fv_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j-hs:j+hs,k,b),d2q_ds2=d2q_ds2(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(3)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         ! call compute_derivative2_fv_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j,k-hs:k+hs,b),d2q_ds2=d2q_ds2(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   endselect
   endassociate
   endsubroutine compute_derivative2_fv

   subroutine compute_derivative4_fd(self, dir, ivar, q, d4q_ds4)
   !< Compute derivative2 of scalar fields, d4q(ivar)/ds4, using finite difference schemes.
   class(prism_common_object), intent(in)    :: self                                            !< The equation.
   integer(I4P),               intent(in)    :: dir                                             !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),               intent(in)    :: ivar                                            !< Start index of variable of q.
   real(R8P),                  intent(in)    :: q(1:,   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),                  intent(inout) :: d4q_ds4(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative4, d4q/ds4.
   integer(I4P)                              :: i,j,k,b                                         !< Counter.
   integer(I4P)                              :: is,js,ks                                        !< Stencils.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz=>self%field%dxyz, &
             hs=>self%numerics%fdv_half_stencil)
   select case(dir)
   case(1)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative4_fd_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i-hs:i+hs,j,k,b),d4q_ds4=d4q_ds4(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(2)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative4_fd_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j-hs:j+hs,k,b),d4q_ds4=d4q_ds4(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   case(3)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative4_fd_centered(s=hs,ds=dxyz(dir,b),q=q(ivar,i,j,k-hs:k+hs,b),d4q_ds4=d4q_ds4(i,j,k,b))
      enddo
      enddo
      enddo
      enddo
   endselect
   endassociate
   endsubroutine compute_derivative4_fd

   subroutine compute_divergence_fd(self, ivar, q, divergence)
   !< Compute divergence of vector fields, div(q(ivar:ivar+2), using finite difference schemes.
   class(prism_common_object), intent(in)    :: self                                               !< The equation.
   integer(I4P),               intent(in)    :: ivar                                               !< Start index of field of q.
   real(R8P),                  intent(in)    :: q(1:,      1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),                  intent(inout) :: divergence(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
   integer(I4P)                              :: i,j,k,b                                            !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz=>self%field%dxyz, &
             hs=>self%numerics%fdv_half_stencils(1))
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      call compute_divergence_fd_centered(s=hs,dxyz=dxyz(:,b),q=q(ivar:ivar+2,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,b),&
                                          divergence=divergence(i,j,k,b))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_divergence_fd

   subroutine compute_divergence_fv(self, ivar, q, divergence)
   !< Compute divergence of vector fields, div(q(ivar:ivar+2), using finite volume schemes.
   class(prism_common_object), intent(in)    :: self                                               !< The equation.
   integer(I4P),               intent(in)    :: ivar                                               !< Start index of field of q.
   real(R8P),                  intent(in)    :: q(1:,      1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),                  intent(inout) :: divergence(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
   integer(I4P)                              :: i,j,k,b                                            !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz=>self%field%dxyz, &
             hs=>self%numerics%fdv_half_stencils(1))
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      call compute_divergence_fv_centered(s=hs,dxyz=dxyz(:,b),q=q(ivar:ivar+2,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,b),&
                                          divergence=divergence(i,j,k,b))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_divergence_fv

   subroutine compute_gradient_fd(self, ivar, q, gradient)
   !< Compute gradient of scalar variable q(ivar), finite difference schemes.
   class(prism_common_object), intent(in)    :: self                                               !< The equation.
   integer(I4P),               intent(in)    :: ivar                                               !< Index of scalar variable of q.
   real(R8P),                  intent(in)    :: q(       1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Field variables.
   real(R8P),                  intent(inout) :: gradient(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Gradient.
   integer(I4P)                              :: i, j, k, b                                         !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz=>self%field%dxyz, &
             hs=>self%numerics%fdv_half_stencils(1))
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      call compute_gradient_fd_centered(s=hs,dxyz=dxyz(:,b),q=q(ivar,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,b),&
                                        gradient=gradient(1:3,i,j,k,b))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_gradient_fd

   subroutine compute_gradient_fv(self, ivar, q, gradient)
   !< Compute gradient of scalar variable q(ivar), finite volume schemes.
   class(prism_common_object), intent(in)    :: self                                               !< The equation.
   integer(I4P),               intent(in)    :: ivar                                               !< Index of scalar variable of q.
   real(R8P),                  intent(in)    :: q(       1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Field variables.
   real(R8P),                  intent(inout) :: gradient(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)!< Gradient.
   integer(I4P)                              :: i, j, k, b                                         !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz=>self%field%dxyz, &
             hs=>self%numerics%fdv_half_stencils(1))
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      call compute_gradient_fv_centered(s=hs,dxyz=dxyz(:,b),q=q(ivar,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,b),&
                                        gradient=gradient(1:3,i,j,k,b))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_gradient_fv

   subroutine compute_laplacian_fd(self, ivar, q, laplacian)
   !< Compute laplacian of scalar variable q(ivar), finite difference schemes.
   class(prism_common_object), intent(in)    :: self                                              !< The equation.
   integer(I4P),               intent(in)    :: ivar                                              !< Index of scalar variable of q.
   real(R8P),                  intent(in)    :: q(     1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),                  intent(inout) :: laplacian(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Gradient.
   integer(I4P)                              :: i, j, k, b                                        !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz=>self%field%dxyz, &
             hs=>self%numerics%fdv_half_stencils(2))
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      call compute_laplacian_fd_centered(s=hs,dxyz=dxyz(:,b),q=q(ivar,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,b),laplacian=laplacian(i,j,k,b))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_laplacian_fd

   subroutine compute_laplacian_fv(self, ivar, q, laplacian)
   !< Compute laplacian of scalar variable q(ivar), finite volume schemes.
   class(prism_common_object), intent(in)    :: self                                              !< The equation.
   integer(I4P),               intent(in)    :: ivar                                              !< Index of scalar variable of q.
   real(R8P),                  intent(in)    :: q(     1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),                  intent(inout) :: laplacian(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Gradient.
   integer(I4P)                              :: i, j, k, b                                        !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz=>self%field%dxyz, &
             hs=>self%numerics%fdv_half_stencils(2))
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
     !call compute_laplacian_fv_centered(s=hs,dxyz=dxyz(:,b),q=q(ivar,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,b),laplacian=laplacian(i,j,k,b))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_laplacian_fv

   ! coils initialization methods
   subroutine compute_coil_current_density_flux(self, n, adjust_amplitude)
   !< Subroutine to adjust current amplitude in order to match the input one
   class(prism_common_object), intent(inout) :: self             !< Cpu object.
   integer(I4P),               intent(in)    :: n                !< Coil number.
   logical,                    intent(in)    :: adjust_amplitude !< If true, adjust amplitude
   real(R8P)                                 :: x_s, y_s, z_s    !< Flux center coordinates
   real(R8P)                                 :: flux, correction !< Computed flux
   integer(I4P)                              :: i_s, j_s, k_s    !< Flux center cell coordinates
   integer(I4P)                              :: i, j, k          !< Counter.

   associate(blocks_number=>self%blocks_number, ni=>self%ni, nj=>self%nj,                     &
             nk=>self%field%grid%nk, ngc=>self%field%grid%ngc, x_c=>self%coil%x_center(n),    &
             lx=>self%coil%lx(n), ly=>self%coil%ly(n), coil_flag =>self%coil%coil_flag,       &
             y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), dx=>self%field%dxyz(1,:),&
             dy=>self%field%dxyz(2,:), dz=>self%field%dxyz(3,:), normal=>self%coil%normal(n), &
             nb=>self%field%nb, x_cell=>self%field%x_cell, y_cell=>self%field%y_cell,         &
             z_cell=>self%field%z_cell, sigma=>self%coil%sigma(n), J_vec=>self%coil%J_vec,    &
             e_min => self%field%grid%domain_emin, e_max => self%field%grid%domain_emax,      &
             q=>self%q)

   !Per ora la imposto per griglia uniforme monoblocco. Vedremo come estendere il problema
   flux = 0.0_R8P
   correction = 1.0_R8P

   if (normal == NORMAL_P_X .or. normal == NORMAL_M_X) then !Valuto corrente lungo z

      x_s = x_c
      y_s = y_c - lx/2
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

      x_s = x_c - lx/2
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

      x_s = x_c - lx/2
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
   print *, 'Valore corrente pre correzione: ', flux
   if (adjust_amplitude) then
      print*, self%coil%A(n), 'Ampiezza A(n) pre correzione'
      correction = (self%coil%A(n)/flux)
      print *, 'Scaling factor ampiezza: ', correction
      self%coil%A(n) = self%coil%A(n)*correction
      print*, self%coil%A(n), 'Ampiezza A(n) post correzione'
   else
      print *, 'Ampiezza A(n) non corretta: ', self%coil%A(n)
      print *, 'Valore corrente: ', flux
   endif

   endassociate
   endsubroutine compute_coil_current_density_flux

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
      endselect
      call self%compute_divergence(ivar=1_I4P,q=self%coil%J_vec(1:3,:,:,:,:,n),divergence=self%divergence(3,:,:,:,:))
      print *, 'Divergenza J vec della spira: ',n, ' pari a: ',maxval(abs(self%divergence(3,:,:,:,:)))
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
   associate(blocks_number=>self%field%blocks_number, ni=>self%field%grid%ni, nj=>self%field%grid%nj, &
            nk=>self%field%grid%nk, ngc=>self%field%grid%ngc, x_c=>self%coil%x_center(n),             &
            lx=>self%coil%lx(n), ly=>self%coil%ly(n), coil_flag =>self%coil%coil_flag,                &
            y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), dx=>self%field%dxyz(1,:),         &
            dy=>self%field%dxyz(2,:), dz=>self%field%dxyz(3,:), normal=>self%coil%normal(n),          &
            nb=>self%field%nb, x_cell=>self%field%x_cell, y_cell=>self%field%y_cell,                  &
            z_cell=>self%field%z_cell, sigma=>self%coil%sigma(n), J_vec=>self%coil%J_vec)

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

   call self%compute_curl(ivar=1_I4P,q=A,curl=J_vec_buffer)
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
   associate(blocks_number=>self%field%blocks_number, ni=>self%field%grid%ni, nj=>self%field%grid%nj, &
            nk=>self%field%grid%nk, ngc=>self%field%grid%ngc, x_c=>self%coil%x_center(n),             &
            lx=>self%coil%lx(n), ly=>self%coil%ly(n), coil_flag =>self%coil%coil_flag,                &
            y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), dx=>self%field%dxyz(1,:),         &
            dy=>self%field%dxyz(2,:), dz=>self%field%dxyz(3,:), normal=>self%coil%normal(n),          &
            nb=>self%field%nb, x_cell=>self%field%x_cell, y_cell=>self%field%y_cell,                  &
            z_cell=>self%field%z_cell, sigma=>self%coil%sigma(n), J_vec=>self%coil%J_vec)

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

   call self%compute_curl(ivar=1_I4P,q=A,curl=J_vec_buffer)
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
   associate(blocks_number=>self%field%blocks_number, ni=>self%field%grid%ni, nj=>self%field%grid%nj, &
            nk=>self%field%grid%nk, ngc=>self%field%grid%ngc, x_c=>self%coil%x_center(n),             &
            lx=>self%coil%lx(n), ly=>self%coil%ly(n), coil_flag =>self%coil%coil_flag,                &
            y_c=>self%coil%y_center(n), z_c=>self%coil%z_center(n), dx=>self%field%dxyz(1,:),         &
            dy=>self%field%dxyz(2,:), dz=>self%field%dxyz(3,:), normal=>self%coil%normal(n),          &
            nb=>self%field%nb, x_cell=>self%field%x_cell, y_cell=>self%field%y_cell,                  &
            z_cell=>self%field%z_cell, sigma=>self%coil%sigma(n), J_vec=>self%coil%J_vec)

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
   call self%compute_curl(ivar=1_I4P,q=A,curl=J_vec_buffer)
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
