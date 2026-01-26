!< ADAM, PRISM (Plasma Research usIng Simulation Methods) equations system class definition, GPU (FNL) backend.

#include "fundal.H"

module adam_prism_fnl_object
!< ADAM, PRISM (Plasma Research usIng Simulation Methods) equations system class definition, GPU (FNL) backend.

! ADAM modules
use :: adam_common_library
! PRSIM modules
use :: adam_prism_common_library
use :: adam_prism_fnl_library
! third party modules
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
   type(ib_fnl_object)    :: ib_gpu    !< IB handler, FNL backend.
   type(rk_fnl_object)    :: rk_gpu    !< RK integrator, FNL backend.
   type(weno_fnl_object)  :: weno_gpu  !< WENO reconstructor, FNL backend.
   ! PRISM library objects
   type(prism_fnl_coil_object)    :: coil_gpu    !< Coil handler.
   type(prism_fnl_fwlayer_object) :: fwlayer_gpu !< fWLayer handler.
   ! device data
   real(R8P), pointer :: q_gpu(:,:,:,:,:)=>null()           !< Field cell centered variables.
   real(R8P), pointer :: dq_gpu(:,:,:,:,:)=>null()          !< Residuals right hand side.
   real(R8P), pointer :: flxyz_c_gpu(:,:,:,:,:,:,:)=>null() !< Fluxes at cell center with +/- decomposition for all directions.
   real(R8P), pointer :: flx_f_gpu(:,:,:,:,:)=>null()       !< Fluxes along x at cell face.
   real(R8P), pointer :: fly_f_gpu(:,:,:,:,:)=>null()       !< Fluxes along y at cell face.
   real(R8P), pointer :: flz_f_gpu(:,:,:,:,:)=>null()       !< Fluxes along z at cell face.
   real(R8P), pointer :: curl_gpu(:,:,:,:,:)=>null()        !< Curl fields.
   real(R8P), pointer :: divergence_gpu(:,:,:,:,:)=>null()  !< Divergence fields.
   !< Pointer (abstract) TBP.
   procedure(compute_curl_interface),       pass(self),pointer :: compute_curl       =>null()!< Compute curl of vector field.
   procedure(compute_derivative1_interface),pass(self),pointer :: compute_derivative1=>null()!< Compute derivative1 of scalar field.
   procedure(compute_derivative2_interface),pass(self),pointer :: compute_derivative2=>null()!< Compute derivative2 of scalar field.
   procedure(compute_derivative4_interface),pass(self),pointer :: compute_derivative4=>null()!< Compute derivative4 of scalar field.
   procedure(compute_divergence_interface), pass(self),pointer :: compute_divergence =>null()!< Compute divergence of vector field.
   procedure(compute_gradient_interface),   pass(self),pointer :: compute_gradient   =>null()!< Compute gradient of scalar field.
   procedure(compute_laplacian_interface),  pass(self),pointer :: compute_laplacian  =>null()!< Compute laplacian of scalar field.
   procedure(compute_residuals_interface),  pass(self),pointer :: compute_residuals  =>null()!< Compute residuals, space operator.
   procedure(integrate_interface),          pass(self),pointer :: integrate          =>null()!< Integrate, time operator.
   contains
      ! auxiliary methods
      procedure, pass(self) :: allocate_gpu !< Allocate GPU data.
      procedure, pass(self) :: copy_cpu_gpu !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu !< Copy data from GPU to CPU.
      procedure, pass(self) :: initialize   !< Initialize the equation.
      ! IO methods
      procedure, pass(self) :: load_restart_files   !< Load restart files.
      procedure, pass(self) :: save_residuals       !< Save residuals history.
      procedure, pass(self) :: save_simulation_data !< Save all simulation data.
      ! IC/BC/sources
      procedure, pass(self) :: apply_fwl_correction    !< Apply fWLayer correction (if present)
      procedure, pass(self) :: compute_coils_current   !< Compute current coils sources.
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions of equation.
      procedure, pass(self) :: set_initial_conditions  !< Set initial conditions of equation.
      procedure, pass(self) :: update_ghost            !< Update ghost cells and set boundary conditions.
      procedure, pass(self) :: update_rk_ghost         !< Update RK stage ghost cells.
      ! numerical methods, FDV operators
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
      ! numerical methods, space operators
      procedure, pass(self) :: compute_residuals_fd_centered !< Compute residuals, centered finite difference schemes.
      ! procedure, pass(self) :: compute_residuals_fv_centered !< Compute residuals, centered finite volume schemes.
      ! procedure, pass(self) :: compute_residuals_weno        !< Compute residuals, WENO schemes.
      ! numerical methods, time operators
      procedure, pass(self) :: integrate_blanesmoan   !< Blanes and Moan scheme.
      procedure, pass(self) :: integrate_cfm          !< Commutator-Free Magnus scheme.
      procedure, pass(self) :: integrate_leapfrog     !< Leapfrog scheme.
      procedure, pass(self) :: integrate_leapfrog_pic !< Leapfrog scheme, PIC version.
      procedure, pass(self) :: integrate_rk_ls        !< RK classical low storage schemes.
      procedure, pass(self) :: integrate_rk_ssp       !< SSP RK schemes.
      procedure, pass(self) :: integrate_rk_yoshida   !< Yoshida schemes.
      ! numerical methods, miscellanea
      procedure, pass(self) :: compute_auxiliary_fields !< Compute auxiliary fields.
      procedure, pass(self) :: compute_dt               !< Compute time step.
      procedure, pass(self) :: compute_energy           !< Compute energy.
      procedure, pass(self) :: compute_energy_error     !< Compute energy error.
      procedure, pass(self) :: impose_ct_correction     !< Impose Constrained Transport correction on q(ivar:ivar+2).
      procedure, pass(self) :: impose_div_free          !< Impose divergence-free property.
      procedure, pass(self) :: simulate                 !< Perform the simulation.
endtype prism_fnl_object

interface
   subroutine compute_curl_interface(self, ivar, q_gpu, curl_gpu)
   !< Compute curl of vector fields, div(q(ivar:ivar+2).
   import :: prism_fnl_object, I4P, R8P
   class(prism_fnl_object), intent(in)    :: self                                                !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(   1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: curl_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Curl.
   endsubroutine compute_curl_interface

   subroutine compute_derivative1_interface(self, dir, ivar, q_gpu, dq_ds_gpu)
   !< Compute derivative1 of scalar fields, dq(ivar)/ds.
   import :: prism_fnl_object, I4P, R8P
   class(prism_fnl_object), intent(in)    :: self                                              !< The equation.
   integer(I4P),            intent(in)    :: dir                                               !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),            intent(in)    :: ivar                                              !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: dq_ds_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative1, dq/ds.
   endsubroutine compute_derivative1_interface

   subroutine compute_derivative2_interface(self, dir, ivar, q_gpu, d2q_ds2_gpu)
   !< Compute derivative2 of scalar fields, d2q(ivar)/ds2.
   import :: prism_fnl_object, I4P, R8P
   class(prism_fnl_object), intent(in)    :: self                                                !< The equation.
   integer(I4P),            intent(in)    :: dir                                                 !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),            intent(in)    :: ivar                                                !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(1:,   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: d2q_ds2_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative2, d2q/ds2.
   endsubroutine compute_derivative2_interface

   subroutine compute_derivative4_interface(self, dir, ivar, q_gpu, d4q_ds4_gpu)
   !< Compute derivative4 of scalar fields, d4q(ivar)/ds4.
   import :: prism_fnl_object, I4P, R8P
   class(prism_fnl_object), intent(in)    :: self                                                !< The equation.
   integer(I4P),            intent(in)    :: dir                                                 !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),            intent(in)    :: ivar                                                !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(1:,   1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: d4q_ds4_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative4, d4q/ds4.
   endsubroutine compute_derivative4_interface

   subroutine compute_divergence_interface(self, ivar, q_gpu, divergence_gpu)
   !< Compute divergence of vector fields, div(q(ivar:ivar+2).
   import :: prism_fnl_object, I4P, R8P
   class(prism_fnl_object), intent(in)    :: self                                                   !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                   !< Start index of field of q.
   real(R8P),               intent(in)    :: q_gpu(1:,      1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: divergence_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
   endsubroutine compute_divergence_interface

   subroutine compute_gradient_interface(self, ivar, q_gpu, gradient_gpu)
   !< Compute gradient of scalar variable q(ivar).
   import :: prism_fnl_object, I4P, R8P
   class(prism_fnl_object), intent(in)    :: self                                                    !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                    !< Index of scalar var of q.
   real(R8P),               intent(in)    :: q_gpu(       1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: gradient_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Gradient.
   endsubroutine compute_gradient_interface

   subroutine compute_laplacian_interface(self, ivar, q_gpu, laplacian_gpu)
   !< Compute laplacian of scalar variable q(ivar).
   import :: prism_fnl_object, I4P, R8P
   class(prism_fnl_object), intent(in)    :: self                                                  !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                  !< Index of scalar variable of q.
   real(R8P),               intent(in)    :: q_gpu(     1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: laplacian_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Laplacian.
   endsubroutine compute_laplacian_interface

   subroutine compute_residuals_interface(self, q_gpu, dq_gpu, s)
   !< Compute residuals of equation, space operator.
   import :: prism_fnl_object, R8P, I4P
   class(prism_fnl_object), intent(inout) :: self   !< The equation.
   real(R8P),               intent(inout) :: q_gpu(1:,     &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1:)  !< Conservative variables.
   real(R8P),               intent(inout) :: dq_gpu(1:,     &
                                                1-self%ngc:,&
                                                1-self%ngc:,&
                                                1-self%ngc:,&
                                                1:) !< Residuals.
   integer(I4P),  optional, intent(in)    :: s      !< Stage counter.
   endsubroutine compute_residuals_interface

   subroutine integrate_interface(self)
   !< Integrate equation, time operator.
   import :: prism_fnl_object, R8P
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   endsubroutine integrate_interface
endinterface

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
   call dev_alloc(fptr_dev=self%flxyz_c_gpu, &
                  ubounds=[nb,3,3,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1,1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   call dev_alloc(fptr_dev=self%flx_f_gpu, &
                  ubounds=[nb,ni,nj,nk,nv], lbounds=[1,0,1,1,1], init_value=0._R8P, ierr=ierr)
   call dev_alloc(fptr_dev=self%fly_f_gpu, &
                  ubounds=[nb,ni,nj,nk,nv], lbounds=[1,1,0,1,1], init_value=0._R8P, ierr=ierr)
   call dev_alloc(fptr_dev=self%flz_f_gpu, &
                  ubounds=[nb,ni,nj,nk,nv], lbounds=[1,1,1,0,1], init_value=0._R8P, ierr=ierr)
   call dev_alloc(fptr_dev=self%curl_gpu, &
                  ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   call dev_alloc(fptr_dev=self%divergence_gpu, &
                  ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   ! call dev_alloc(fptr_dev=self%si_gpu,  ubounds=[3,3], init_value=0_I4P,  ierr=ierr)
   ! call dev_alloc(fptr_dev=self%sir_gpu, ubounds=[3,3], init_value=0._R8P, ierr=ierr)
   ! call dev_assign_to_device(src=ER,   dst=self%ER_GPU  )
   ! call dev_assign_to_device(src=EL,   dst=self%EL_GPU  )
   ! call dev_assign_to_device(src=IERL, dst=self%IERL_GPU)
   endassociate
   call self%mpih%print_message('prism_fnl_object%allocate_gpu finish')
   endsubroutine allocate_gpu

   subroutine copy_cpu_gpu(self)
   !< Copy data from CPU to GPU.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   call self%field_gpu%copy_transpose_cpu_gpu(nv=self%nv, q_cpu=self%q,          q_gpu=self%q_gpu            )
   call self%field_gpu%copy_transpose_cpu_gpu(nv=self%nv, q_cpu=self%curl,       q_gpu=self%curl_gpu         )
   call self%field_gpu%copy_transpose_cpu_gpu(nv=self%nv, q_cpu=self%divergence, q_gpu=self%divergence_gpu   )
   call self%field_gpu%copy_transpose_cpu_gpu(nv=3,       q_cpu=self%fWLayer%f,  q_gpu=self%fwlayer_gpu%f_gpu)
   call self%field_gpu%copy_cpu_gpu(verbose=.false.)
   call self%coil_gpu%copy_cpu_gpu
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self, compute_copy_q_aux, copy_phi)
   !< Copy data from GPU to CPU.
   class(prism_fnl_object), intent(inout)        :: self               !< The equation.
   logical,                 intent(in), optional :: compute_copy_q_aux !< Flag to compute auxiliary variables.
   logical,                 intent(in), optional :: copy_phi           !< Copy also phi.

   call self%field_gpu%copy_transpose_gpu_cpu(nv=self%nv, q_gpu=self%q_gpu,             q_cpu=self%q         )
   call self%field_gpu%copy_transpose_gpu_cpu(nv=self%nv, q_gpu=self%curl_gpu,          q_cpu=self%curl      )
   call self%field_gpu%copy_transpose_gpu_cpu(nv=self%nv, q_gpu=self%divergence_gpu,    q_cpu=self%divergence)
   call self%field_gpu%copy_transpose_gpu_cpu(nv=3,       q_gpu=self%fwlayer_gpu%f_gpu, q_cpu=self%fWLayer%f )
   call self%coil_gpu%copy_gpu_cpu
   endsubroutine copy_gpu_cpu

   subroutine initialize(self, filename)
   !< Initialize the equation.
   class(prism_fnl_object), intent(inout) :: self     !< The equation.
   character(*),            intent(in)    :: filename !< Input file name.

   call self%mpih_gpu%initialize(do_mpi_init=.true., do_device_init=.true., verbose=.true.)
   call self%mpih_gpu%print_message('prism_fnl_object%initialize start')
   call self%initialize_common(field = self%adam%field, filename=filename, memory_avail=real(self%mpih_gpu%dev_memory_avail,R8P), &
                               verbose=.true.)
   call self%field_gpu%initialize(field=self%adam%field, verbose=.true.)
   call self%ib_gpu%initialize(ib=self%ib, field_gpu=self%field_gpu)
   call self%rk_gpu%initialize(rk=self%rk, nb=self%nb, ngc=self%ngc, ni=self%ni, nj=self%nj, nk=self%nk, nv=self%nv)
   call self%weno_gpu%initialize(weno=self%weno)
   call self%allocate_gpu(q_gpu=self%field_gpu%q_gpu)
   call self%coil_gpu%initialize(coil=self%coil,&
                                 blocks_number=self%blocks_number,nb=self%nb,ngc=self%ngc,ni=self%ni,nj=self%nj,nk=self%nk)
   call self%fwlayer_gpu%initialize(field=self%field)
   ! call set_sir_dev(si_gpu=self%si_gpu, sir_gpu=self%sir_gpu)

   ! set pointer (abstract) TBP
   if (self%physics%physical_model == EM_PHYSICAL_MODEL) then
      select case(self%numerics%scheme_time)
      case(NUM_SCHEME_TIME_BLANES_MOAN)
         self%integrate => integrate_blanesmoan
      case(NUM_SCHEME_TIME_CFM)
         self%integrate => integrate_cfm
      case(NUM_SCHEME_TIME_LEAPFROG)
         self%integrate => integrate_leapfrog
      case(NUM_SCHEME_TIME_RUNGE_KUTTA)
         select case(self%rk%scheme)
         case(RK_1, RK_2, RK_3)
            self%integrate => integrate_rk_ls
         case(RK_SSP_22, RK_SSP_33, RK_SSP_54)
            self%integrate => integrate_rk_ssp
         case(RK_YOSHIDA)
            self%integrate => integrate_rk_yoshida
         endselect
      endselect
   elseif (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
      select case(self%numerics%scheme_time)
      case(NUM_SCHEME_TIME_LEAPFROG)
         select case(self%pic%scheme_time)
         case(NUM_SCHEME_TIME_PIC_LEAPFROG)
            self%integrate => integrate_leapfrog_pic
         case(NUM_SCHEME_TIME_PIC_RUNGE_KUTTA)
            !self%integrate =>
         endselect
      case(NUM_SCHEME_TIME_RUNGE_KUTTA)
         select case(self%pic%scheme_time)
         case(NUM_SCHEME_TIME_PIC_LEAPFROG)
            self%integrate => integrate_leapfrog_pic
         case(NUM_SCHEME_TIME_PIC_RUNGE_KUTTA)
            !self%integrate =>
         endselect
      endselect
   endif

   select case(self%numerics%scheme_space)
   case(NUM_SCHEME_SPACE_WENO)
      self%compute_curl        => compute_curl_fv
      self%compute_derivative1 => compute_derivative1_fv
      self%compute_derivative2 => compute_derivative2_fv
      self%compute_divergence  => compute_divergence_fv
      self%compute_gradient    => compute_gradient_fv
      self%compute_laplacian   => compute_laplacian_fv
      ! self%compute_residuals   => compute_residuals_weno
   case(NUM_SCHEME_SPACE_FD_CENTERED)
      self%compute_curl        => compute_curl_fd
      self%compute_derivative1 => compute_derivative1_fd
      self%compute_derivative2 => compute_derivative2_fd
      self%compute_divergence  => compute_divergence_fd
      self%compute_gradient    => compute_gradient_fd
      self%compute_laplacian   => compute_laplacian_fd
      self%compute_residuals   => compute_residuals_fd_centered
   case(NUM_SCHEME_SPACE_FV_CENTERED)
      self%compute_curl        => compute_curl_fv
      self%compute_derivative1 => compute_derivative1_fv
      self%compute_derivative2 => compute_derivative2_fv
      self%compute_divergence  => compute_divergence_fv
      self%compute_gradient    => compute_gradient_fv
      self%compute_laplacian   => compute_laplacian_fv
      ! self%compute_residuals   => compute_residuals_fv_centered
   endselect

   call external_fields_initialize_dev(external_fields=self%external_fields)

   call self%mpih%print_message('prism_fnl_object%initialize finish')
   endsubroutine initialize

   ! IO methods
   subroutine load_restart_files(self, t, time)
   !< Save restart files.
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P),            intent(out)   :: t    !< Time iteration.
   real(R8P),               intent(out)   :: time !< Time.

   call self%prism_common_object%load_restart_files(t=t, time=time)
   call self%copy_cpu_gpu
   endsubroutine load_restart_files

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

   subroutine save_simulation_data(self)
   !< Save all simulation data.
   class(prism_fnl_object), intent(inout) :: self      !< The equation.

   if ((self%time%is_to_save(it_save=self%io%it_save)).or.      &
       (self%time%is_to_save(it_save=self%io%restart_save)).or. &
       (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max))) then
      call self%update_ghost(q_gpu=self%q_gpu)
      call self%copy_gpu_cpu
      call self%compute_auxiliary_fields

      if (self%time%is_to_save(it_save=self%io%it_save)) call self%save_xh5f(with_ghost=.true.)
      if (mod(self%time%it,self%io%restart_save)==0) call self%save_restart_files
      if (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max)) then
         call self%slices%save_mat(basename=self%io%output_basename, &
                                   it=self%time%it,                  &
                                   it_max=self%time%it_max,          &
                                   time=self%time%time,              &
                                   time_max=self%time%time_max,      &
                                   adam=self%adam,                   &
                                   q=self%q,                         &
                                   q_name=self%q_name)

      endif
   endif
   endsubroutine save_simulation_data

   ! IC/BC/sources
   subroutine apply_fwl_correction(self)
   !< Apply correction if a fWL is present.
   class(prism_fnl_object), intent(inout) :: self                    !< The equation.
   integer(I4P)                           :: ni(2,6),nj(2,6),nk(2,6) !< Dimensions of FWL domain.
   real(R8P)                              :: s2(6)                   !< Side coefficient.
   integer(I4P)                           :: n(6)                    !< FWL f function index.
   integer(I4P)                           :: alfa_D(6), beta_D(6)    !< Corrected var index of D (Barbas' notation).
   integer(I4P)                           :: alfa_B(6), beta_B(6)    !< Corrected var index of D (Barbas' notation).
   integer(I4P)                           :: face                    !< Counter.

   associate(C=>self%fWLayer%C, layer=>self%fWLayer%layer)
   if (C>0) then
      ! below arrays should be initialized elsewhere...
        ni(1,1)=1_I4P   ;   ni(1,2)=self%ni-C ;   ni(1,3)=1_I4P   ;   ni(1,4)=1_I4P     ;   ni(1,5)=1_I4P   ;   ni(1,6)=1_I4P
        ni(2,1)=C       ;   ni(2,2)=self%ni   ;   ni(2,3)=self%ni ;   ni(2,4)=self%ni   ;   ni(2,5)=self%ni ;   ni(2,6)=self%ni
        nj(1,1)=1_I4P   ;   nj(1,2)=1_I4P     ;   nj(1,3)=1_I4P   ;   nj(1,4)=self%nj-C ;   nj(1,5)=1_I4P   ;   nj(1,6)=1_I4P
        nj(2,1)=self%nj ;   nj(2,2)=self%nj   ;   nj(2,3)=C       ;   nj(2,4)=self%nj   ;   nj(2,5)=self%nj ;   nj(2,6)=self%nj
        nk(1,1)=1_I4P   ;   nk(1,2)=1_I4P     ;   nk(1,3)=1_I4P   ;   nk(1,4)=1_I4P     ;   nk(1,5)=1_I4P   ;   nk(1,6)=self%nk-C
        nk(2,1)=self%nk ;   nk(2,2)=self%nk   ;   nk(2,3)=self%nk ;   nk(2,4)=self%nk   ;   nk(2,5)=C       ;   nk(2,6)=self%nk
           n(1)=1_I4P   ;      n(2)= 1_I4P    ;      n(3)=2_I4P   ;      n(4)= 2_I4P    ;      n(5)=3_I4P   ;      n(6)= 3_I4P
          s2(1)=1.0_R8P ;     s2(2)=-1.0_R8P  ;     s2(3)=1.0_R8P ;     s2(4)=-1.0_R8P  ;     s2(5)=1.0_R8P ;     s2(6)=-1.0_R8P
      alfa_D(1)=2_I4P   ; alfa_D(2)= 2_I4P    ; alfa_D(3)=3_I4P   ; alfa_D(4)= 3_I4P    ; alfa_D(5)=1_I4P   ; alfa_D(6)= 1_I4P
      beta_D(1)=3_I4P   ; beta_D(2)= 3_I4P    ; beta_D(3)=1_I4P   ; beta_D(4)= 1_I4P    ; beta_D(5)=2_I4P   ; beta_D(6)= 2_I4P
      alfa_B(1)=5_I4P   ; alfa_B(2)= 5_I4P    ; alfa_B(3)=6_I4P   ; alfa_B(4)= 6_I4P    ; alfa_B(5)=4_I4P   ; alfa_B(6)= 4_I4P
      beta_B(1)=6_I4P   ; beta_B(2)= 6_I4P    ; beta_B(3)=4_I4P   ; beta_B(4)= 4_I4P    ; beta_B(5)=5_I4P   ; beta_B(6)= 5_I4P
      do face=1, 6
         if (layer(face)) call apply_fwl_correction_dev(blocks_number=self%blocks_number,ngc=self%ngc,&
                                                        ni1   =  ni(1,face),                          &
                                                        ni2   =  ni(2,face),                          &
                                                        nj1   =  nj(1,face),                          &
                                                        nj2   =  nj(2,face),                          &
                                                        nk1   =  nk(1,face),                          &
                                                        nk2   =  nk(2,face),                          &
                                                        n     =     n(face),                          &
                                                        s2    =    s2(face),                          &
                                                        alfa_D=alfa_D(face),                          &
                                                        beta_D=beta_D(face),                          &
                                                        alfa_B=alfa_B(face),                          &
                                                        beta_B=beta_B(face),                          &
                                                        f_gpu=self%fwlayer_gpu%f_gpu,q_gpu=self%q_gpu)
      enddo
   endif
   endassociate
   endsubroutine apply_fwl_correction

   subroutine compute_coils_current(self, gamma)
   !< Compute current coils sources.
   class(prism_fnl_object), intent(inout)        :: self   !< The equation.
   real(R8P),               intent(in), optional :: gamma  !< RK coefficient.
   real(R8P)                                     :: time_s !< Local time.

   if (self%coil%total_coils_number >= 1_I4P) then
      time_s = self%time%time ; if (present(gamma)) time_s = self%time%time + self%time%dt*gamma
      call compute_coils_current_dev(ni=self%ni, nj=self%nj, nk=self%nk, ngc=self%ngc, blocks_number=self%blocks_number,     &
                                     time_s=time_s, td=self%coil%td,                                                         &
                                     A_gpu=self%coil_gpu%A_gpu, f_gpu=self%coil_gpu%f_gpu, phase_gpu=self%coil_gpu%phase_gpu,&
                                     coil_flag_gpu=self%coil_gpu%coil_flag_gpu, j_vec_gpu=self%coil_gpu%j_vec_gpu,           &
                                     var_Jx=self%physics%var_Jx, var_Jy=self%physics%var_Jy, var_Jz=self%physics%var_Jz,     &
                                     q_gpu=self%q_gpu)
   endif
   endsubroutine compute_coils_current

   subroutine set_boundary_conditions(self, q_gpu)
   !< Set boundary conditions of equation.
   class(prism_fnl_object), intent(in)    :: self                  !< The equation.
   real(R8P),               intent(inout) :: q_gpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,1:) !< Conservative variables.
   integer(I4P)                           :: b                     !< Counter.
   integer(I4P)                           :: c, i, j, k, v         !< Counter.
   integer(I4P)                           :: idelta                !< IJK i delta step for extrapolation.
   integer(I4P)                           :: jdelta                !< IJK j delta step for extrapolation.
   integer(I4P)                           :: kdelta                !< IJK k delta step for extrapolation.
   integer(I4P)                           :: bc_type               !< Boundary condition type.
   integer(I4P)                           :: crown                 !< Crown counter.

   associate(local_map_bc_crown_gpu=>self%field_gpu%maps%local_map_bc_crown_gpu, &
             nv=>self%nv, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk)
   if (associated(self%field_gpu%maps%local_map_bc_crown_gpu)) then
      do crown=1, ngc
         !$acc parallel loop independent gang vector DEVICEVAR(local_map_bc_crown_gpu, q_gpu)
         do c=1, size(local_map_bc_crown_gpu, dim=1)
            b = local_map_bc_crown_gpu(c, 1 ,crown)
            if (b>0) then
               i       = local_map_bc_crown_gpu(c, 2 ,crown)
               j       = local_map_bc_crown_gpu(c, 3 ,crown)
               k       = local_map_bc_crown_gpu(c, 4 ,crown)
               idelta  = local_map_bc_crown_gpu(c, 5 ,crown)
               jdelta  = local_map_bc_crown_gpu(c, 6 ,crown)
               kdelta  = local_map_bc_crown_gpu(c, 7 ,crown)
               bc_type = local_map_bc_crown_gpu(c, 8 ,crown)
               if (bc_type == BC_EXTRAPOLATION) then
                  do v=1, nv
                     q_gpu(b,i,j,k,v) = q_gpu(b,i-idelta,j-jdelta,k-kdelta,v)
                  enddo
               elseif (bc_type == BC_NEUMANN) then
                  do v=1, nv
                     q_gpu(b,i,j,k,v) = q_gpu(b,i+abs(idelta)*(-2*i+1+(idelta+1)*ni),&
                                                j+abs(jdelta)*(-2*j+1+(jdelta+1)*nj),&
                                                k+abs(kdelta)*(-2*k+1+(kdelta+1)*nk),v)
                  enddo
               elseif (bc_type == BC_SILVER_MULLER) then
                  ! to be impelmented
               elseif (bc_type == BC_DIRICHLET) then
                  do v=1, nv
                     q_gpu(b,i,j,k,v) = 0._R8P
                  enddo
               elseif (bc_type == BC_PERIOD) then
                  ! to be impelmented
               endif
            endif
         enddo
      enddo
   endif
   endassociate
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self)
   !< Set initial conditions of field.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   call self%ic%set_initial_conditions(physics=self%physics, field=self%field, q=self%q)
   call self%coil%set_coils(physics=self%physics, field=self%field)
   call self%copy_cpu_gpu
   endsubroutine set_initial_conditions

   subroutine update_ghost(self, q_gpu, step, s)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(prism_fnl_object), intent(inout)        :: self            !< The equation.
   real(R8P),               intent(inout)        :: q_gpu(1:,         &
                                                          1-self%ngc:,&
                                                          1-self%ngc:,&
                                                          1-self%ngc:,&
                                                          1:)       !< Conservative variables.
   integer(I4P),            intent(in), optional :: step            !< Step to be perfordmed in asyncronous comp.
   integer(I4P),            intent(in), optional :: s               !< Stage counter.
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
   if (do_set_bc)       call self%apply_fwl_correction
   if (do_set_bc)       call self%set_boundary_conditions(q_gpu=q_gpu)
   endsubroutine update_ghost

   subroutine update_rk_ghost(self, dt, phi_gpu)
   !< Update RK ghost cells.
   class(prism_fnl_object), intent(inout)        :: self                 !< RK object.
   real(R8P),               intent(in)           :: dt                   !< Current time step.
   real(R8P),               intent(in), optional :: phi_gpu(1:,          &
                                                            1-self%ngc:, &
                                                            1-self%ngc:, &
                                                            1-self%ngc:, &
                                                            1:)          !< IB distance.
   ! to be implemented
   endsubroutine update_rk_ghost

   ! numerical methods, FDV operators
   subroutine compute_curl_fd(self, ivar, q_gpu, curl_gpu)
   !< Compute curl of vector fields, div(q(ivar:ivar+2), using finite difference schemes.
   class(prism_fnl_object), intent(in)    :: self                                                !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(   1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: curl_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Curl.
   integer(I4P)                           :: i,j,k,b                                             !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_gpu%dxyz_gpu,&
             hs=>self%numerics%fdv_half_stencils(1))
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,curl_gpu,dxyz_gpu)
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      call compute_curl_fd_centered(s=hs,dxyz=dxyz_gpu(b,:),q=q_gpu(b,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,ivar:ivar+2),&
                                    curl=curl_gpu(b,i,j,k,ivar:))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_curl_fd

   subroutine compute_curl_fv(self, ivar, q_gpu, curl_gpu)
   !< Compute curl of vector fields, div(q(ivar:ivar+2), using finite volume schemes.
   class(prism_fnl_object), intent(in)    :: self                                                !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(   1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: curl_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Curl.
   integer(I4P)                           :: i,j,k,b                                             !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_gpu%dxyz_gpu,&
             hs=>self%numerics%fdv_half_stencils(1))
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,curl_gpu,dxyz_gpu)
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      call compute_curl_fv_centered(s=hs,dxyz=dxyz_gpu(b,:),q=q_gpu(b,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,ivar:ivar+2),&
                                    curl=curl_gpu(b,i,j,k,ivar:))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_curl_fv

   subroutine compute_derivative1_fd(self, dir, ivar, q_gpu, dq_ds_gpu)
   !< Compute derivative1 of scalar fields, dq(ivar)/ds, using finite difference schemes.
   class(prism_fnl_object), intent(in)    :: self                                              !< The equation.
   integer(I4P),            intent(in)    :: dir                                               !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),            intent(in)    :: ivar                                              !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: dq_ds_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative1, dq/ds.
   integer(I4P)                           :: i,j,k,b                                           !< Counter.
   integer(I4P)                           :: is,js,ks                                          !< Stencils.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_gpu%dxyz_gpu,&
             hs=>self%numerics%fdv_half_stencils(1))
   select case(dir)
   case(1)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,dq_ds_gpu,dxyz_gpu)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative1_fd_centered(s=hs,ds=dxyz_gpu(b,dir),q=q_gpu(b,i-hs:i+hs,j,k,ivar),dq_ds=dq_ds_gpu(b,i,j,k))
      enddo
      enddo
      enddo
      enddo
   case(2)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,dq_ds_gpu,dxyz_gpu)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative1_fd_centered(s=hs,ds=dxyz_gpu(b,dir),q=q_gpu(b,i,j-hs:j+hs,k,ivar),dq_ds=dq_ds_gpu(b,i,j,k))
      enddo
      enddo
      enddo
      enddo
   case(3)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,dq_ds_gpu,dxyz_gpu)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative1_fd_centered(s=hs,ds=dxyz_gpu(b,dir),q=q_gpu(b,i,j,k-hs:k+hs,ivar),dq_ds=dq_ds_gpu(b,i,j,k))
      enddo
      enddo
      enddo
      enddo
   endselect
   endassociate
   endsubroutine compute_derivative1_fd

   subroutine compute_derivative1_fv(self, dir, ivar, q_gpu, dq_ds_gpu)
   !< Compute derivative1 of scalar fields, dq(ivar)/ds, using finite volume schemes.
   class(prism_fnl_object), intent(in)    :: self                                              !< The equation.
   integer(I4P),            intent(in)    :: dir                                               !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),            intent(in)    :: ivar                                              !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: dq_ds_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative1, dq/ds.
   integer(I4P)                           :: i,j,k,b                                           !< Counter.
   integer(I4P)                           :: is,js,ks                                          !< Stencils.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_gpu%dxyz_gpu,&
             hs=>self%numerics%fdv_half_stencils(1))
   select case(dir)
   case(1)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,dq_ds_gpu,dxyz_gpu)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative1_fv_centered(s=hs,ds=dxyz_gpu(b,dir),q=q_gpu(b,i-hs:i+hs,j,k,ivar),dq_ds=dq_ds_gpu(b,i,j,k))
      enddo
      enddo
      enddo
      enddo
   case(2)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,dq_ds_gpu,dxyz_gpu)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative1_fv_centered(s=hs,ds=dxyz_gpu(b,dir),q=q_gpu(b,i,j-hs:j+hs,k,ivar),dq_ds=dq_ds_gpu(b,i,j,k))
      enddo
      enddo
      enddo
      enddo
   case(3)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,dq_ds_gpu,dxyz_gpu)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative1_fv_centered(s=hs,ds=dxyz_gpu(b,dir),q=q_gpu(b,i,j,k-hs:k+hs,ivar),dq_ds=dq_ds_gpu(b,i,j,k))
      enddo
      enddo
      enddo
      enddo
   endselect
   endassociate
   endsubroutine compute_derivative1_fv

   subroutine compute_derivative2_fd(self, dir, ivar, q_gpu, d2q_ds2_gpu)
   !< Compute derivative2 of scalar fields, d2q(ivar)/ds2, using finite difference schemes.
   class(prism_fnl_object), intent(in)    :: self                                                !< The equation.
   integer(I4P),            intent(in)    :: dir                                                 !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),            intent(in)    :: ivar                                                !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)   !< Field variables.
   real(R8P),               intent(inout) :: d2q_ds2_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative2, d2q/ds2.
   integer(I4P)                           :: i,j,k,b                                             !< Counter.
   integer(I4P)                           :: is,js,ks                                            !< Stencils.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_gpu%dxyz_gpu,&
             hs=>self%numerics%fdv_half_stencils(1))
   select case(dir)
   case(1)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,d2q_ds2_gpu,dxyz_gpu)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative2_fd_centered(s=hs,ds=dxyz_gpu(b,dir),q=q_gpu(b,i-hs:i+hs,j,k,ivar),d2q_ds2=d2q_ds2_gpu(b,i,j,k))
      enddo
      enddo
      enddo
      enddo
   case(2)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,d2q_ds2_gpu,dxyz_gpu)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative2_fd_centered(s=hs,ds=dxyz_gpu(b,dir),q=q_gpu(b,i,j-hs:j+hs,k,ivar),d2q_ds2=d2q_ds2_gpu(b,i,j,k))
      enddo
      enddo
      enddo
      enddo
   case(3)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,d2q_ds2_gpu,dxyz_gpu)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative2_fd_centered(s=hs,ds=dxyz_gpu(b,dir),q=q_gpu(b,i,j,k-hs:k+hs,ivar),d2q_ds2=d2q_ds2_gpu(b,i,j,k))
      enddo
      enddo
      enddo
      enddo
   endselect
   endassociate
   endsubroutine compute_derivative2_fd

   subroutine compute_derivative2_fv(self, dir, ivar, q_gpu, d2q_ds2_gpu)
   !< Compute derivative2 of scalar fields, d2q(ivar)/ds2, using finite volume schemes.
   class(prism_fnl_object), intent(in)    :: self                                                !< The equation.
   integer(I4P),            intent(in)    :: dir                                                 !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),            intent(in)    :: ivar                                                !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)   !< Field variables.
   real(R8P),               intent(inout) :: d2q_ds2_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative2, d2q/ds2.
   integer(I4P)                           :: i,j,k,b                                             !< Counter.
   integer(I4P)                           :: is,js,ks                                            !< Stencils.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_gpu%dxyz_gpu,&
             hs=>self%numerics%fdv_half_stencils(1))
   select case(dir)
   case(1)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,d2q_ds2_gpu,dxyz_gpu)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative2_fv_centered(s=hs,ds=dxyz_gpu(b,dir),q=q_gpu(b,i-hs:i+hs,j,k,ivar),d2q_ds2=d2q_ds2_gpu(b,i,j,k))
      enddo
      enddo
      enddo
      enddo
   case(2)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,d2q_ds2_gpu,dxyz_gpu)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative2_fv_centered(s=hs,ds=dxyz_gpu(b,dir),q=q_gpu(b,i,j-hs:j+hs,k,ivar),d2q_ds2=d2q_ds2_gpu(b,i,j,k))
      enddo
      enddo
      enddo
      enddo
   case(3)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,d2q_ds2_gpu,dxyz_gpu)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative2_fv_centered(s=hs,ds=dxyz_gpu(b,dir),q=q_gpu(b,i,j,k-hs:k+hs,ivar),d2q_ds2=d2q_ds2_gpu(b,i,j,k))
      enddo
      enddo
      enddo
      enddo
   endselect
   endassociate
   endsubroutine compute_derivative2_fv

   subroutine compute_derivative4_fd(self, dir, ivar, q_gpu, d4q_ds4_gpu)
   !< Compute derivative4 of scalar fields, d4q(ivar)/ds4, using finite difference schemes.
   class(prism_fnl_object), intent(in)    :: self                                                !< The equation.
   integer(I4P),            intent(in)    :: dir                                                 !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P),            intent(in)    :: ivar                                                !< Start index of variable of q.
   real(R8P),               intent(in)    :: q_gpu(1:, 1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)   !< Field variables.
   real(R8P),               intent(inout) :: d4q_ds4_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Derivative4, d4q/ds4.
   integer(I4P)                           :: i,j,k,b                                             !< Counter.
   integer(I4P)                           :: is,js,ks                                            !< Stencils.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_gpu%dxyz_gpu,&
             hs=>self%numerics%fdv_half_stencils(1))
   select case(dir)
   case(1)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,d4q_ds4_gpu,dxyz_gpu)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative4_fd_centered(s=hs,ds=dxyz_gpu(b,dir),q=q_gpu(b,i-hs:i+hs,j,k,ivar),d4q_ds4=d4q_ds4_gpu(b,i,j,k))
      enddo
      enddo
      enddo
      enddo
   case(2)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,d4q_ds4_gpu,dxyz_gpu)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative4_fd_centered(s=hs,ds=dxyz_gpu(b,dir),q=q_gpu(b,i,j-hs:j+hs,k,ivar),d4q_ds4=d4q_ds4_gpu(b,i,j,k))
      enddo
      enddo
      enddo
      enddo
   case(3)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,d4q_ds4_gpu,dxyz_gpu)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         call compute_derivative4_fd_centered(s=hs,ds=dxyz_gpu(b,dir),q=q_gpu(b,i,j,k-hs:k+hs,ivar),d4q_ds4=d4q_ds4_gpu(b,i,j,k))
      enddo
      enddo
      enddo
      enddo
   endselect
   endassociate
   endsubroutine compute_derivative4_fd

   subroutine compute_divergence_fd(self, ivar, q_gpu, divergence_gpu)
   !< Compute divergence of vector fields, div(q(ivar:ivar+2), using finite difference schemes.
   class(prism_fnl_object), intent(in)    :: self                                                   !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                   !< Start index of field of q.
   real(R8P),               intent(in)    :: q_gpu(1:,      1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: divergence_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
   integer(I4P)                           :: i,j,k,b                                                !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_gpu%dxyz_gpu,&
             hs=>self%numerics%fdv_half_stencils(1))
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,divergence_gpu,dxyz_gpu)
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      call compute_divergence_fd_centered(s=hs,dxyz=dxyz_gpu(b,:),q=q_gpu(b,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,ivar:ivar+2),&
                                          divergence=divergence_gpu(b,i,j,k))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_divergence_fd

   subroutine compute_divergence_fv(self, ivar, q_gpu, divergence_gpu)
   !< Compute divergence of vector fields, div(q(ivar:ivar+2), using finite volume schemes.
   class(prism_fnl_object), intent(in)    :: self                                                   !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                   !< Start index of field of q.
   real(R8P),               intent(in)    :: q_gpu(1:,      1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: divergence_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Divergence.
   integer(I4P)                           :: i,j,k,b                                                !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_gpu%dxyz_gpu,&
             hs=>self%numerics%fdv_half_stencils(1))
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,divergence_gpu,dxyz_gpu)
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      call compute_divergence_fv_centered(s=hs,dxyz=dxyz_gpu(b,:),q=q_gpu(b,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,ivar:ivar+2),&
                                          divergence=divergence_gpu(b,i,j,k))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_divergence_fv

   subroutine compute_gradient_fd(self, ivar, q_gpu, gradient_gpu)
   !< Compute gradient of scalar variable q(ivar), finite difference schemes.
   class(prism_fnl_object), intent(in)    :: self                                                    !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                    !< Index of scalar var of q.
   real(R8P),               intent(in)    :: q_gpu(       1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: gradient_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Gradient.
   integer(I4P)                           :: i, j, k, b                                              !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_gpu%dxyz_gpu,&
             hs=>self%numerics%fdv_half_stencils(1))
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,gradient_gpu,dxyz_gpu)
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      call compute_gradient_fd_centered(s=hs,dxyz=dxyz_gpu(b,:),q=q_gpu(b,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,ivar),&
                                        gradient=gradient_gpu(b,i,j,k,1:3))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_gradient_fd

   subroutine compute_gradient_fv(self, ivar, q_gpu, gradient_gpu)
   !< Compute gradient of scalar variable q(ivar), finite volume schemes.
   class(prism_fnl_object), intent(in)    :: self                                                    !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                    !< Index of scalar var of q.
   real(R8P),               intent(in)    :: q_gpu(       1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: gradient_gpu(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Gradient.
   integer(I4P)                           :: i, j, k, b                                              !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_gpu%dxyz_gpu,&
             hs=>self%numerics%fdv_half_stencils(1))
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,gradient_gpu,dxyz_gpu)
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      call compute_gradient_fv_centered(s=hs,dxyz=dxyz_gpu(b,:),q=q_gpu(b,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,ivar),&
                                        gradient=gradient_gpu(b,i,j,k,1:3))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_gradient_fv

   subroutine compute_laplacian_fd(self, ivar, q_gpu, laplacian_gpu)
   !< Compute laplacian of scalar variable q(ivar), finite difference schemes.
   class(prism_fnl_object), intent(in)    :: self                                                  !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                  !< Index of scalar variable of q.
   real(R8P),               intent(in)    :: q_gpu(     1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: laplacian_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Gradient.
   integer(I4P)                           :: i, j, k, b                                            !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_gpu%dxyz_gpu,&
             hs=>self%numerics%fdv_half_stencils(2))
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,laplacian_gpu,dxyz_gpu)
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      call compute_laplacian_fd_centered(s=hs,dxyz=dxyz_gpu(b,:),q=q_gpu(b,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,ivar),&
                                         laplacian=laplacian_gpu(b,i,j,k))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_laplacian_fd

   subroutine compute_laplacian_fv(self, ivar, q_gpu, laplacian_gpu)
   !< Compute laplacian of scalar variable q(ivar), finite volume schemes.
   class(prism_fnl_object), intent(in)    :: self                                                  !< The equation.
   integer(I4P),            intent(in)    :: ivar                                                  !< Index of scalar variable of q.
   real(R8P),               intent(in)    :: q_gpu(     1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   real(R8P),               intent(inout) :: laplacian_gpu(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Gradient.
   integer(I4P)                           :: i, j, k, b                                            !< Counter.

   associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,dxyz_gpu=>self%field_gpu%dxyz_gpu,&
             hs=>self%numerics%fdv_half_stencils(2))
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,laplacian_gpu,dxyz_gpu)
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      call compute_laplacian_fv_centered(s=hs,dxyz=dxyz_gpu(b,:),q=q_gpu(b,i-hs:i+hs,j-hs:j+hs,k-hs:k+hs,ivar),&
                                         laplacian=laplacian_gpu(b,i,j,k))
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine compute_laplacian_fv

   ! numerical methods, space operators
   subroutine compute_residuals_fd_centered(self, q_gpu, dq_gpu, s)
   !< Compute residuals of equation, space operator, centered finite difference schemes.
   class(prism_fnl_object), intent(inout) :: self               !< The equation.
   real(R8P),               intent(inout) :: q_gpu(1:,     &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1:)              !< Conservative variables.
   real(R8P),               intent(inout) :: dq_gpu(1:,     &
                                                1-self%ngc:,&
                                                1-self%ngc:,&
                                                1-self%ngc:,&
                                                1:)             !< Residuals.
   integer(I4P),  optional, intent(in)    :: s                  !< Stage counter.
   integer(I4P)                           :: i,j,k,b            !< Counter
   real(R8P)                              :: curlD(3), curlB(3) !< Residuals components.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv_c=>self%nv_c,blocks_number=>self%blocks_number, &
             dxyz_gpu=>self%field_gpu%dxyz_gpu,                                                                       &
             s1=>self%numerics%fdv_half_stencils(1),                                                                  &
             s4=>self%numerics%fdv_half_stencils(4),                                                                  &
             var_Jx=>self%physics%var_Jx, var_Jy=>self%physics%var_Jy, var_Jz=>self%physics%var_Jz)
   if (blocks_number > 0) then
      call self%update_ghost(q_gpu=q_gpu, s=s)
      ! compute RHS dD/dt = curl(B/MU0) - J, dB/dt = -curl(D/EPS0)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,dq_gpu)
      do b=1,blocks_number
      do k=1,nk
      do j=1,nj
      do i=1,ni
         call compute_curl_fd_centered(s=s1,dxyz=dxyz_gpu(b,1:3),                             &
                                       q=q_gpu(b,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,VAR_DX:VAR_DZ),&
                                       curl=curlD)
         call compute_curl_fd_centered(s=s1,dxyz=dxyz_gpu(1:3,b),                             &
                                       q=q_gpu(b,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,VAR_BX:VAR_BZ),&
                                       curl=curlB)

         dq_gpu(b,i,j,k,VAR_DX) =  curlB(1)/MU0 - q_gpu(b,i,j,k,var_Jx)
         dq_gpu(b,i,j,k,VAR_DY) =  curlB(2)/MU0 - q_gpu(b,i,j,k,var_Jy)
         dq_gpu(b,i,j,k,VAR_DZ) =  curlB(3)/MU0 - q_gpu(b,i,j,k,var_Jz)
         dq_gpu(b,i,j,k,VAR_BX) = -curlD(1)/EPS0
         dq_gpu(b,i,j,k,VAR_BY) = -curlD(2)/EPS0
         dq_gpu(b,i,j,k,VAR_BZ) = -curlD(3)/EPS0
      enddo
      enddo
      enddo
      enddo
   endif
   endassociate
   endsubroutine compute_residuals_fd_centered

   ! numerical methods, time operators
   subroutine integrate_blanesmoan(self)
   !< Integrate equation, time operator, Blanes and Moan scheme.
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   associate(nc=>self%blanesmoan%nc,a=>self%blanesmoan%a,b=>self%blanesmoan%b)
   call self%compute_coils_current
   ! do s=1, nc
   !    call self%compute_residuals(q=self%q, dq=self%dq)
   !    if (s==1) call self%save_residuals
   !    self%q(VAR_BX,:,:,:,:) = self%q(VAR_BX,:,:,:,:) + b(s) * self%time%dt * self%dq(VAR_BX,:,:,:,:)
   !    self%q(VAR_BY,:,:,:,:) = self%q(VAR_BY,:,:,:,:) + b(s) * self%time%dt * self%dq(VAR_BY,:,:,:,:)
   !    self%q(VAR_BZ,:,:,:,:) = self%q(VAR_BZ,:,:,:,:) + b(s) * self%time%dt * self%dq(VAR_BZ,:,:,:,:)
   !    call self%compute_residuals(q=self%q, dq=self%dq)
   !    self%q(VAR_DX,:,:,:,:) = self%q(VAR_DX,:,:,:,:) + a(s) * self%time%dt * self%dq(VAR_DX,:,:,:,:)
   !    self%q(VAR_DY,:,:,:,:) = self%q(VAR_DY,:,:,:,:) + a(s) * self%time%dt * self%dq(VAR_DY,:,:,:,:)
   !    self%q(VAR_DZ,:,:,:,:) = self%q(VAR_DZ,:,:,:,:) + a(s) * self%time%dt * self%dq(VAR_DZ,:,:,:,:)
   ! enddo
   ! call self%impose_div_free
   endassociate
   endsubroutine integrate_blanesmoan

   subroutine integrate_cfm(self)
   !< Integrate equation, time operator, Commutator-Free Magnus integrator.
   class(prism_fnl_object), intent(inout) :: self             !< The equation.
   real(R8P), parameter                   :: toll=1.0e-14_R8P !< CFM coefficients tollerance.
   integer(I4P)                           :: s,ss             !< Counter.

   ! call self%compute_coils_current
   ! associate(dt=>self%time%dt,s_coeffs=>self%cfm%s_coeffs,e_coeffs=>self%cfm%e_coeffs)
   ! self%cfm%q = self%q
   ! call self%compute_residuals(q=self%cfm%q, dq=self%cfm%dq(:,:,:,:,:,1))
   ! do s=2, self%cfm%n_stages
   !    do ss=1, s-1
   !       if (abs(s_coeffs(s,ss))>toll) &
   !          call self%cfm%compute_exponential_update(alpha=dt*s_coeffs(s,ss),dq=self%cfm%dq(:,:,:,:,:,ss),q=self%cfm%q)
   !    enddo
   !    call self%compute_residuals(q=self%cfm%q, dq=self%cfm%dq(:,:,:,:,:,s))
   ! enddo
   ! self%cfm%q = self%q
   ! do s=1, self%cfm%n_stages
   !    if (abs(self%cfm%e_coeffs(s))>toll) &
   !       call self%cfm%compute_exponential_update(alpha=dt*e_coeffs(s),dq=self%cfm%dq(:,:,:,:,:,s),q=self%cfm%q)
   ! enddo
   ! self%q = self%cfm%q
   ! endassociate
   ! call self%impose_div_free
   endsubroutine integrate_cfm

   subroutine integrate_leapfrog(self)
   !< Integrate equation, time operator, leapfrog scheme.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   ! call self%compute_coils_current
   ! call self%compute_residuals(q=self%q, dq=self%dq)
   ! call self%save_residuals
   ! call self%leapfrog%integrate(dt=self%time%dt, q=self%q, dq=self%dq)
   ! call self%impose_div_free
   endsubroutine integrate_leapfrog

   subroutine integrate_leapfrog_pic(self)
   !< Integrate equation, time operator, leapfrog scheme.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   !!Fai check su come si parlano questa subroutine e quella che calcola la corrente associata alle particelle
   !call self%compute_coils_current

   !! qua ci va la chiamata alla subroutine che aggiorna la neighbour list delle particelle
   !! qua ci va la chiamata alla subroutine che calcola la corrente associata alle particelle

   !call self%compute_residuals(q=self%q, dq=self%dq) !< Calcolo i residui relativi ai campi E e B
   !call self%save_residuals

   !!Qua ci va la chiamata alla subroutine che calcola i resiudi delle particellle dq_pic

   !call self%leapfrog%integrate(dt=self%time%dt, q=self%q, dq=self%dq)

   !!Qua ci va la chiamata alla subroutine che integra aggiornando le velocita e le posizioni delle particelle
   !            !decidi se fare qui all'inizio del tempo successivo l'aggiornamento della neighbour list delle particelle

   !call self%impose_div_free
   !!call self%apply_fWL_correction
   endsubroutine integrate_leapfrog_pic

   subroutine integrate_rk_ls(self)
   !< Integrate equation, time operator, RK classical low storage schemes.
   !< Low storage RK working on q_rk(:,:,:,:,:,1)/q as stages, update q in place.
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   ! call self%compute_coils_current
   ! call self%rk%initialize_stages(q=self%q)
   ! do s=1, self%rk%nrk
   !    call self%compute_residuals(q=self%q, dq=self%dq)
   !    if (s==1) call self%save_residuals
   !    if (self%ib%solids_number>0) then
   !       call self%rk%compute_stage_ls(s=s,dt=self%time%dt,phi=self%ib%phi,dq=self%dq,q=self%q)
   !    else
   !       call self%rk%compute_stage_ls(s=s,dt=self%time%dt,dq=self%dq,q=self%q)
   !    endif
   ! enddo
   ! call self%impose_div_free
   ! call self%apply_fWL_correction
   endsubroutine integrate_rk_ls

   subroutine integrate_rk_ssp(self)
   !< Integrate equation, time operator, SSP RK schemes.
   !< SSP RK working on q_rk as stages.
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   if (self%external_fields%ef_type/=EF_TYPE_NONE) &
      call sub_external_fields_dev(external_fields=self%external_fields, field_gpu=self%field_gpu, &
                                   dt=self%time%dt, time=self%time%time, q_gpu=self%q_gpu)
   call self%rk_gpu%initialize_stages(q_gpu=self%q_gpu)
   do s=1, self%rk%nrk
      call self%compute_coils_current(gamma=self%rk%gamm(s))
      if (self%ib%solids_number>0) then
         call self%rk_gpu%compute_stage(s=s, dt=self%time%dt, phi_gpu=self%ib_gpu%phi_gpu)
      else
         call self%rk_gpu%compute_stage(s=s, dt=self%time%dt)
      endif
      call self%compute_residuals(q_gpu=self%rk_gpu%q_rk_gpu(:,:,:,:,:,s), dq_gpu=self%dq_gpu, s=s)
      if (s==1) call self%save_residuals
      if (self%ib%solids_number>0) then
         call self%rk_gpu%assign_stage(s=s, q_gpu=self%dq_gpu, phi_gpu=self%ib_gpu%phi_gpu)
      else
         call self%rk_gpu%assign_stage(s=s, q_gpu=self%dq_gpu)
      endif
   enddo
   if (self%ib%solids_number>0) then
      call self%rk_gpu%update_q(dt=self%time%dt, phi_gpu=self%ib_gpu%phi_gpu, q_gpu=self%q_gpu)
      call self%update_rk_ghost(dt=self%time%dt, phi_gpu=self%ib_gpu%phi_gpu)
   else
      call self%rk_gpu%update_q(dt=self%time%dt, q_gpu=self%q_gpu)
      call self%update_rk_ghost(dt=self%time%dt)
   endif
   call self%impose_div_free
   ! call self%apply_fwl_correction  ! to be removed, probably
   if (self%external_fields%ef_type/=EF_TYPE_NONE) &
      call add_external_fields_dev(external_fields=self%external_fields, field_gpu=self%field_gpu, &
                                   dt=self%time%dt, time=self%time%time, q_gpu=self%q_gpu)
   endsubroutine integrate_rk_ssp

   subroutine integrate_rk_yoshida(self)
   !< Integrate equation, time operator, Yoshida RK scheme.
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   ! call self%compute_coils_current
   ! do s=1, self%rk%nrk - 1
   !    call self%compute_residuals(q=self%q, dq=self%dq)
   !    if (s==1) call self%save_residuals
   !    self%q(VAR_BX,:,:,:,:) = self%q(VAR_BX,:,:,:,:) + self%rk%ssa(s) * self%time%dt * self%dq(VAR_BX,:,:,:,:)
   !    self%q(VAR_BY,:,:,:,:) = self%q(VAR_BY,:,:,:,:) + self%rk%ssa(s) * self%time%dt * self%dq(VAR_BY,:,:,:,:)
   !    self%q(VAR_BZ,:,:,:,:) = self%q(VAR_BZ,:,:,:,:) + self%rk%ssa(s) * self%time%dt * self%dq(VAR_BZ,:,:,:,:)
   !    call self%compute_residuals(q=self%q, dq=self%dq)
   !    self%q(VAR_DX,:,:,:,:) = self%q(VAR_DX,:,:,:,:) + self%rk%ssb(s) * self%time%dt * self%dq(VAR_DX,:,:,:,:)
   !    self%q(VAR_DY,:,:,:,:) = self%q(VAR_DY,:,:,:,:) + self%rk%ssb(s) * self%time%dt * self%dq(VAR_DY,:,:,:,:)
   !    self%q(VAR_DZ,:,:,:,:) = self%q(VAR_DZ,:,:,:,:) + self%rk%ssb(s) * self%time%dt * self%dq(VAR_DZ,:,:,:,:)
   ! enddo
   ! call self%compute_residuals(q=self%q, dq=self%dq)
   ! self%q(VAR_BX,:,:,:,:) = self%q(VAR_BX,:,:,:,:) + self%rk%ssa(self%rk%nrk) * self%time%dt * self%dq(VAR_BX,:,:,:,:)
   ! self%q(VAR_BY,:,:,:,:) = self%q(VAR_BY,:,:,:,:) + self%rk%ssa(self%rk%nrk) * self%time%dt * self%dq(VAR_BY,:,:,:,:)
   ! self%q(VAR_BZ,:,:,:,:) = self%q(VAR_BZ,:,:,:,:) + self%rk%ssa(self%rk%nrk) * self%time%dt * self%dq(VAR_BZ,:,:,:,:)
   ! call self%impose_div_free
   endsubroutine integrate_rk_yoshida

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
         !if (self%ib%solids_number > 0) call self%compute_phi()
         !call self%amr_update
      enddo
      call self%set_initial_conditions
      self%time%time = 0._R8P
      self%time%it = 0
      call self%mpih_gpu%print_message('impose initial conditions finish')
   endif
   !if (self%ib%solids_number > 0) call self%compute_phi()
   ! call self%amr_update
   call self%compute_divergence(ivar=1,q_gpu=self%q_gpu,divergence_gpu=self%divergence_gpu(1,:,:,:,:))
   call self%compute_divergence(ivar=4,q_gpu=self%q_gpu,divergence_gpu=self%divergence_gpu(2,:,:,:,:))
   call self%compute_divergence(ivar=7,q_gpu=self%q_gpu,divergence_gpu=self%divergence_gpu(3,:,:,:,:))
   call self%save_simulation_data
   call self%compute_energy
   call self%save_energy_error(is_to_open=.true.)
   call self%io%open_file_residuals(nv=self%nv)

   if (self%numerics%scheme_time==NUM_SCHEME_TIME_LEAPFROG) then
      ! to be implemented leapfrog on device
      ! call self%leapfrog%assign_step(s=1, q=self%q)
      ! call self%compute_dt
      ! call self%compute_residuals(q=self%q, dq=self%dq)
      ! self%q = self%q + self%time%dt * self%dq
   endif

   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
      if (self%pic%scheme_time==NUM_SCHEME_TIME_PIC_LEAPFROG) then
         ! to be implemented
      endif
   endif

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
         !call self%amr_update
         call self%mpih_gpu%barrier(tictoc=.true.)
      endif

      call self%compute_dt
      if ((self%time%it_max <= 0).and.(self%time%time+self%time%dt > self%time%time_max)) &
         self%time%dt=self%time%time_max-self%time%time

      call self%integrate

      self%time%time = self%time%time + self%time%dt
      call self%time%print_progress(nodes_number=self%adam%tree%nodes_number)

      call self%save_simulation_data
      call self%compute_energy
      call self%save_energy_error

      if (((self%time%it_max <= 0).and.(self%time%time >= self%time%time_max)).or.&
         ((self%time%it>=self%time%it_max).and.(self%time%it_max > 0))) exit integration

      call self%mpih_gpu%barrier(tictoc=.true., timing=timing_step(2), single=.true.)
   enddo integration
   call self%mpih_gpu%barrier(tictoc=.true., timing=timing(2), single=.true.)
   call self%compute_energy_error
   call self%save_simulation_data
   call self%io%close_file_residuals
   call self%save_energy_error(is_to_close=.true.)
   call self%mpih_gpu%print_message('Initial/final energy of D field: '//trim(str(sqrt(self%energy_D(1))))//' '//&
                                                                         trim(str(sqrt(self%energy_D(size(self%energy_D))))))
   call self%mpih_gpu%print_message('Initial/final energy of B field: '//trim(str(sqrt(self%energy_B(1))))//' '//&
                                                                         trim(str(sqrt(self%energy_D(size(self%energy_B))))))
   call self%mpih_gpu%print_message('RMS Error of D field: '//trim(str(self%rms_energy_error_D)))
   call self%mpih_gpu%print_message('RMS Error of B field: '//trim(str(self%rms_energy_error_B)))
   call self%mpih_gpu%finalize
   endsubroutine simulate

   ! numerical methods, miscellanea
   subroutine compute_auxiliary_fields(self)
   !< Compute auxiliary fields.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   if (self%io%save_divergence_fields) then
      call self%compute_divergence(ivar=VAR_DX,q_gpu=self%q_gpu,divergence_gpu=self%divergence_gpu(1,:,:,:,:))
      call self%compute_divergence(ivar=VAR_BX,q_gpu=self%q_gpu,divergence_gpu=self%divergence_gpu(2,:,:,:,:))
      ! call self%compute_divergence(ivar=7,q=self%q,divergence=self%divergence(3,:,:,:,:))
   endif
   if (self%io%save_curl_fields) then
      call self%compute_curl(ivar=VAR_DX,q_gpu=self%q_gpu,curl_gpu=self%curl_gpu(1:3,:,:,:,:))
      call self%compute_curl(ivar=VAR_BX,q_gpu=self%q_gpu,curl_gpu=self%curl_gpu(4:6,:,:,:,:))
      ! call self%compute_curl(ivar=7,q=self%q,curl=self%curl(7:9,:,:,:,:))
   endif
   endsubroutine compute_auxiliary_fields

   subroutine compute_dt(self)
   !< Compute maximum time step accordingly to CFL stabilty criterion.
   class(prism_fnl_object), intent(inout) :: self     !< The equation.
   real(R8P)                              :: dxyz_min !< Minimum space step.
   integer(I4P)                           :: b        !< Counter.

   associate(blocks_number=>self%blocks_number, dxyz_gpu=>self%field_gpu%dxyz_gpu, evmax=>self%physics%evmax)
   dxyz_min = huge(0._R8P)
   !$acc parallel loop independent gang vector DEVICEVAR(dxyz_gpu) reduction(min: dxyz_min)
   do b=1, blocks_number
      dxyz_min = min(dxyz_min, dxyz_gpu(b,1), dxyz_gpu(b,2), dxyz_gpu(b,3))
   enddo
   dxyz_min = dxyz_min * 0.5_R8P
   self%time%dt = self%time%CFL*dxyz_min / evmax
   call MPI_ALLREDUCE(MPI_IN_PLACE, self%time%dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, self%mpih_gpu%error)
   endassociate
   endsubroutine compute_dt

   subroutine compute_energy(self)
   !< Compute energy.
   class(prism_fnl_object), intent(inout) :: self     !< The equation.
   real(R8P)                              :: energy_D !< Energy of D field.
   real(R8P)                              :: energy_B !< Energy of B field.

   call compute_e(ivar=VAR_DX, ngc=self%ngc, q_gpu=self%q_gpu, energy=energy_D)
   call compute_e(ivar=VAR_BX, ngc=self%ngc, q_gpu=self%q_gpu, energy=energy_B)
   call MPI_ALLREDUCE(MPI_IN_PLACE, energy_D, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, self%mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, energy_B, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, self%mpih%error)
   if (allocated(self%energy_D).and.allocated(self%energy_B)) then
      self%energy_D = [self%energy_D, energy_D]
      self%energy_B = [self%energy_B, energy_B]
   else
      allocate(self%energy_D(1:self%time%it))
      allocate(self%energy_B(1:self%time%it))
      self%energy_D = energy_D
      self%energy_B = energy_B
   endif
   contains
      subroutine compute_e(ivar, ngc, q_gpu, energy)
      !< Compute energy of vector field starting from ivar.
      integer(I4P), intent(in)  :: ivar    !< Starting position of vector field.
      integer(I4P), intent(in)  :: ngc     !< Number of ghost cells.
      real(R8P),    intent(in)  :: q_gpu(1:,&
                                     1-ngc:,&
                                     1-ngc:,&
                                     1-ngc:,&
                                     1:)   !< Conservative variables.
      real(R8P),    intent(out) :: energy  !< Energy of the vector field starting from ivar.
      integer(I4P)              :: i,j,k,b !< Counter.

      energy = 0.0_R8P
      associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,blocks_number=>self%blocks_number)
      !$acc parallel loop independent gang vector DEVICEVAR(q_gpu) reduction(+: energy)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         energy = energy + 0.5_R8P * (q_gpu(b,i,j,k,ivar  )*q_gpu(b,i,j,k,ivar  ) + &
                                      q_gpu(b,i,j,k,ivar+1)*q_gpu(b,i,j,k,ivar+1) + &
                                      q_gpu(b,i,j,k,ivar+2)*q_gpu(b,i,j,k,ivar+2))
      enddo
      enddo
      enddo
      enddo
      endassociate
      endsubroutine compute_e
   endsubroutine compute_energy

   subroutine compute_energy_error(self)
   !< Compute energy error.
   class(prism_fnl_object), intent(inout) :: self       !< The equation.
   real(R8P)                              :: energy_D0  !< Initial energy of D field.
   real(R8P)                              :: energy_B0  !< Initial energy of B field.
   real(R8P), allocatable                 :: error_D(:) !< Error (normalized) of energy of D field.
   real(R8P), allocatable                 :: error_B(:) !< Error (normalized) of energy of B field.

   if (allocated(self%energy_D).and.allocated(self%energy_B)) then
      energy_D0 = self%energy_D(1)
      energy_B0 = self%energy_B(1)

      error_D = abs(self%energy_D - energy_D0) / abs(energy_D0)
      error_B = abs(self%energy_B - energy_B0) / abs(energy_B0)

      self%rms_energy_error_D = sqrt(sum(error_D)/size(error_D))
      self%rms_energy_error_B = sqrt(sum(error_B)/size(error_B))
   endif
   endsubroutine compute_energy_error

   subroutine impose_ct_correction(self, ivar)
   !< Impose Constrained Transport Correction on vectorial variable q(ivar:ivar+2).
   !< Note that self%divergence memory is used as buffer, be carefull.
   class(prism_fnl_object), intent(inout) :: self      !< The equation.
   integer(I4P),            intent(in)    :: ivar      !< Variable (start) index in q.
   real(R8P)                              :: dq_max    !< Maximum residual.
   integer(I4P)                           :: iter      !< Counter.
   integer(I4P)                           :: i,j,k,b,v !< Counter.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, blocks_number=>self%blocks_number, buffer=>self%divergence_gpu,&
             q_gpu=>self%q_gpu)
   if (blocks_number>0) then
      call self%compute_divergence(ivar=ivar,q_gpu=q_gpu,divergence_gpu=buffer(:,:,:,:,4))
      do iter=1, self%flail%iterations
         ! call compute_smoothing_multigrid(ni=ni,nj=nj,nk=nk,ngc=ngc,nv=1_I4P,blocks_number=blocks_number, &
         !                                  dxyz=self%field%dxyz,                                           &
         !                                  f=-buffer(4:4,:,:,:,:),                                         &
         !                                  q=buffer(7:7,:,:,:,:),                                          &
         !                                  dq_max=dq_max,                                                  &
         !                                  dq=buffer(5:5,:,:,:,:),                                         &
         !                                  iterations_init=self%flail%iterations_init,                     &
         !                                  iterations_fine=self%flail%iterations_fine,                     &
         !                                  iterations_coarse=self%flail%iterations_coarse)
         if (dq_max < self%flail%tolerance) exit
      enddo
      call self%mpih_gpu%print_message('FLAIL convergence reached at iteration '//trim(str(iter,.true.)))
      call self%compute_gradient(ivar=1,q_gpu=buffer(:,:,:,:,7:7),gradient_gpu=buffer(:,:,:,:,4:6))
      !$acc parallel loop independent gang vector collapse(5) DEVICEVAR(q_gpu,buffer)
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, nj
                  do v=1, 3
                     q_gpu(b,i,j,k,ivar+v-1) = q_gpu(b,i,j,k,ivar+v-1) + buffer(b,i,j,k,3+v)
                  enddo
               enddo
            enddo
         enddo
      enddo
   endif
   endassociate
   endsubroutine impose_ct_correction

   subroutine impose_div_free(self)
   !< Impose divergence-free property.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   associate(constrained_transport_D=>self%numerics%constrained_transport_D,&
             constrained_transport_B=>self%numerics%constrained_transport_B,div_corr_var=>self%numerics%div_corr_var)
   if (constrained_transport_D.and.div_corr_var==DIV_CORR_VAR_POISS) call self%impose_ct_correction(ivar=1_I4P)
   if (constrained_transport_B.and.div_corr_var==DIV_CORR_VAR_POISS) call self%impose_ct_correction(ivar=4_I4P)
   ! here should go also other corrections...
   endassociate
   endsubroutine impose_div_free
endmodule adam_prism_fnl_object
