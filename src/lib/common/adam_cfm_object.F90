!< ADAM, Commutator-Free Magnus integrators class definition.
module adam_cfm_object
!< ADAM, Commutator-Free Magnus integrators class definition.

use adam_field_object
use adam_global_mpih, only: mpih
use adam_global_grid, only: grid
use finer
use penf

implicit none
save
private
public :: cfm_object
public :: CFM_4
public :: CFM_6
public :: CFM_8

character(len=5), parameter :: CFM_4="CFM-4" !< Parameter of time scheme, Commutator-Free Magnus 4th order.
character(len=5), parameter :: CFM_6="CFM-6" !< Parameter of time scheme, Commutator-Free Magnus 4th order.
character(len=5), parameter :: CFM_8="CFM-8" !< Parameter of time scheme, Commutator-Free Magnus 4th order.

character(len=22), parameter :: INI_SECTION_NAME="commutator_free_magnus" !< INI (config) file section name containing configs.

type :: cfm_object
   !< Commutator-Free Magnus integratos class definition.
   integer(I4P)              :: order           !< Order of the scheme.
   integer(I4P)              :: n_stages        !< Number of stages.
   integer(I4P)              :: n_exponentials  !< Number of exponential integration combinations.
   real(R8P),    allocatable :: s_coeffs(:,:)   !< Stage coefficients.
   real(R8P),    allocatable :: e_coeffs(:)     !< Final exponential coefficients.
   character(:), allocatable :: scheme          !< CFM scheme name.
   real(R8P),    allocatable :: dq(:,:,:,:,:,:) !< Field cell centered variables residuals, CFM stages.
   real(R8P),    allocatable ::  q(:,:,:,:,:)   !< Field cell centered variables, CFM buffer.
   ! grid/field data replica for easy handling
   type(field_object), pointer :: field=>null()         !< The field.
   integer(I4P),       pointer :: ngc=>null()           !< Number of ghost cells.
   integer(I4P),       pointer :: ni=>null()            !< Number of cells in i direction.
   integer(I4P),       pointer :: nj=>null()            !< Number of cells in j direction.
   integer(I4P),       pointer :: nk=>null()            !< Number of cells in k direction.
   integer(I4P),       pointer :: nb=>null()            !< Total blocks number for MPI.
   integer(I4P),       pointer :: blocks_number=>null() !< Actual blocks number.
   integer(I4P),       pointer :: ns=>null()            !< Number of fluids specie.
   integer(I4P),       pointer :: nv=>null()            !< Number of conservative variables.
   contains
      ! public methods
      procedure, pass(self) :: compute_exponential_update !< Compute exponential update.
      procedure, pass(self) :: description                !< Return pretty-printed object description.
      procedure, pass(self) :: initialize                 !< Initialize class.
      procedure, pass(self) :: initialize_CFM_4           !< Initialize CFM 4th order.
      procedure, pass(self) :: initialize_CFM_6           !< Initialize CFM 6th order.
      procedure, pass(self) :: initialize_CFM_8           !< Initialize CFM 8th order.
      procedure, pass(self) :: load_from_file             !< Load config from file.
endtype cfm_object
contains
   ! public methods
   subroutine compute_exponential_update(self, alpha, dq, q)
   !< Compute exponential update q(n+1) = exp(alpha * R(q(n))) * q(n).
   !< For linear PDE q(n+1) = q(n) + alpha * R(q(n))
   class(cfm_object), intent(in)    :: self   !< CFM object.
   real(R8P),         intent(in)    :: alpha  !< Pseudo time step.
   real(R8P),         intent(in)    :: dq(1:,          &
                                          1-self%ngc:, &
                                          1-self%ngc:, &
                                          1-self%ngc:, &
                                          1:) !< Conservative variables residuals.
   real(R8P),         intent(inout) :: q(1:,          &
                                         1-self%ngc:, &
                                         1-self%ngc:, &
                                         1-self%ngc:, &
                                         1:)  !< Conservative variables.

   q = q + alpha * dq
   endsubroutine compute_exponential_update

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(cfm_object), intent(in) :: self             !< CFM object.
   character(len=:), allocatable :: desc             !< Description.
   character(len=1), parameter   :: NL=new_line('a') !< New line character.

   desc =       mpih%myrankstr//'Commutator-Free Magnus scheme main data'                   //NL
   desc = desc//mpih%myrankstr//'  scheme name:            '//trim(    self%scheme         )//NL
   desc = desc//mpih%myrankstr//'  order:                  '//trim(str(self%order         ))//NL
   desc = desc//mpih%myrankstr//'  number of stages:       '//trim(str(self%n_stages      ))//NL
   desc = desc//mpih%myrankstr//'  number of exponentials: '//trim(str(self%n_exponentials))
   endfunction description

   subroutine initialize(self, file_parameters, scheme, field)
   !< Initialize class.
   class(cfm_object),   intent(inout)        :: self            !< CFM object.
   type(file_ini),      intent(in), optional :: file_parameters !< Simulation parameters ini file handler.
   character(*),        intent(in), optional :: scheme          !< CFM scheme.
   type(field_object),  intent(in), target   :: field           !< The field.

   call mpih%print_message('cfm_object%initialize start')
   call associate_adam_data(field=field)
   if (present(file_parameters)) then
      call self%load_from_file(file_parameters=file_parameters)
   elseif (present(scheme)) then
      self%scheme = trim(adjustl(scheme))
   else
      call mpih%error_stop(msg=': failed to initialize cfm object, one between file parameters and scheme must be passed')
   endif
   select case(self%scheme)
   case(CFM_4)
      call self%initialize_CFM_4
   case(CFM_6)
      call mpih%error_stop(msg=': 6th order CFM scheme "'//self%scheme //'" not yet implemented')
      call self%initialize_CFM_6
   case(CFM_8)
      call mpih%error_stop(msg=': 8th order CFM scheme "'//self%scheme //'" not yet implemented')
      call self%initialize_CFM_8
   case default
      call mpih%error_stop(msg=': CFM scheme "'//self%scheme //'" unknown')
   endselect
   print '(A)', self%description()
   call mpih%print_message('cfm_object%initialize finish')
   contains
      subroutine associate_adam_data(field)
      !< Associate objects data for easy handling.
      type(field_object), intent(in), target :: field !< The field.

      self%field         => field
      self%blocks_number => field%blocks_number
      self%ni            => grid%ni
      self%nj            => grid%nj
      self%nk            => grid%nk
      self%ngc           => grid%ngc
      self%nb            => field%nb
      self%nv            => field%nv
      endsubroutine associate_adam_data
   endsubroutine initialize

   subroutine initialize_CFM_4(self)
   !< Initialize CFM 4th order.
   class(cfm_object), intent(inout) :: self                !< CFM object.
   real(R8P), parameter             :: sqrt3=sqrt(3.0_R8P) !< Square root of 3.

   self%order          = 4_I4P
   self%n_stages       = 4_I4P
   self%n_exponentials = 4_I4P
   associate(ns=>self%n_stages,ne=>self%n_exponentials,ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,nv=>self%nv,nb=>self%nb)
   if (allocated(self%s_coeffs)) deallocate(self%s_coeffs) ; allocate(self%s_coeffs(ns,ns))
   if (allocated(self%e_coeffs)) deallocate(self%e_coeffs) ; allocate(self%e_coeffs(ne   ))
   call allocate_variable(var=self%dq,               &
                          ulb=reshape([1,nv,         &
                                       1-ngc,ni+ngc, &
                                       1-ngc,nj+ngc, &
                                       1-ngc,nk+ngc, &
                                       1,nb,         &
                                       1,ns],[2,6]), &
                          msg=mpih%myrankstr//'cfm_object%initialize allocate dq')
   call allocate_variable(var=self%q,                &
                          ulb=reshape([1,nv,         &
                                       1-ngc,ni+ngc, &
                                       1-ngc,nj+ngc, &
                                       1-ngc,nk+ngc, &
                                       1,nb],[2,5]), &
                          msg=mpih%myrankstr//'cfm_object%initialize allocate q')
   endassociate
   self%q        = 0.0_R8P
   self%s_coeffs = 0.0_R8P
   self%e_coeffs = 0.0_R8P
   ! final exponential coefficients (symmetric)
   self%e_coeffs(1) = (3.0_R8P - 2.0_R8P*sqrt3) / 12.0_R8P
   self%e_coeffs(2) = (3.0_R8P + 2.0_R8P*sqrt3) / 12.0_R8P
   self%e_coeffs(3) = (3.0_R8P + 2.0_R8P*sqrt3) / 12.0_R8P
   self%e_coeffs(4) = (3.0_R8P - 2.0_R8P*sqrt3) / 12.0_R8P
   ! stage coefficients
   self%s_coeffs(2,1) = 0.5_R8P - sqrt3/6.0_R8P
   self%s_coeffs(3,1) =           sqrt3/3.0_R8P
   self%s_coeffs(3,2) = 0.5_R8P - sqrt3/6.0_R8P
   self%s_coeffs(4,1) = 0.5_R8P - sqrt3/6.0_R8P
   self%s_coeffs(4,2) =           sqrt3/3.0_R8P
   self%s_coeffs(4,3) = 0.5_R8P - sqrt3/6.0_R8P
   endsubroutine initialize_CFM_4

   subroutine initialize_CFM_6(self)
   !< Initialize CFM 6th order.
   class(cfm_object), intent(inout) :: self !< CFM object.

   ! to be implemented
   endsubroutine initialize_CFM_6

   subroutine initialize_CFM_8(self)
   !< Initialize CFM 8th order.
   class(cfm_object), intent(inout) :: self !< CFM object.

   ! to be implemented
   endsubroutine initialize_CFM_8

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(cfm_object), intent(inout)        :: self            !< CFM object.
   type(file_ini),    intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,           intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                 :: go_on_fail_     !< Go on if load fails.
   character(99)                           :: buff_c          !< Character buffer.
   integer(I4P)                            :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='scheme', val=buff_c, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(scheme)')
   self%scheme = trim(adjustl(buff_c))
   endsubroutine load_from_file
endmodule adam_cfm_object
