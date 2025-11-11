!< ADAM, PRISM (Plasma Research usIng Simulation Methods) numerics class definition, common backend.
module adam_prism_numerics_object
!< ADAM, PRISM (Plasma Research usIng Simulation Methods) numerics class definition, common backend.
!<

! ADAM modules
use :: adam_mpih_object, only : mpih_object
! PRISM modules
use :: adam_prism_parameters
! third party modules
use :: finer, only : file_ini
use :: penf, only : I4P, R8P, str

implicit none
private
public :: prism_numerics_object

character(len=8), parameter :: INI_SECTION_NAME='numerics' !< INI file section name containing fluid numerics.

character(11), parameter, public :: NUM_SCHEME_TIME_BLANES_MOAN='BLANES_MOAN' !< Blanes-Moan numerical scheme for time operator.
character(8),  parameter, public :: NUM_SCHEME_TIME_LEAPFROG='LEAPFROG'       !< Leapfrog numerical scheme for time operator.
character(11), parameter, public :: NUM_SCHEME_TIME_RUNGE_KUTTA='RUNGE_KUTTA' !< Runge-Kutta numerical scheme for time operator.
character(11), parameter, public :: NUM_SCHEME_SPACE_FD_CENTERED='FD_CENTERED'!< FD centered scheme for space operator.
character(11), parameter, public :: NUM_SCHEME_SPACE_FV_CENTERED='FV_CENTERED'!< FV centered scheme for space operator.
character(4),  parameter, public :: NUM_SCHEME_SPACE_WENO='WENO'              !< WENO numerical scheme for space operator.
character(12), parameter, public :: RECONSTRUCTION_VARS_CONS='CONSERVATIVE'   !< High-order reconstruction on conservative vars.
character(15), parameter, public :: RECONSTRUCTION_VARS_CHAR='CHARACTERISTICS'!< High-order reconstruction on characteristics vars.
character(7),  parameter, public :: DIV_CORR_VAR_POISS='POISSON'              !< Poisson divergence correction.
character(10), parameter, public :: DIV_CORR_VAR_HYPER='HYPERBOLIC'           !< Hyperbolic divergence correction.

type :: prism_numerics_object
   !< PRISM numerics class definition.
   type(mpih_object)         :: mpih                            !< MPI handler.
   character(:), allocatable :: scheme_time                     !< Numerical scheme for time operator [runge_kutta, leapfrog,...].
   character(:), allocatable :: scheme_space                    !< Numerical scheme for space operator [weno, centered].
   integer(I4P)              :: fdv_order=2_I4P                 !< Order of finite difference/volume schemes.
   integer(I4P)              :: fdv_half_stencil=1_I4P          !< Half stencil length of finite difference/volume schemes.
   character(:), allocatable :: reconstruction_vars             !< Type of WENO reconstruction variables (cons., charct.,...).
   logical                   :: constrained_transport_D=.false. !< Enable Constrained Transport Correction on D.
   logical                   :: constrained_transport_B=.false. !< Enable Constrained Transport Correction on D.
   ! logical                   :: d_divergence_cleaner=.false. !< Enable electric field divergence cleaning.
   ! logical                   :: b_divergence_cleaner=.false. !< Enable magnetic field divergence cleaning.
   character(:), allocatable :: div_corr_var                    !< Type of divergence correction variables (poisson, hyperbolic,...).
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize numerics.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype prism_numerics_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_numerics_object), intent(in) :: self             !< numerics.
   character(len=:), allocatable            :: desc             !< Description.
   character(len=1), parameter              :: NL=new_line('a') !< New line character.

   desc =       self%mpih%myrankstr//'numerics main data:'                                                                    //NL
   desc = desc//self%mpih%myrankstr//'  Numerical scheme for time operator:               '//self%scheme_time                 //NL
   desc = desc//self%mpih%myrankstr//'  Numerical scheme for space operator:              '//self%scheme_space                //NL
   desc = desc//self%mpih%myrankstr//'  Order of finite difference/volume schemes:        '//trim(str(self%fdv_order))        //NL
   desc = desc//self%mpih%myrankstr//'  Half stencil of finite difference/volume schemes: '//trim(str(self%fdv_half_stencil)) //NL
   desc = desc//self%mpih%myrankstr//'  Reconstruction variables:                         '//self%reconstruction_vars         //NL
   desc = desc//self%mpih%myrankstr//'  Enable Constrained Transport Correction on D:     '//str(self%constrained_transport_D)//NL
   desc = desc//self%mpih%myrankstr//'  Enable Constrained Transport Correction on B:     '//str(self%constrained_transport_B)//NL
   desc = desc//self%mpih%myrankstr//'  Divergence correction:                            '//self%div_corr_var
   ! desc = desc//self%mpih%myrankstr//'  D divergence correction:                     '//trim(str(self%d_divergence_cleaner))//NL
   ! desc = desc//self%mpih%myrankstr//'  B divergence correction:                     '//trim(str(self%b_divergence_cleaner))//NL
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize the equation.
   class(prism_numerics_object), intent(inout) :: self            !< numerics.
   type(file_ini),               intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'prism_numerics_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'prism_numerics_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(prism_numerics_object), intent(inout)        :: self            !< numerics.
   type(file_ini),               intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                      intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                            :: go_on_fail_     !< Go on if load fails.
   integer(I4P)                                       :: s               !< Counter.
   integer(I4P)                                       :: error           !< Error status.
   character(99)                                      :: buff            !< Character buffer.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='scheme_time', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(scheme_time)')
   select case(trim(adjustl(buff)))
   case('BLANES_MOAN', 'blanes_moan', 'Blanes_Moan')
      self%scheme_time = NUM_SCHEME_TIME_BLANES_MOAN
   case('LEAPFROG', 'leapfrog', 'Leapfrog')
      self%scheme_time = NUM_SCHEME_TIME_LEAPFROG
   case('RUNGE_KUTTA', 'runge_kutta', 'Runge_Kutta')
      self%scheme_time = NUM_SCHEME_TIME_RUNGE_KUTTA
   case default
      call self%mpih%print_message(msg='warning: numerical scheme "'//trim(adjustl(buff))//'" unknown. Revert back to RK scheme')
      self%scheme_time = NUM_SCHEME_TIME_RUNGE_KUTTA
   endselect

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='scheme_space', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(scheme_space)')
   select case(trim(adjustl(buff)))
   case('WENO', 'weno', 'Weno')
      self%scheme_space = NUM_SCHEME_SPACE_WENO
   case('FD_CENTERED', 'fd_centered', 'FD_Centered')
      self%scheme_space = NUM_SCHEME_SPACE_FD_CENTERED
   case('FV_CENTERED', 'fv_centered', 'FV_Centered')
      self%scheme_space = NUM_SCHEME_SPACE_FV_CENTERED
   case default
      call self%mpih%print_message(msg='warning: numerical scheme "'//trim(adjustl(buff))//'" unknown. Revert back to WENO scheme')
      self%scheme_space = NUM_SCHEME_SPACE_WENO
   endselect

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='fdv_order', val=self%fdv_order,error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(fdv_order)')
   self%fdv_half_stencil = self%fdv_order / 2 ! valid only for centered schems

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='reconstruction_variables', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(reconstruction_variables)')
   select case(trim(adjustl(buff)))
   case('CONSERVATIVE', 'conservative', 'Conservative')
      self%reconstruction_vars = RECONSTRUCTION_VARS_CONS
   case('CHARACTERISTICS', 'characteristics', 'Characteristics')
      self%reconstruction_vars = RECONSTRUCTION_VARS_CHAR
   case default
      call self%mpih%print_message(msg='warning: WENO reconstruction variable "'//trim(adjustl(buff))// &
                                   '" unknown. Revert back to conservative variables WENO reconstruction')
      self%reconstruction_vars = RECONSTRUCTION_VARS_CONS
   endselect

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='constrained_transport', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(constrained_transport)')
   select case(trim(adjustl(buff)))
   case('NO', 'no', 'No', 'nO')
      self%constrained_transport_D = .false.
      self%constrained_transport_B = .false.
   case('D', 'd')
      self%constrained_transport_D = .true.
      self%constrained_transport_B = .false.
   case('B', 'b')
      self%constrained_transport_D = .false.
      self%constrained_transport_B = .true.
   case('DB', 'db', 'Db', 'dB', 'BD', 'bd', 'bD', 'Bd')
      self%constrained_transport_D = .true.
      self%constrained_transport_B = .true.
   case default
      call self%mpih%print_message(msg='warning: constrained transport config "'//trim(adjustl(buff))//'" ambiguos, disabled!')
      self%constrained_transport_D = .false.
      self%constrained_transport_B = .false.
   endselect

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='divergence_correction', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(divergence_correction)')
   select case(trim(adjustl(buff)))
   case('POISSON', 'poisson', 'Poisson')
      self%div_corr_var = DIV_CORR_VAR_POISS
   case('HYPERBOLIC', 'hyperbolic', 'Hyperbolic')
      self%div_corr_var = DIV_CORR_VAR_HYPER
   case default
      call self%mpih%print_message(msg='warning: divergence correction variable not activated: divergence correction disabled!')
      self%div_corr_var = 'No'
      self%constrained_transport_D = .false.
      self%constrained_transport_B = .false.
   endselect

   ! call file_parameters%get(section_name=INI_SECTION_NAME, option_name='d_divergence_cleaner', &
   !                          val=self%d_divergence_cleaner,error=error)
   ! if (.not.go_on_fail_.and.error>0) &
   !    call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(d_divergence_cleaner)')
   ! call file_parameters%get(section_name=INI_SECTION_NAME, option_name='b_divergence_cleaner', &
   !                          val=self%b_divergence_cleaner,error=error)
   ! if (.not.go_on_fail_.and.error>0) &
   !    call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(b_divergence_cleaner)')
   ! call file_parameters%get(section_name=INI_SECTION_NAME, option_name='chi', val=self%chi, error=error)
   ! if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(chi)')
   ! call file_parameters%get(section_name=INI_SECTION_NAME, option_name='eta', val=self%eta, error=error)
   ! if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(eta)')
   ! ! increase variables number if q if cleaners are used
   ! if     (self%d_divergence_cleaner.and..not.self%b_divergence_cleaner.and.self%div_corr_var==DIV_CORR_VAR_HYPER)then
   !    self%nv = self%nv + 1_I4P
   !    self%nv_cl = 1_I4P
   ! elseif (self%d_divergence_cleaner.and.self%b_divergence_cleaner.and.self%div_corr_var == DIV_CORR_VAR_HYPER) then
   !    self%nv = self%nv + 2_I4P
   !    self%nv_cl = 2_I4P
   ! endif
   ! ! consistency check
   ! if     (.not.self%d_divergence_cleaner.and..not.self%b_divergence_cleaner.and.self%nv/=9_I4P) then
   !    call self%mpih%error_stop(msg=': D and B divergence cleaners are false but nv /= 9')
   ! elseif (.not.self%d_divergence_cleaner.and.self%b_divergence_cleaner) then
   !    call self%mpih%error_stop(msg=': D divergence cleaner is false but B divergence cleaner is true')
   ! elseif (self%d_divergence_cleaner.and.self%chi<=1._R8P) then
   !    call self%mpih%error_stop(msg=': if D divergence cleaner is true chi cannot be lower than 1.0')
   ! endif
   ! self%evmax = sqrt(1._R8P/(EPS0*MU0))
   ! if (self%d_divergence_cleaner .and. self%div_corr_var == DIV_CORR_VAR_HYPER) then
   !    self%evmax = self%chi*sqrt(1._R8P/(EPS0*MU0))
   ! endif
   endsubroutine load_from_file
endmodule adam_prism_numerics_object
