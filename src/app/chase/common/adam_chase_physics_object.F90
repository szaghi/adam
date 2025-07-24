!< ADAM, CHASE physics class definition, common CPU backend.
module adam_chase_physics_object
!< ADAM, CHASE physics class definition, common CPU backend.
!<
!< Field variables are arranged as follows:
!<```
!< q(1): r
!< q(2): r * u
!< q(3): r * v
!< q(4): r * w
!< q(5): r * E
!<```
!< Field variables fluxes are:
!<```
!< Fx(1) = r * u       Fy(1) = r * v       Fz(1) = r * w
!< Fx(2) = r * u^2 + p Fy(2) = r * u * v   Fz(2) = r * u * w
!< Fx(3) = r * v * u   Fy(3) = r * v^2 + p Fz(3) = r * v * w
!< Fx(4) = r * w * u   Fy(4) = r * w * v   Fz(4) = r * w^2 + p
!< Fx(5) = r * u * H   Fy(5) = r * v * H   Fz(5) = r * w * H
!<```
!< Field auxiliary variables are arranged as follows:
!<```
!< q_aux(1 ): r  (density)              |
!< q_aux(2 ): u  (x velocity component) |
!< q_aux(3 ): v  (y velocity component) | primitive variables
!< q_aux(4 ): w  (z velocity component) |
!< q_aux(5 ): P  (pressure)             |
!< q_aux(6 ): T  (temperature)
!< q_aux(7 ): E  (total energy)
!< q_aux(8 ): H  (total entalpy)
!< q_aux(9 ): a  (speed of sound)
!< q_aux(10): cv (specific heat at constant volume)
!< q_aux(11): Rf (fluid constant, cp-cv)
!< q_aux(12): g  (specific heats ratio, cp/cv)
!<```
!< Note that the first nv=5 auxiliary variables are also the primitive ones.

! ADAM modules
use :: adam_mpih_object, only : mpih_object
use :: adam_eos_ic_object, only : eos_ic_object
! CHASE modules
use :: adam_chase_parameters
! third party modules
use :: finer, only : file_ini
use :: penf, only : I4P, R8P, str

implicit none
private
public :: chase_physics_object
public :: compute_auxiliary
public :: compute_conservative
public :: compute_fluxes_conservative

character(len=7), parameter :: INI_SECTION_NAME='physics' !< INI file section name containing fluid physics.

integer(I4P),  parameter, public :: VAR_R =1_I4P                       !< Auxiliary/primitive/conservative variable 1, r, density.
integer(I4P),  parameter, public :: VAR_U =2_I4P                       !< Auxiliary/primitive variable 2, u, x velocity component.
integer(I4P),  parameter, public :: VAR_V =3_I4P                       !< Auxiliary/primitive variable 3, v, y velocity component.
integer(I4P),  parameter, public :: VAR_W =4_I4P                       !< Auxiliary/primitive variable 4, w, z velocity component.
integer(I4P),  parameter, public :: VAR_P =5_I4P                       !< Auxiliary/primitive variable 5, P, pressure.
integer(I4P),  parameter, public :: VAR_A =6_I4P                       !< Auxiliary           variable 6, a, speed of sound.
integer(I4P),  parameter, public :: VAR_T =7_I4P                       !< Auxiliary           variable 7, T, temperature.
integer(I4P),  parameter, public :: VAR_E =8_I4P                       !< Auxiliary           variable 8, E, total specific energy.
integer(I4P),  parameter, public :: VAR_H =9_I4P                       !< Auxiliary           variable 9, H, total specific entalpy.
integer(I4P),  parameter, public :: VAR_RF=10_I4P                      !< Auxiliary           variable 10, Rf, cp-cv.
integer(I4P),  parameter, public :: VAR_G =11_I4P                      !< Auxiliary           variable 11, g, cp/cv.
integer(I4P),  parameter, public :: VAR_RU=2_I4P                       !< Conservative variable 2, r*u, x momemtum component.
integer(I4P),  parameter, public :: VAR_RV=3_I4P                       !< Conservative variable 3, r*v, y momemtum component.
integer(I4P),  parameter, public :: VAR_RW=4_I4P                       !< Conservative variable 4, r*w, z momemtum component.
integer(I4P),  parameter, public :: VAR_RE=5_I4P                       !< Conservative variable 5, r*E, total energy.
character(12), parameter, public :: WENO_REC_VAR_CONS='CONSERVATIVE'   !< WENO reconstruction on conservative variables.
character(15), parameter, public :: WENO_REC_VAR_CHAR='CHARACTERISTICS'!< WENO reconstruction on characteristics variables.

type :: chase_physics_object
   !< CHASE physics class definition.
   type(mpih_object)         :: mpih               !< MPI handler.
   integer(I4P)              :: nv=5_I4P           !< Number of conservative/primitive variables.
   integer(I4P)              :: nv_aux=11_I4P      !< Number of auxiliary variables.
   character(:), allocatable :: weno_rec_var       !< Type of WENO reconstruction variables (cons., charct.,...).
   real(R8P), pointer        :: erw(:,:,:)=>null() !< Right eigenvectors for WENO reconstruction.
   real(R8P), pointer        :: elw(:,:,:)=>null() !< Left  eigenvectors for WENO reconstruction.
   type(eos_ic_object)       :: eos                !< Equations of state of ideal compressibile fluid.
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize physics.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype chase_physics_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(chase_physics_object), intent(in) :: self             !< Physics.
   character(len=:), allocatable           :: desc             !< Description.
   character(len=1), parameter             :: NL=new_line('a') !< New line character.

   desc =       self%mpih%myrankstr//'Physics main data:'                                                //NL
   desc = desc//self%mpih%myrankstr//'  number of conservative variables (nv):  '//trim(str(self%nv    ))//NL
   desc = desc//self%mpih%myrankstr//'  number of auxiliary variables (nv_aux): '//trim(str(self%nv_aux))//NL
   desc = desc//self%mpih%myrankstr//'  WENO reconstruction variables:          '//self%weno_rec_var
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize the equation.
   class(chase_physics_object), intent(inout) :: self            !< Physics.
   type(file_ini),              intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'chase_physics_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'chase_physics_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(chase_physics_object), intent(inout)        :: self            !< Physics.
   type(file_ini),              intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                     intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                           :: go_on_fail_     !< Go on if load fails.
   integer(I4P)                                      :: s               !< Counter.
   integer(I4P)                                      :: error           !< Error status.
   character(99)                                     :: buff            !< Character buffer.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='weno_rec_var', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(weno_rec_var)')
   select case(trim(adjustl(buff)))
   case('CONSERVATIVE', 'conservative', 'Conservative')
      self%weno_rec_var = WENO_REC_VAR_CONS
      self%erw => IERL
      self%elw => IERL
   case('CHARACTERISTICS', 'characteristics', 'Characteristics')
      self%weno_rec_var = WENO_REC_VAR_CHAR
      self%erw => ER
      self%elw => EL
   case default
      call self%mpih%print_message(msg='warning: WENO reconstruction variable "'//trim(adjustl(buff))// &
                                   '" unknown. Revert back to conservative variables WENO reconstruction')
      self%weno_rec_var = WENO_REC_VAR_CONS
      self%erw => ER
      self%elw => EL
   endselect
   call self%eos%initialize(file_parameters=file_parameters, s=1, verbose=.true.)
   endsubroutine load_from_file

   ! public non TBP
   pure subroutine compute_auxiliary(cv, Rf, g, conservative, auxiliary)
   !< Compute auxiliary variables from conservative ones.
   real(R8P), intent(in)    :: cv               !< Specific heat at constant volume.
   real(R8P), intent(in)    :: Rf               !< Fluid constant (cp-cv).
   real(R8P), intent(in)    :: g                !< Specific heats ration (cp/cv).
   real(R8P), intent(in)    :: conservative(1:) !< Conservative variables.
   real(R8P), intent(inout) :: auxiliary(1:)    !< Auxiliary variables.
   real(R8P)                :: r, u, v, w, e, t !< State variables.

   r = conservative(VAR_R )
   u = conservative(VAR_RU)/r
   v = conservative(VAR_RV)/r
   w = conservative(VAR_RW)/r
   e = conservative(VAR_RE)/r
   t = (e-0.5_R8P*(u*u+v*v+w*w))/cv
   auxiliary(VAR_R) = r            ! density
   auxiliary(VAR_U) = u            ! velocity x
   auxiliary(VAR_V) = v            ! velocity y
   auxiliary(VAR_W) = w            ! velocity z
   auxiliary(VAR_P) = Rf*r*t       ! pressure
   auxiliary(VAR_T) = t            ! temperature
   auxiliary(VAR_E) = e            ! total energy
   auxiliary(VAR_H) = e+Rf*t       ! total entalpy
   auxiliary(VAR_A) = sqrt(g*Rf*t) ! sound speed
   auxiliary(VAR_RF)= Rf           ! fluid constant, cp-cv
   auxiliary(VAR_G) = g            ! specific heats ration, cp/cv
   endsubroutine compute_auxiliary

   pure subroutine compute_conservative(auxiliary, conservative)
   !< Compute convervative variables from auxiliary ones, scalar input.
   real(R8P), intent(in)    :: auxiliary(1:)    !< Auxiliary variables.
   real(R8P), intent(inout) :: conservative(1:) !< Conservative variables.

   associate(r=>auxiliary(VAR_R),u=>auxiliary(VAR_U),v=>auxiliary(VAR_V),w=>auxiliary(VAR_W),&
             p=>auxiliary(VAR_P),H=>auxiliary(VAR_H))
   conservative(VAR_R ) = r
   conservative(VAR_RU) = r*u
   conservative(VAR_RV) = r*v
   conservative(VAR_RW) = r*w
   conservative(VAR_RE) = r*H - p
   endassociate
   endsubroutine compute_conservative

   pure subroutine compute_fluxes_conservative(sir,auxiliary,fluxes)
   !< Compute fluxes of convervative variables from auxiliary ones.
   real(R8P), intent(in)    :: sir(3)        !< Directional (1=x,2=y,3=z) increment.
   real(R8P), intent(in)    :: auxiliary(1:) !< Auxiliary variables.
   real(R8P), intent(inout) :: fluxes(1:)    !< Conservative fluxes.

   associate(r=>auxiliary(VAR_R),u=>auxiliary(VAR_U),v=>auxiliary(VAR_V),w=>auxiliary(VAR_W),&
             p=>auxiliary(VAR_P),H=>auxiliary(VAR_H))
   fluxes(VAR_R ) = r*u*sir(1) + r*v*sir(2) + r*w*sir(3)
   fluxes(VAR_RU) = fluxes(VAR_R)*u + p*sir(1)
   fluxes(VAR_RV) = fluxes(VAR_R)*v + p*sir(2)
   fluxes(VAR_RW) = fluxes(VAR_R)*w + p*sir(3)
   fluxes(VAR_RE) = fluxes(VAR_R)*H
   endassociate
   endsubroutine compute_fluxes_conservative
endmodule adam_chase_physics_object
