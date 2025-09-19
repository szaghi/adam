!< ADAM, PRISM (Plasma Research usIng Simulation Methods) physics class definition, common backend.
module adam_prism_physics_object
!< ADAM, PRISM (Plasma Research usIng Simulation Methods) physics class definition, common backend.
!<
!< Field variables are arranged as follows:
!<```
!< q(1): Dx
!< q(2): Dy
!< q(3): Dz
!< q(4): Bx
!< q(5): By
!< q(6): Bz
!< q(7): Jx
!< q(8): Jy
!< q(9): Jz
!<```
!< Field variables fluxes are:
!<```
!< Fx(1) =  0       Fy(1) = -Bz/muz  Fz(1) =  By/muy
!< Fx(2) =  Bz/muz  Fy(2) =  0       Fz(2) = -Bx/mux
!< Fx(3) = -By/muy  Fy(3) =  Bx/mux  Fz(3) =  0
!< Fx(4) =  0       Fy(4) =  Dz/epsz Fz(4) = -Dy/epsy
!< Fx(5) = -Dz/epsz Fy(5) =  0       Fz(5) =  Dx/epsx
!< Fx(6) =  Dy/epsy Fy(6) = -Dx/epsx Fz(6) =  0
!< Fx(7) =  0       Fy(7) =  0       Fz(7) =  0
!< Fx(8) =  0       Fy(8) =  0       Fz(8) =  0
!< Fx(9) =  0       Fy(9) =  0       Fz(9) =  0
!<```

! ADAM modules
use :: adam_mpih_object, only : mpih_object
! PRISM modules
use :: adam_prism_parameters
! third party modules
use :: finer, only : file_ini
use :: penf, only : I4P, R8P, str

implicit none
private
public :: prism_physics_object

character(len=7), parameter :: INI_SECTION_NAME='physics' !< INI file section name containing fluid physics.

integer(I4P),  parameter, public :: VAR_DX = 1_I4P                      !< Conservative variable 1, Dx.
integer(I4P),  parameter, public :: VAR_DY = 2_I4P                      !< Conservative variable 2, Dy.
integer(I4P),  parameter, public :: VAR_DZ = 3_I4P                      !< Conservative variable 3, Dz.
integer(I4P),  parameter, public :: VAR_BX = 4_I4P                      !< Conservative variable 4, Bx.
integer(I4P),  parameter, public :: VAR_BY = 5_I4P                      !< Conservative variable 5, By.
integer(I4P),  parameter, public :: VAR_BZ = 6_I4P                      !< Conservative variable 6, Bz.
integer(I4P),  parameter, public :: VAR_JX = 7_I4P                      !< Source variable 1, Jx.
integer(I4P),  parameter, public :: VAR_JY = 8_I4P                      !< Source variable 2, Jy.
integer(I4P),  parameter, public :: VAR_JZ = 9_I4P                      !< Source variable 3, Jz.
character(12), parameter, public :: WENO_REC_VAR_CONS='CONSERVATIVE'    !< WENO reconstruction on conservative variables.
character(15), parameter, public :: WENO_REC_VAR_CHAR='CHARACTERISTICS' !< WENO reconstruction on characteristics variables.
character(7),  parameter, public :: DIV_CORR_VAR_POISS='POISSON'        !< Poisson divergence correction.
character(10), parameter, public :: DIV_CORR_VAR_HYPER='HYPERBOLIC'     !< Hyperbolic divergence correction.

type :: prism_physics_object
   !< PRISM physics class definition.
   type(mpih_object)         :: mpih                         !< MPI handler.
   integer(I4P)              :: nv    = 9_I4P                !< Number of variables in q vector (nv=nv_c+nv_s+nv_cl).
   integer(I4P)              :: nv_c  = 6_I4P                !< Number of conservative variables in q vector.
   integer(I4P)              :: nv_s  = 3_I4P                !< Number of source variables in q vector.
   integer(I4P)              :: nv_cl = 0_I4P                !< Number of divergence cleaning variables in q vector.
   logical                   :: d_divergence_cleaner=.false. !< Enable electric field divergence cleaning.
   logical                   :: b_divergence_cleaner=.false. !< Enable magnetic field divergence cleaning.
   real(R8P)                 :: chi                          !< Coefficiente for D div-cleaning.
   real(R8P)                 :: eta                          !< Coefficiente for B div-cleaning.
   real(R8P)                 :: evmax                        !< Maximum signal speed (eigenvalue).
   character(:), allocatable :: weno_rec_var                 !< Type of WENO reconstruction variables (cons., charct.,...).
   character(:), allocatable :: div_corr_var                 !< Type of divergence correction variables (poisson, hyperbolic,...).
   real(R8P), pointer        :: erw(:,:,:)=>null()           !< Right eigenvectors for WENO reconstruction.
   real(R8P), pointer        :: elw(:,:,:)=>null()           !< Left  eigenvectors for WENO reconstruction.
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize physics.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype prism_physics_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_physics_object), intent(in) :: self             !< Physics.
   character(len=:), allocatable           :: desc             !< Description.
   character(len=1), parameter             :: NL=new_line('a') !< New line character.

   desc =       self%mpih%myrankstr//'Physics main data:'                                                                    //NL
   desc = desc//self%mpih%myrankstr//'  number of variables in q (nv):                '//trim(str(self%nv                  ))//NL
   desc = desc//self%mpih%myrankstr//'  number of conservative variables in q (nv_c): '//trim(str(self%nv_c                ))//NL
   desc = desc//self%mpih%myrankstr//'  WENO reconstruction variables:                '//self%weno_rec_var                   //NL
   desc = desc//self%mpih%myrankstr//'  Divergence correction:                        '//self%div_corr_var                   //NL    
   desc = desc//self%mpih%myrankstr//'  D divergence correction:                      '//trim(str(self%d_divergence_cleaner))//NL
   desc = desc//self%mpih%myrankstr//'  B divergence correction:                      '//trim(str(self%b_divergence_cleaner))//NL
   desc = desc//self%mpih%myrankstr//'  Chi:                                          '//trim(str(self%chi                 ))//NL
   desc = desc//self%mpih%myrankstr//'  Eta:                                          '//trim(str(self%eta                 ))//NL
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize the equation.
   class(prism_physics_object), intent(inout) :: self            !< Physics.
   type(file_ini),              intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'prism_physics_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'prism_physics_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(prism_physics_object), intent(inout)        :: self            !< Physics.
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

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='divergence_correction', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(divergence_correction)')
   select case(trim(adjustl(buff)))
   case('POISSON', 'poisson', 'Poisson')
      self%div_corr_var = DIV_CORR_VAR_POISS
   case('HYPERBOLIC', 'hyperbolic', 'Hyperbolic')
      self%div_corr_var = DIV_CORR_VAR_HYPER
   case default
      call self%mpih%print_message(msg='warning: divergence correction variable not activated')
      self%div_corr_var = 'No'
   endselect   

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='d_divergence_cleaner', &
                            val=self%d_divergence_cleaner,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(d_divergence_cleaner)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='b_divergence_cleaner', &
                            val=self%b_divergence_cleaner,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(b_divergence_cleaner)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='chi', val=self%chi, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(chi)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='eta', val=self%eta, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(eta)')
   ! increase variables number if q if cleaners are used
   if     (self%d_divergence_cleaner.and..not.self%b_divergence_cleaner.and.self%div_corr_var==DIV_CORR_VAR_HYPER)then
      self%nv = self%nv + 1_I4P
      self%nv_cl = 1_I4P
   elseif (self%d_divergence_cleaner.and.self%b_divergence_cleaner.and.self%div_corr_var == DIV_CORR_VAR_HYPER) then
      self%nv = self%nv + 2_I4P
      self%nv_cl = 2_I4P
   endif
   ! consistency check
   if     (.not.self%d_divergence_cleaner.and..not.self%b_divergence_cleaner.and.self%nv/=9_I4P) then
      call self%mpih%error_stop(msg=': D and B divergence cleaners are false but nv /= 9')
   elseif (.not.self%d_divergence_cleaner.and.self%b_divergence_cleaner) then
      call self%mpih%error_stop(msg=': D divergence cleaner is false but B divergence cleaner is true')
   elseif (self%d_divergence_cleaner.and.self%chi<=1._R8P) then
      call self%mpih%error_stop(msg=': if D divergence cleaner is true chi cannot be lower than 1.0')
   endif
   self%evmax = sqrt(1._R8P/(EPS0*MU0))
   if (self%d_divergence_cleaner .and. self%div_corr_var == DIV_CORR_VAR_HYPER) then 
      self%evmax = self%chi*sqrt(1._R8P/(EPS0*MU0))
   endif
   endsubroutine load_from_file
endmodule adam_prism_physics_object
