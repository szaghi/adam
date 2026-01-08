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
use :: adam_prism_numerics_object, only : RECONSTRUCTION_VARS_CONS, RECONSTRUCTION_VARS_CHAR, &
                                          DIV_CORR_VAR_HYPER, DIV_CORR_VAR_POISS
use :: adam_prism_parameters
! third party modules
use :: finer, only : file_ini
use :: penf, only : I4P, R8P, str

implicit none
private
public :: prism_physics_object
public :: EM_PHYSICAL_MODEL
public :: PIC_PHYSICAL_MODEL

character(len=7),  parameter :: INI_SECTION_NAME=   'physics'         !< INI file section name containing fluid physics.
character(len=15), parameter :: EM_PHYSICAL_MODEL=  'electromagnetic' !< Electromagnetic physical model.
character(len=3),  parameter :: PIC_PHYSICAL_MODEL= 'PIC'             !< PIC physical model.

integer(I4P),  parameter, public :: VAR_DX = 1_I4P !< Conservative variable 1, Dx.
integer(I4P),  parameter, public :: VAR_DY = 2_I4P !< Conservative variable 2, Dy.
integer(I4P),  parameter, public :: VAR_DZ = 3_I4P !< Conservative variable 3, Dz.
integer(I4P),  parameter, public :: VAR_BX = 4_I4P !< Conservative variable 4, Bx.
integer(I4P),  parameter, public :: VAR_BY = 5_I4P !< Conservative variable 5, By.
integer(I4P),  parameter, public :: VAR_BZ = 6_I4P !< Conservative variable 6, Bz.
!integer(I4P),  parameter, public :: VAR_JX = 7_I4P !< Source variable 1, Jx.
!integer(I4P),  parameter, public :: VAR_JY = 8_I4P !< Source variable 2, Jy.
!integer(I4P),  parameter, public :: VAR_JZ = 9_I4P !< Source variable 3, Jz.

type :: prism_physics_object
   !< PRISM physics class definition.
   type(mpih_object)           :: mpih                   !< MPI handler.
   character(len=99)           :: physical_model         !< Physical model.
   integer(I4P)                :: nv    = 9_I4P          !< Number of variables in q vector (nv=nv_c+nv_s+nv_cl).
   integer(I4P)                :: nv_c  = 6_I4P          !< Number of conservative variables in q vector.
   integer(I4P)                :: nv_s  = 3_I4P          !< Number of source variables in q vector.
   integer(I4P)                :: nv_cl = 0_I4P          !< Number of divergence cleaning variables in q vector.
   integer(I4P)                :: nv_pic = 0_I4P         !< Number of PIC variables in q vector.
   integer(I4P)                :: var_Jx, var_Jy, var_Jz !< Indices of current density components in q vector.
   real(R8P)                   :: chi                    !< Speed coefficient for D & B div-cleaning.
   !real(R8P)                   :: eta                   !< Coefficiente for B div-cleaning.
   real(R8P)                   :: evmax                  !< Maximum signal speed (eigenvalue).
   real(R8P), pointer          :: erw(:,:,:)=>null()     !< Right eigenvectors for high order reconstruction.
   real(R8P), pointer          :: elw(:,:,:)=>null()     !< Left  eigenvectors for high order reconstruction.
   real(R8P)                   :: EV_D(7)                !< Eigenvalues with D divergence cleaning.
   real(R8P)                   :: EV_B(7)                !< Eigenvalues with B divergence cleaning.
   real(R8P)                   :: EV_D_B(8)              !< Eigenvalues with D and B divergence cleaning.
   real(R8P)                   :: ER_D(7,7,3)            !< Right eigenvectors with D divergence cleaning.
   real(R8P)                   :: ER_B(7,7,3)            !< Right eigenvectors with B divergence cleaning.
   real(R8P)                   :: ER_D_B(8,8,3)          !< Right eigenvectors with D and B divergence cleaning.
   real(R8P)                   :: EL_D(7,7,3)            !< Left eigenvectors with D divergence cleaning.
   real(R8P)                   :: EL_B(7,7,3)            !< Left eigenvectors with B divergence cleaning.
   real(R8P)                   :: EL_D_B(8,8,3)          !< Left eigenvectors with D and B divergence cleaning.

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

   desc =       self%mpih%myrankstr//'Physics main data:'                                                         //NL
   desc = desc//self%mpih%myrankstr//'  physical model:                               '//trim(self%physical_model)//NL
   desc = desc//self%mpih%myrankstr//'  number of variables in q (nv):                '//trim(str(self%nv       ))//NL
   desc = desc//self%mpih%myrankstr//'  number of conservative variables in q (nv_c): '//trim(str(self%nv_c     ))//NL
   desc = desc//self%mpih%myrankstr//'  number of PIC variables in q (nv_PIC):        '//trim(str(self%nv_PIC   ))//NL
   desc = desc//self%mpih%myrankstr//'  Chi:                                          '//trim(str(self%chi      ))!//NL
   !desc = desc//self%mpih%myrankstr//'  Eta:                                          '//trim(str(self%eta ))
   endfunction description

   subroutine initialize(self, file_parameters, reconstruction_vars, div_corr_var, &
                         constrained_transport_D, constrained_transport_B)
   !< Initialize the equation.
   class(prism_physics_object), target, intent(inout)        :: self                    !< Physics.
   type(file_ini),                      intent(in)           :: file_parameters         !< Simulation parameters ini file handler.
   character(*),                        intent(in), optional :: reconstruction_vars     !< Variables used for HO reconstruction.
   character(*),                        intent(in), optional :: div_corr_var            !< Type of divergence correction variables
                                                                                        !<  (poisson, hyperbolic,...).
   logical,                             intent(in), optional :: constrained_transport_D !< Enable CT-Correction on D.
   logical,                             intent(in), optional :: constrained_transport_B !< Enable CT-Correction on B.
   character(:), allocatable                                 :: reconstruction_vars_    !< Variables used for HO reconstruction.
   character(:), allocatable                                 :: div_corr_var_           !< Type of divergence correction variables.
   real(R8P)                                                 :: ch                      !< Divergence cleaning propagation speed.

   call self%load_from_file(file_parameters=file_parameters, div_corr_var=div_corr_var)
   ch = self%evmax

   self%EV_D =   [0._R8P,ch,-ch,C0,C0,-C0,-C0]
   self%EV_B =   [0._R8P,ch,-ch,C0,C0,-C0,-C0]
   self%EV_D_B = [ch,ch,-ch,-ch,C0,C0,-C0,-C0]

   self%ER_D   = reshape([                                                                                        &
                 0._R8P, 1._R8P/ch, -1._R8P/ch, 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         , &
                 0._R8P, 0._R8P   , 0._R8P    , 0._R8P         , sqrt(EPS0/MU0), 0._R8P        , -sqrt(EPS0/MU0), &
                 0._R8P, 0._R8P   , 0._R8P    , -sqrt(EPS0/MU0), 0._R8P        , sqrt(EPS0/MU0), 0._R8P         , &
                 1._R8P, 0._R8P   , 0._R8P    , 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         , &
                 0._R8P, 0._R8P   , 0._R8P    , 1._R8P         , 0._R8P        , 1._R8P        , 0._R8P         , &
                 0._R8P, 0._R8P   , 0._R8P    , 0._R8P         , 1._R8P        , 0._R8P        , 1._R8P         , &
                 0._R8P, 1._R8P   , 1._R8P    , 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         , &

                 0._R8P, 0._R8P   , 0._R8P    , 0._R8P        , -sqrt(EPS0/MU0), 0._R8P         , sqrt(EPS0/MU0), &
                 0._R8P, 1._R8P/ch, -1._R8P/ch, 0._R8P        , 0._R8P         , 0._R8P         , 0._R8P        , &
                 0._R8P, 0._R8P   , 0._R8P    , sqrt(EPS0/MU0), 0._R8P         , -sqrt(EPS0/MU0), 0._R8P        , &
                 0._R8P, 0._R8P   , 0._R8P    , 1._R8P        , 0._R8P         , 1._R8P         , 0._R8P        , &
                 1._R8P, 0._R8P   , 0._R8P    , 0._R8P        , 0._R8P         , 0._R8P         , 0._R8P        , &
                 0._R8P, 0._R8P   , 0._R8P    , 0._R8P        , 1._R8P         , 0._R8P         , 1._R8P        , &
                 0._R8P, 1._R8P   , 1._R8P    , 0._R8P        , 0._R8P         , 0._R8P         , 0._R8P        , &

                 0._R8P, 0._R8P   , 0._R8P    , 0._R8P         , sqrt(EPS0/MU0), 0._R8P        , -sqrt(EPS0/MU0), &
                 0._R8P, 0._R8P   , 0._R8P    , -sqrt(EPS0/MU0), 0._R8P        , sqrt(EPS0/MU0), 0._R8P         , &
                 0._R8P, 1._R8P/ch, -1._R8P/ch, 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         , &
                 0._R8P, 0._R8P   , 0._R8P    , 1._R8P         , 0._R8P        , 1._R8P        , 0._R8P         , &
                 0._R8P, 0._R8P   , 0._R8P    , 0._R8P         , 1._R8P        , 0._R8P        , 1._R8P         , &
                 1._R8P, 0._R8P   , 0._R8P    , 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         , &
                 0._R8P, 1._R8P   , 1._R8P    , 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         ] &
                 , [7,7,3])

   self%ER_B   = reshape([                                                                                        &
                 1._R8P, 0._R8P   , 0._R8P    , 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         , &
                 0._R8P, 0._R8P   , 0._R8P    , 0._R8P         , sqrt(EPS0/MU0), 0._R8P        , -sqrt(EPS0/MU0), &
                 0._R8P, 0._R8P   , 0._R8P    , -sqrt(EPS0/MU0), 0._R8P        , sqrt(EPS0/MU0), 0._R8P         , &
                 0._R8P, 1._R8P/ch, -1._R8P/ch, 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         , &
                 0._R8P, 0._R8P   , 0._R8P    , 1._R8P         , 0._R8P        , 1._R8P        , 0._R8P         , &
                 0._R8P, 0._R8P   , 0._R8P    , 0._R8P         , 1._R8P        , 0._R8P        , 1._R8P         , &
                 0._R8P, 1._R8P   , 1._R8P    , 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         , &

                 0._R8P, 0._R8P   , 0._R8P    , 0._R8P        , -sqrt(EPS0/MU0), 0._R8P         , sqrt(EPS0/MU0), &
                 1._R8P, 0._R8P   , 0._R8P    , 0._R8P        , 0._R8P         , 0._R8P         , 0._R8P        , &
                 0._R8P, 0._R8P   , 0._R8P    , sqrt(EPS0/MU0), 0._R8P         , -sqrt(EPS0/MU0), 0._R8P        , &
                 0._R8P, 0._R8P   , 0._R8P    , 1._R8P        , 0._R8P         , 1._R8P         , 0._R8P        , &
                 0._R8P, 1._R8P/ch, -1._R8P/ch, 0._R8P        , 0._R8P         , 0._R8P         , 0._R8P        , &
                 0._R8P, 0._R8P   , 0._R8P    , 0._R8P        , 1._R8P         , 0._R8P         , 1._R8P        , &
                 0._R8P, 1._R8P   , 1._R8P    , 0._R8P        , 0._R8P         , 0._R8P         , 0._R8P        , &

                 0._R8P, 0._R8P   , 0._R8P    , 0._R8P         , sqrt(EPS0/MU0), 0._R8P        , -sqrt(EPS0/MU0), &
                 0._R8P, 0._R8P   , 0._R8P    , -sqrt(EPS0/MU0), 0._R8P        , sqrt(EPS0/MU0), 0._R8P         , &
                 1._R8P, 0._R8P   , 0._R8P    , 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         , &
                 0._R8P, 0._R8P   , 0._R8P    , 1._R8P         , 0._R8P        , 1._R8P        , 0._R8P         , &
                 0._R8P, 0._R8P   , 0._R8P    , 0._R8P         , 1._R8P        , 0._R8P        , 1._R8P         , &
                 0._R8P, 1._R8P/ch, -1._R8P/ch, 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         , &
                 0._R8P, 1._R8P   , 1._R8P    , 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         ] &
                 , [7,7,3])

   self%ER_D_B = reshape([                                                                                                       &
                 1._R8P/ch, 0._R8P,    -1._R8P/ch, 0._R8P    , 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         , &
                 0._R8P   , 0._R8P,    0._R8P    , 0._R8P    , 0._R8P         , sqrt(EPS0/MU0), 0._R8P        , -sqrt(EPS0/MU0), &
                 0._R8P   , 0._R8P,    0._R8P    , 0._R8P    , -sqrt(EPS0/MU0), 0._R8P        , sqrt(EPS0/MU0), 0._R8P         , &
                 0._R8P   , 1._R8P/ch, 0._R8P    , -1._R8P/ch, 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         , &
                 0._R8P   , 0._R8P,    0._R8P    , 0._R8P    , 1._R8P         , 0._R8P        , 1._R8P        , 0._R8P         , &
                 0._R8P   , 0._R8P,    0._R8P    , 0._R8P    , 0._R8P         , 1._R8P        , 0._R8P        , 1._R8P         , &
                 1._R8P   , 0._R8P,    1._R8P    , 0._R8P    , 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         , &
                 0._R8P   , 1._R8P,    0._R8P    , 1._R8P    , 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         , &
                 !
                 0._R8P   , 0._R8P,    0._R8P    , 0._R8P    , 0._R8P        , -sqrt(EPS0/MU0), 0._R8P         , sqrt(EPS0/MU0), &
                 1._R8P/ch, 0._R8P,    -1._R8P/ch, 0._R8P    , 0._R8P        , 0._R8P         , 0._R8P         , 0._R8P        , &
                 0._R8P   , 0._R8P,    0._R8P    , 0._R8P    , sqrt(EPS0/MU0), 0._R8P         , -sqrt(EPS0/MU0), 0._R8P        , &
                 0._R8P   , 0._R8P,    0._R8P    , 0._R8P    , 1._R8P        , 0._R8P         , 1._R8P         , 0._R8P        , &
                 0._R8P   , 1._R8P/ch, 0._R8P    , -1._R8P/ch, 0._R8P        , 0._R8P         , 0._R8P         , 0._R8P        , &
                 0._R8P   , 0._R8P,    0._R8P    , 0._R8P    , 0._R8P        , 1._R8P         , 0._R8P         , 1._R8P        , &
                 1._R8P   , 0._R8P,    1._R8P    , 0._R8P    , 0._R8P        , 0._R8P         , 0._R8P         , 0._R8P        , &
                 0._R8P   , 1._R8P,    0._R8P    , 1._R8P    , 0._R8P        , 0._R8P         , 0._R8P         , 0._R8P        , &
                 !
                 0._R8P   , 0._R8P   , 0._R8P    , 0._R8P    , 0._R8P         , sqrt(EPS0/MU0), 0._R8P        , -sqrt(EPS0/MU0), &
                 0._R8P   , 0._R8P   , 0._R8P    , 0._R8P    , -sqrt(EPS0/MU0), 0._R8P        , sqrt(EPS0/MU0), 0._R8P         , &
                 1._R8P/ch, 0._R8P   , -1._R8P/ch, 0._R8P    , 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         , &
                 0._R8P   , 0._R8P   , 0._R8P    , 0._R8P    , 1._R8P         , 0._R8P        , 1._R8P        , 0._R8P         , &
                 0._R8P   , 0._R8P   , 0._R8P    , 0._R8P    , 0._R8P         , 1._R8P        , 0._R8P        , 1._R8P         , &
                 0._R8P   , 1._R8P/ch, 0._R8P    , -1._R8P/ch, 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         , &
                 1._R8P   , 0._R8P   , 1._R8P    , 0._R8P    , 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         , &
                 0._R8P   , 1._R8P   , 0._R8P    , 1._R8P    , 0._R8P         , 0._R8P        , 0._R8P        , 0._R8P         ] &
                 , [8,8,3]) ! < Right eigenvectors, D and B divergence cleaning.

   self%EL_D   = reshape([                                                                                                        &
                 0._R8P    , 0._R8P                , 0._R8P                , 1._R8P, 0._R8P       , 0._R8P       , 0._R8P       , &
                 ch/2._R8P , 0._R8P                , 0._R8P                , 0._R8P, 0._R8P       , 0._R8P       , 1._R8P/2._R8P, &
                 -ch/2._R8P, 0._R8P                , 0._R8P                , 0._R8P, 0._R8P       , 0._R8P       , 1._R8P/2._R8P, &
                 0._R8P    , 0._R8P                , -sqrt(MU0/EPS0)/2._R8P, 0._R8P, 1._R8P/2._R8P, 0._R8P       , 0._R8P       , &
                 0._R8P    , sqrt(MU0/EPS0)/2._R8P , 0._R8P                , 0._R8P, 0._R8P       , 1._R8P/2._R8P, 0._R8P       , &
                 0._R8P    , 0._R8P                , sqrt(MU0/EPS0)/2._R8P , 0._R8P, 1._R8P/2._R8P, 0._R8P       , 0._R8P       , &
                 0._R8P    , -sqrt(MU0/EPS0)/2._R8P, 0._R8P                , 0._R8P, 0._R8P       , 1._R8P/2._R8P, 0._R8P       , &

                 0._R8P                , 0._R8P    , 0._R8P                , 0._R8P       , 1._R8P, 0._R8P       , 0._R8P       , &
                 0._R8P                , ch/2._R8P , 0._R8P                , 0._R8P       , 0._R8P, 0._R8P       , 1._R8P/2._R8P, &
                 0._R8P                , -ch/2._R8P, 0._R8P                , 0._R8P       , 0._R8P, 0._R8P       , 1._R8P/2._R8P, &
                 0._R8P                , 0._R8P    , sqrt(MU0/EPS0)/2._R8P , 1._R8P/2._R8P, 0._R8P, 0._R8P       , 0._R8P       , &
                 -sqrt(MU0/EPS0)/2._R8P, 0._R8P    , 0._R8P                , 0._R8P       , 0._R8P, 1._R8P/2._R8P, 0._R8P       , &
                 0._R8P                , 0._R8P    , -sqrt(MU0/EPS0)/2._R8P, 1._R8P/2._R8P, 0._R8P, 0._R8P       , 0._R8P       , &
                 sqrt(MU0/EPS0)/2._R8P , 0._R8P    , 0._R8P                , 0._R8P       , 0._R8P, 1._R8P/2._R8P, 0._R8P       , &

                 0._R8P                , 0._R8P                , 0._R8P    , 0._R8P       , 0._R8P       , 1._R8P, 0._R8P       , &
                 0._R8P                , 0._R8P                , ch/2._R8P , 0._R8P       , 0._R8P       , 0._R8P, 1._R8P/2._R8P, &
                 0._R8P                , 0._R8P                , -ch/2._R8P, 0._R8P       , 0._R8P       , 0._R8P, 1._R8P/2._R8P, &
                 0._R8P                , -sqrt(MU0/EPS0)/2._R8P, 0._R8P    , 1._R8P/2._R8P, 0._R8P       , 0._R8P, 0._R8P       , &
                 sqrt(MU0/EPS0)/2._R8P , 0._R8P                , 0._R8P    , 0._R8P       , 1._R8P/2._R8P, 0._R8P, 0._R8P       , &
                 0._R8P                , sqrt(MU0/EPS0)/2._R8P , 0._R8P    , 1._R8P/2._R8P, 0._R8P       , 0._R8P, 0._R8P       , &
                 -sqrt(MU0/EPS0)/2._R8P, 0._R8P                , 0._R8P    , 0._R8P       , 1._R8P/2._R8P, 0._R8P, 0._R8P       ] &
                 , [7,7,3])

   self%EL_B   = reshape([                                                                                                        &
                 1._R8P, 0._R8P                , 0._R8P                , 0._R8P    , 0._R8P       , 0._R8P       , 0._R8P       , &
                 0._R8P, 0._R8P                , 0._R8P                , ch/2._R8P , 0._R8P       , 0._R8P       , 1._R8P/2._R8P, &
                 0._R8P, 0._R8P                , 0._R8P                , -ch/2._R8P, 0._R8P       , 0._R8P       , 1._R8P/2._R8P, &
                 0._R8P, 0._R8P                , -sqrt(MU0/EPS0)/2._R8P, 0._R8P    , 1._R8P/2._R8P, 0._R8P       , 0._R8P       , &
                 0._R8P, sqrt(MU0/EPS0)/2._R8P , 0._R8P                , 0._R8P    , 0._R8P       , 1._R8P/2._R8P, 0._R8P       , &
                 0._R8P, 0._R8P                , sqrt(MU0/EPS0)/2._R8P , 0._R8P    , 1._R8P/2._R8P, 0._R8P       , 0._R8P       , &
                 0._R8P, -sqrt(MU0/EPS0)/2._R8P, 0._R8P                , 0._R8P    , 0._R8P       , 1._R8P/2._R8P, 0._R8P       , &

                 0._R8P                , 1._R8P, 0._R8P                , 0._R8P       , 0._R8P    , 0._R8P       , 0._R8P       , &
                 0._R8P                , 0._R8P, 0._R8P                , 0._R8P       , ch/2._R8P , 0._R8P       , 1._R8P/2._R8P, &
                 0._R8P                , 0._R8P, 0._R8P                , 0._R8P       , -ch/2._R8P, 0._R8P       , 1._R8P/2._R8P, &
                 0._R8P                , 0._R8P, sqrt(MU0/EPS0)/2._R8P , 1._R8P/2._R8P, 0._R8P    , 0._R8P       , 0._R8P       , &
                 -sqrt(MU0/EPS0)/2._R8P, 0._R8P, 0._R8P                , 0._R8P       , 0._R8P    , 1._R8P/2._R8P, 0._R8P       , &
                 0._R8P                , 0._R8P, -sqrt(MU0/EPS0)/2._R8P, 1._R8P/2._R8P, 0._R8P    , 0._R8P       , 0._R8P       , &
                 sqrt(MU0/EPS0)/2._R8P , 0._R8P, 0._R8P                , 0._R8P       , 0._R8P    , 1._R8P/2._R8P, 0._R8P       , &

                 0._R8P                , 0._R8P                , 1._R8P, 0._R8P       , 0._R8P       , 0._R8P    , 0._R8P       , &
                 0._R8P                , 0._R8P                , 0._R8P, 0._R8P       , 0._R8P       , ch/2._R8P , 1._R8P/2._R8P, &
                 0._R8P                , 0._R8P                , 0._R8P, 0._R8P       , 0._R8P       , -ch/2._R8P, 1._R8P/2._R8P, &
                 0._R8P                , -sqrt(MU0/EPS0)/2._R8P, 0._R8P, 1._R8P/2._R8P, 0._R8P       , 0._R8P    , 0._R8P       , &
                 sqrt(MU0/EPS0)/2._R8P , 0._R8P                , 0._R8P, 0._R8P       , 1._R8P/2._R8P, 0._R8P    , 0._R8P       , &
                 0._R8P                , sqrt(MU0/EPS0)/2._R8P , 0._R8P, 1._R8P/2._R8P, 0._R8P       , 0._R8P    , 0._R8P       , &
                 -sqrt(MU0/EPS0)/2._R8P, 0._R8P                , 0._R8P, 0._R8P       , 1._R8P/2._R8P, 0._R8P    , 0._R8P       ] &
                 , [7,7,3])

   self%EL_D_B = reshape([                                                                                                        &
      ch/2._R8P ,0._R8P                ,0._R8P                ,0._R8P    ,0._R8P       ,0._R8P       ,1._R8P/2._R8P,0._R8P       ,&
      0._R8P    ,0._R8P                ,0._R8P                ,ch/2._R8P ,0._R8P       ,0._R8P       ,0._R8P       ,1._R8P/2._R8P,&
      -ch/2._R8P,0._R8P                ,0._R8P                ,0._R8P    ,0._R8P       ,0._R8P       ,1._R8P/2._R8P,0._R8P       ,&
      0._R8P    ,0._R8P                ,0._R8P                ,-ch/2._R8P,0._R8P       ,0._R8P       ,0._R8P       ,1._R8P/2._R8P,&
      0._R8P    ,0._R8P                ,-sqrt(MU0/EPS0)/2._R8P,0._R8P    ,1._R8P/2._R8P,0._R8P       ,0._R8P       ,0._R8P       ,&
      0._R8P    ,sqrt(MU0/EPS0)/2._R8P ,0._R8P                ,0._R8P    ,0._R8P       ,1._R8P/2._R8P,0._R8P       ,0._R8P       ,&
      0._R8P    ,0._R8P                ,sqrt(MU0/EPS0)/2._R8P ,0._R8P    ,1._R8P/2._R8P,0._R8P       ,0._R8P       ,0._R8P       ,&
      0._R8P    ,-sqrt(MU0/EPS0)/2._R8P,0._R8P                ,0._R8P    ,0._R8P       ,1._R8P/2._R8P,0._R8P       ,0._R8P       ,&

      0._R8P                ,ch/2._R8P ,0._R8P                ,0._R8P       ,0._R8P    ,0._R8P       ,1._R8P/2._R8P,0._R8P       ,&
      0._R8P                ,0._R8P    ,0._R8P                ,0._R8P       ,ch/2._R8P ,0._R8P       ,0._R8P       ,1._R8P/2._R8P,&
      0._R8P                ,-ch/2._R8P,0._R8P                ,0._R8P       ,0._R8P    ,0._R8P       ,1._R8P/2._R8P,0._R8P       ,&
      0._R8P                ,0._R8P    ,0._R8P                ,0._R8P       ,-ch/2._R8P,0._R8P       ,0._R8P       ,1._R8P/2._R8P,&
      0._R8P                ,0._R8P    ,sqrt(MU0/EPS0)/2._R8P ,1._R8P/2._R8P,0._R8P    ,0._R8P       ,0._R8P       ,0._R8P       ,&
      -sqrt(MU0/EPS0)/2._R8P,0._R8P    ,0._R8P                ,0._R8P       ,0._R8P    ,1._R8P/2._R8P,0._R8P       ,0._R8P       ,&
      0._R8P                ,0._R8P    ,-sqrt(MU0/EPS0)/2._R8P,1._R8P/2._R8P,0._R8P    ,0._R8P       ,0._R8P       ,0._R8P       ,&
      sqrt(MU0/EPS0)/2._R8P ,0._R8P    ,0._R8P                ,0._R8P       ,0._R8P    ,1._R8P/2._R8P,0._R8P       ,0._R8P       ,&

      0._R8P                ,0._R8P                ,ch/2._R8P ,0._R8P       ,0._R8P       ,0._R8P    ,1._R8P/2._R8P,0._R8P       ,&
      0._R8P                ,0._R8P                ,0._R8P    ,0._R8P       ,0._R8P       ,ch/2._R8P ,0._R8P       ,1._R8P/2._R8P,&
      0._R8P                ,0._R8P                ,-ch/2._R8P,0._R8P       ,0._R8P       ,0._R8P    ,1._R8P/2._R8P,0._R8P       ,&
      0._R8P                ,0._R8P                ,0._R8P    ,0._R8P       ,0._R8P       ,-ch/2._R8P,0._R8P       ,1._R8P/2._R8P,&
      0._R8P                ,-sqrt(MU0/EPS0)/2._R8P,0._R8P    ,1._R8P/2._R8P,0._R8P       ,0._R8P    ,0._R8P       ,0._R8P       ,&
      sqrt(MU0/EPS0)/2._R8P ,0._R8P                ,0._R8P    ,0._R8P       ,1._R8P/2._R8P,0._R8P    ,0._R8P       ,0._R8P       ,&
      0._R8P                , sqrt(MU0/EPS0)/2._R8P,0._R8P    ,1._R8P/2._R8P,0._R8P       ,0._R8P    ,0._R8P       ,0._R8P       ,&
      -sqrt(MU0/EPS0)/2._R8P,0._R8P                ,0._R8P    ,0._R8P       ,1._R8P/2._R8P,0._R8P    ,0._R8P       ,0._R8P       ]&
      , [8,8,3])

   reconstruction_vars_ = RECONSTRUCTION_VARS_CONS
   if (present(reconstruction_vars)) reconstruction_vars_ = trim(adjustl(reconstruction_vars))
   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'prism_physics_object%initialize start'

   div_corr_var_ = 'No'
   if (present(div_corr_var)) div_corr_var_ = trim(adjustl(div_corr_var))
   select case(div_corr_var_)
   case(DIV_CORR_VAR_POISS)
      self%nv_cl = 0_I4P
      self%var_Jx = 7_I4P
      self%var_Jy = 8_I4P
      self%var_Jz = 9_I4P
      select case(reconstruction_vars_)
         case(RECONSTRUCTION_VARS_CONS)
            self%erw => IERL
            self%elw => IERL
         case(RECONSTRUCTION_VARS_CHAR)
            self%erw => ER
            self%elw => EL
         case default
            call self%mpih%print_message(msg='warning: HO reconstruction variable "'//reconstruction_vars_// &
                                         '" unknown. Revert back to conservative variables reconstruction')
            self%erw => IERL
            self%elw => IERL
      endselect
   case(DIV_CORR_VAR_HYPER)
      if (constrained_transport_D .and. .not.constrained_transport_B) then
         self%nv_cl  = 1_I4P
         self%var_Jx = 8_I4P
         self%var_Jy = 9_I4P
         self%var_Jz = 10_I4P
         select case(reconstruction_vars_)
            case(RECONSTRUCTION_VARS_CONS)
               self%erw => IERL_D
               self%elw => IERL_D
            case(RECONSTRUCTION_VARS_CHAR)
               self%erw => self%ER_D
               self%elw => self%EL_D
            case default
               call self%mpih%print_message(msg='warning: HO reconstruction variable "'//reconstruction_vars_// &
                                            '" unknown. Revert back to conservative variables reconstruction')
               self%erw => IERL_D
               self%elw => IERL_D
            endselect
      elseif (.not.constrained_transport_D .and. constrained_transport_B) then
         self%nv_cl  = 1_I4P
         self%var_Jx = 8_I4P
         self%var_Jy = 9_I4P
         self%var_Jz = 10_I4P
         select case(reconstruction_vars_)
            case(RECONSTRUCTION_VARS_CONS)
               self%erw => IERL_B
               self%elw => IERL_B
            case(RECONSTRUCTION_VARS_CHAR)
               self%erw => self%ER_B
               self%elw => self%EL_B
            case default
               call self%mpih%print_message(msg='warning: HO reconstruction variable "'//reconstruction_vars_// &
                                            '" unknown. Revert back to conservative variables reconstruction')
               self%erw => IERL_B
               self%elw => IERL_B
         endselect
      elseif (constrained_transport_D .and. constrained_transport_B) then
         self%nv_cl  = 2_I4P
         self%var_Jx = 9_I4P
         self%var_Jy = 10_I4P
         self%var_Jz = 11_I4P
         select case(reconstruction_vars_)
            case(RECONSTRUCTION_VARS_CONS)
               self%erw => IERL_D_B
               self%elw => IERL_D_B
            case(RECONSTRUCTION_VARS_CHAR)
               self%erw => self%ER_D_B
               self%elw => self%EL_D_B
            case default
               call self%mpih%print_message(msg='warning: HO reconstruction variable "'//reconstruction_vars_// &
                                            '" unknown. Revert back to conservative variables reconstruction')
               self%erw => IERL_D_B
               self%elw => IERL_D_B
         endselect
      end if
   case default
      self%nv_cl  = 0_I4P
      self%var_Jx = 7_I4P
      self%var_Jy = 8_I4P
      self%var_Jz = 9_I4P
      select case(reconstruction_vars_)
      case(RECONSTRUCTION_VARS_CONS)
         self%erw => IERL
         self%elw => IERL
      case(RECONSTRUCTION_VARS_CHAR)
         self%erw => ER
         self%elw => EL
      case default
         call self%mpih%print_message(msg='warning: HO reconstruction variable "'//reconstruction_vars_// &
                                      '" unknown. Revert back to conservative variables reconstruction')
         self%erw => IERL
         self%elw => IERL
      endselect
   endselect
   self%nv = self%nv_c + self%nv_s + self%nv_cl
   self%nv_c = self%nv_c + self%nv_cl
   if (self%physical_model == PIC_PHYSICAL_MODEL) then
      self%nv_pic = 1_I4P
      self%nv = self%nv + self%nv_pic ! Per il PIC aggiungo nel vettore di stato la densita' di carica
   elseif (self%physical_model == EM_PHYSICAL_MODEL) then
      self%nv_pic = 0_I4P
   end if

   !print*, 'cazzo'
   !print*, self%erw(1,:,1)
   !print*, self%erw(2,:,1)
   !print*, self%erw(3,:,1)
   !print*, self%erw(4,:,1)
   !print*, self%erw(5,:,1)
   !print*, self%erw(6,:,1)
   !print*, self%erw(7,:,1)
   !print*, self%erw(1,:,2)
   !print*, self%erw(2,:,2)
   !print*, self%erw(3,:,2)
   !print*, self%erw(4,:,2)
   !print*, self%erw(5,:,2)
   !print*, self%erw(6,:,2)
   !print*, self%erw(7,:,2)
   !print*, self%erw(1,:,3)
   !print*, self%erw(2,:,3)
   !print*, self%erw(3,:,3)
   !print*, self%erw(4,:,3)
   !print*, self%erw(5,:,3)
   !print*, self%erw(6,:,3)
   !print*, self%erw(7,:,3)
!
   !print*, self%elw(1,:,1)
   !print*, self%elw(2,:,1)
   !print*, self%elw(3,:,1)
   !print*, self%elw(4,:,1)
   !print*, self%elw(5,:,1)
   !print*, self%elw(6,:,1)
   !print*, self%elw(7,:,1)
   !print*, self%elw(1,:,2)
   !print*, self%elw(2,:,2)
   !print*, self%elw(3,:,2)
   !print*, self%elw(4,:,2)
   !print*, self%elw(5,:,2)
   !print*, self%elw(6,:,2)
   !print*, self%elw(7,:,2)
   !print*, self%elw(1,:,3)
   !print*, self%elw(2,:,3)
   !print*, self%elw(3,:,3)
   !print*, self%elw(4,:,3)
   !print*, self%elw(5,:,3)
   !print*, self%elw(6,:,3)
   !print*, self%elw(7,:,3)

   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'prism_physics_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail, div_corr_var)
   !< Load config from file.
   class(prism_physics_object), intent(inout)        :: self            !< Physics.
   type(file_ini),              intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                     intent(in), optional :: go_on_fail      !< Go on if load fails.
   character(*),                intent(in), optional :: div_corr_var    !< Type of divergence correction variables (poisson, hyperbolic,...).
   logical                                           :: go_on_fail_     !< Go on if load fails.
   integer(I4P)                                      :: error           !< Error status.
   character(99)                                     :: buff            !< Character buffer.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='physical_model', val=buff, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(physical_model)')
   select case(trim(adjustl(buff)))
	case('Electromagnetic', 'Maxwell', 'EM', 'em', 'maxwell', 'electromagnetic')
		self%physical_model = EM_PHYSICAL_MODEL
	case('PIC', 'pic', 'ParticleInCell', 'particleincell')
      self%physical_model = PIC_PHYSICAL_MODEL
   case default
      call self%mpih%error_stop(msg=': unknown physical model "'//trim(buff)//'" in ['//INI_SECTION_NAME//'].(physical_model)')
	endselect

   if (div_corr_var == DIV_CORR_VAR_HYPER) then
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='chi', val=self%chi, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(chi)')
      !call file_parameters%get(section_name=INI_SECTION_NAME, option_name='eta', val=self%eta, error=error)
      !if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(eta)')

      !In case of hyperbolic divergence cleaning, set evmax to physical max eigenvalue
      self%evmax = self%chi*C0
   else
      !In case of no divergence cleaning or poisson cleaning, set evmax to light speed value
      self%chi = 0.0_R8P
      self%evmax = C0
   end if
   endsubroutine load_from_file
endmodule adam_prism_physics_object
