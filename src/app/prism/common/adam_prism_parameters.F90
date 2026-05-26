!< PRISM, general parameters.
module adam_prism_parameters
!< PRISM, general parameters.

! third party modules
use :: penf

implicit none
private
public :: NV_MAX
public :: MU0
public :: EPS0
public :: MU0_SQ
public :: MU0_SQ_I2
public :: EPS0_SQ
public :: EPS0_SQ_I2
public :: C0
public :: E_CHARGE
public :: E_MASS
public :: K_B
public :: PI
public :: EV
public :: ER
public :: EL
public :: IEV
public :: IEV_D
public :: IEV_B
public :: IEV_D_B
public :: IERL
public :: IERL_D
public :: IERL_B
public :: IERL_D_B

integer(I4P), parameter :: NV_MAX     = 11_I4P                        !< Maximum number of variables for static arrays dimensioning.
real(R8P),    parameter :: MU0        = 1.256637e-6_R8P               !< Vacuum magnetic permeability.
real(R8P),    parameter :: EPS0       = 8.854187e-12_R8P              !< Vacuum dielectric constant.
real(R8P),    parameter :: C0         = 1._R8P/sqrt(MU0*EPS0)         !< Vacuum speed of light.
real(R8P),    parameter :: E_CHARGE   = -1.6e-19_R8P                  !< Elettron charge value.
real(R8P),    parameter :: E_MASS     = 9.11e-31_R8P                  !< Elettron mass value.
real(R8P),    parameter :: Q_OVER_M   = E_CHARGE/E_MASS               !< Q/m value for elettron.
real(R8P),    parameter :: K_B        = 1.380649e-23_R8P              !< Boltzmann coefficient.
real(R8P),    parameter :: T_E        = 300000_R8P                    !< Plasma Temperature.
real(R8P),    parameter :: VTH        = sqrt(k_B * T_e / e_mass)      !< Thermal Velocity value.
real(R8P),    parameter :: PI         = 3.141592653589793_R8P         !< Pi value.
real(R8P),    parameter :: MU0_SQ     = MU0**0.5_R8P                  !< Mu0 square root
real(R8P),    parameter :: MU0_SQ_I2  = 1._R8P/(2._R8P*MU0**0.5_R8P)  !< Half of MU0 square root inverse
real(R8P),    parameter :: EPS0_SQ    = EPS0**0.5_R8P                 !< EPS0 square root
real(R8P),    parameter :: EPS0_SQ_I2 = 1._R8P/(2._R8P*EPS0**0.5_R8P) !< Half of EPS0 square root inverse

real(R8P), target       :: EV(6)    = [0._R8P,0._R8P,C0,C0,-C0,-C0]  !< Eigenvalues.
real(R8P), target       :: ER(6,6,3)=reshape([1._R8P,0._R8P,0._R8P         ,0._R8P        ,0._R8P        ,0._R8P         ,&! x
                                              0._R8P,0._R8P,0._R8P         ,sqrt(EPS0/MU0),0._R8P        ,-sqrt(EPS0/MU0),&
                                              0._R8P,0._R8P,-sqrt(EPS0/MU0),0._R8P        ,sqrt(EPS0/MU0),0._R8P         ,&
                                              0._R8P,1._R8P,0._R8P         ,0._R8P        ,0._R8P        ,0._R8P         ,&
                                              0._R8P,0._R8P,1._R8P         ,0._R8P        ,1._R8P        ,0._R8P         ,&
                                              0._R8P,0._R8P,0._R8P         ,1._R8P        ,0._R8P        ,1._R8P         ,&
                                              !
                                              0._R8P,0._R8P,0._R8P        ,-sqrt(EPS0/MU0),0._R8P         ,sqrt(EPS0/MU0),&! y
                                              1._R8P,0._R8P,0._R8P        ,0._R8P         ,0._R8P         ,0._R8P        ,&
                                              0._R8P,0._R8P,sqrt(EPS0/MU0),0._R8P         ,-sqrt(EPS0/MU0),0._R8P        ,&
                                              0._R8P,0._R8P,1._R8P        ,0._R8P         ,1._R8P         ,0._R8P        ,&
                                              0._R8P,1._R8P,0._R8P        ,0._R8P         ,0._R8P         ,0._R8P        ,&
                                              0._R8P,0._R8P,0._R8P        ,1._R8P         ,0._R8P         ,1._R8P        ,&
                                              !
                                              0._R8P,0._R8P,0._R8P         ,sqrt(EPS0/MU0),0._R8P        ,-sqrt(EPS0/MU0),&! z
                                              0._R8P,0._R8P,-sqrt(EPS0/MU0),0._R8P        ,sqrt(EPS0/MU0),0._R8P         ,&
                                              1._R8P,0._R8P,0._R8P         ,0._R8P        ,0._R8P        ,0._R8P         ,&
                                              0._R8P,0._R8P,1._R8P         ,0._R8P        ,1._R8P        ,0._R8P         ,&
                                              0._R8P,0._R8P,0._R8P         ,1._R8P        ,0._R8P        ,1._R8P         ,&
                                              0._R8P,1._R8P,0._R8P         ,0._R8P        ,0._R8P        ,0._R8P],        &
                                             [6,6,3]) !< Right eigenvectors.
real(R8P), target       :: EL(6,6,3)=reshape([1._R8P,0._R8P                 , 0._R8P                ,0._R8P,0._R8P ,0._R8P ,&! x
                                              0._R8P,0._R8P                 , 0._R8P                ,1._R8P,0._R8P ,0._R8P ,&
                                              0._R8P,0._R8P                 ,-0.5_R8P*sqrt(MU0/EPS0),0._R8P,0.5_R8P,0._R8P ,&
                                              0._R8P,0.5_R8P*sqrt(MU0/EPS0) , 0._R8P                ,0._R8P,0._R8P ,0.5_R8P,&
                                              0._R8P,0._R8P                 , 0.5_R8P*sqrt(MU0/EPS0),0._R8P,0.5_R8P,0._R8P ,&
                                              0._R8P,-0.5_R8P*sqrt(MU0/EPS0),0._R8P                 ,0._R8P,0._R8P ,0.5_R8P,&
                                              !
                                              0._R8P                 ,1._R8P,0._R8P                 ,0._R8P ,0._R8P,0._R8P ,&! y
                                              0._R8P                 ,0._R8P,0._R8P                 ,0._R8P ,1._R8P,0._R8P ,&
                                              0._R8P                 ,0._R8P,0.5_R8P*sqrt(MU0/EPS0) ,0.5_R8P,0._R8P,0._R8P ,&
                                              -0.5_R8P*sqrt(MU0/EPS0),0._R8P,0._R8P                 ,0._R8P ,0._R8P,0.5_R8P,&
                                              0._R8P                 ,0._R8P,-0.5_R8P*sqrt(MU0/EPS0),0.5_R8P,0._R8P,0._R8P ,&
                                              0.5_R8P*sqrt(MU0/EPS0) ,0._R8P,0._R8P                 ,0._R8P ,0._R8P,0.5_R8P,&
                                              !
                                              0._R8P                 ,0._R8P                 ,1._R8P,0._R8P ,0._R8P ,0._R8P, &! z
                                              0._R8P                 ,0._R8P                 ,0._R8P,0._R8P ,0._R8P ,1._R8P, &
                                              0._R8P                 ,-0.5_R8P*sqrt(MU0/EPS0),0._R8P,0.5_R8P,0._R8P ,0._R8P, &
                                              0.5_R8P*sqrt(MU0/EPS0) ,0._R8P                 ,0._R8P,0._R8P ,0.5_R8P,0._R8P, &
                                              0._R8P                 ,0.5_R8P*sqrt(MU0/EPS0) ,0._R8P,0.5_R8P,0._R8P ,0._R8P, &
                                              -0.5_R8P*sqrt(MU0/EPS0),0._R8P                 ,0._R8P,0._R8P ,0.5_R8P,0._R8P],&
                                             [6,6,3]) !< Right eigenvectors.
real(R8P), target       :: IEV(6)         = [1._R8P,1._R8P,1._R8P,1._R8P,1._R8P,1._R8P]                !< Identity eigenvalues.
real(R8P), target       :: IEV_D(7)       = [1._R8P,1._R8P,1._R8P,1._R8P,1._R8P,1._R8P,1._R8P]
                                                                                                       !< Identity eigenvalues with
                                                                                                       !< D divergence cleaning.
real(R8P), target       :: IEV_B(7)       = [1._R8P,1._R8P,1._R8P,1._R8P,1._R8P,1._R8P,1._R8P]
                                                                                                       !< Identity eigenvalues with
                                                                                                       !< B divergence cleaning.
real(R8P), target       :: IEV_D_B(8)     = [1._R8P,1._R8P,1._R8P,1._R8P,1._R8P,1._R8P,1._R8P,1._R8P]
                                                                                                       !< Identity eigenvalues with
                                                                                                       !< D & B divergence cleaning.
real(R8P), target       :: IERL(6,6,3)    =reshape([1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &       ! x
                                                    0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P, &
                                                    1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &       ! y
                                                    0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P, &
                                                    1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &       ! z
                                                    0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P],&
                                                 [6,6,3])   !< Identity eigenvectors.
real(R8P), target      :: IERL_D(7,7,3)   =reshape([1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, & ! x
                                                    0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P, &
                                                    1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, & ! y
                                                    0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P, &
                                                    1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, & ! z
                                                    0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P],&
                                                 [7,7,3])   !< Identity eigenvectors, D divergence cleaning.
real(R8P), target      :: IERL_B(7,7,3)   =reshape([1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, & ! x
                                                    0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P, &
                                                    1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, & ! y
                                                    0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P, &
                                                    1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, & ! z
                                                    0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P],&
                                                 [7,7,3])   !< Identity eigenvectors, B divergence cleaning.
real(R8P), target      :: IERL_D_B(8,8,3) =reshape([1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, & ! x
                                                    0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P, &
                                                    1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, & ! y
                                                    0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P, &
                                                    1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, & ! z
                                                    0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P,0._R8P, &
                                                    0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,0._R8P,1._R8P],&
                                                 [8,8,3])   !< Identity eigenvectors, D & B divergence cleaning.
endmodule adam_prism_parameters
