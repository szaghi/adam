!< PRISM, general parameters.
module adam_prism_parameters
!< PRISM, general parameters.

use PENF

implicit none
save
private
public :: NV_MAX, MU0, EPS0, C0, PI, EV, ER, EL

integer(I4P), parameter :: NV_MAX   = 11_I4P                        !< Maximum number of variables for static arrays dimensioning.
real(R8P),    parameter :: MU0      = 1.256637e-6_R8P               !< Vacuum magnetic permeability.
real(R8P),    parameter :: EPS0     = 8.854187e-12_R8P              !< Vacuum dielectric constant.
real(R8P),    parameter :: C0       = 1._R8P/sqrt(MU0*EPS0)         !< Vacuum speed of light.
real(R8P),    parameter :: E_CHARGE = -1.6e-19_R8P                  !< Elettron charge value.
real(R8P),    parameter :: E_MASS   = 9.11e-31_R8P                  !< Elettron mass value.
real(R8P),    parameter :: Q_OVER_M = E_CHARGE/E_MASS               !< Q/m value for elettron.
real(R8P),    parameter :: K_B      = 1.380649e-23_R8P              !< Boltzmann coefficient.
real(R8P),    parameter :: T_E      = 300000_R8P                    !< Plasma Temperature.
real(R8P),    parameter :: VTH      = sqrt(k_B * T_e / e_mass)      !< Thermal Velocity value.
real(R8P),    parameter :: PI       = 3.141592653589793_R8P         !< Pi value.
real(R8P),    parameter :: EV(6)    = [0._R8P,0._R8P,C0,C0,-C0,-C0] !< Eigenvalues.
real(R8P),    parameter :: ER(6,6,3)=reshape([1._R8P,0._R8P,0._R8P         ,0._R8P        ,0._R8P        ,0._R8P         ,&! x
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
                                             [6,6,3])          !< Right eigenvectors.
real(R8P),    parameter :: EL(6,6,3)=reshape([1._R8P,0._R8P                 , 0._R8P                ,0._R8P,0._R8P ,0._R8P ,&! x
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
                                              0._R8P                ,0._R8P                 ,1._R8P,0._R8P ,0._R8P ,0._R8P, &! z
                                              0._R8P                ,0._R8P                 ,0._R8P,0._R8P ,0._R8P ,1._R8P, &
                                              0._R8P                ,-0.5_R8P*sqrt(MU0/EPS0),0._R8P,0.5_R8P,0._R8P ,0._R8P, &
                                              0.5_R8P*sqrt(MU0/EPS0),0._R8P                 ,0._R8P,0._R8P ,0.5_R8P,0._R8P, &
                                              0._R8P                ,0.5_R8P*sqrt(MU0/EPS0) ,0._R8P,0.5_R8P,0._R8P ,0._R8P, &
                                              0.5_R8P*sqrt(MU0/EPS0),0._R8P                 ,0._R8P,0._R8P ,0.5_R8P,0._R8P],&
                                             [6,6,3])          !< Right eigenvectors.
endmodule adam_prism_parameters
