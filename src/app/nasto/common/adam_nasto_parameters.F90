!< ADAM, Navier-Stokes equations system parameters.
module adam_nasto_parameters
!< ADAM, Navier-Stokes equations system parameters.

use PENF, only : I4P

implicit none
private
public :: IC_UNIFORM
public :: IC_LEFTRIGHT
public :: IC_FLAME
public :: IC_VORTEX
public :: IC_VARS_NUMBER
public :: IC_VARS_NUMBER_MAX
public :: BC_EXTRAPOLATION
public :: BC_INFLOW
public :: BC_VARS_NUMBER
public :: BC_VARS_NUMBER_MAX
public :: BCS_VISCOUS
public :: BCS_EULER
public :: BCS_VARS_NUMBER
public :: BCS_VARS_NUMBER_MAX
public :: IWENO_FROM_SCHEME
public :: IRHO
public :: IU
public :: IV
public :: IW
public :: IYA
public :: ITEM
public :: IPRES
public :: IENTA
public :: ISOUN

! initial conditions
integer(I4P), parameter :: IC_UNIFORM         = 1_I4P
integer(I4P), parameter :: IC_LEFTRIGHT       = 2_I4P
integer(I4P), parameter :: IC_FLAME           = 3_I4P
integer(I4P), parameter :: IC_VORTEX          = 4_I4P
integer(I4P), parameter :: IC_VARS_NUMBER(4)  = [6, 11, 12, 3]
integer(I4P), parameter :: IC_VARS_NUMBER_MAX = 12 !maxval(IC_VARS_NUMBER) !< Maximum number of variables needed for IC.

! boundary conditions
integer(I4P), parameter :: BC_EXTRAPOLATION   = 1_I4P
integer(I4P), parameter :: BC_INFLOW          = 2_I4P
integer(I4P), parameter :: BC_VARS_NUMBER(2)  = [0, 5]
integer(I4P), parameter :: BC_VARS_NUMBER_MAX = 5 !maxval(BC_VARS_NUMBER) !< Maximum number of variables needed for BC.

! boundary conditions IB
integer(I4P), parameter :: BCS_VISCOUS         = 1_I4P
integer(I4P), parameter :: BCS_EULER           = 2_I4P
integer(I4P), parameter :: BCS_VARS_NUMBER(2)  = [3, 0]
integer(I4P), parameter :: BCS_VARS_NUMBER_MAX = 3 !maxval(BCS_VARS_NUMBER) !< Maximum number of variables needed for BCS.

! WENO
integer(I4P), parameter :: IWENO_FROM_SCHEME(6) = [1,2,3,4,3,3]

! conservative variables are arranged as follows:
! q(1): rho
! q(2): rho * u
! q(3): rho * v
! q(4): rho * w
! q(5): rho * E
! q(6): rho * Ya, specific density of specie
! auxiliary variables are arranged as follows:
! q_aux(1): rho
! q_aux(2): u
! q_aux(3): v
! q_aux(4): w
! q_aux(5): ya
! q_aux(6): tem
! q_aux(7): pressure
! q_aux(8): entalpy
! q_aux(9): sound speed

integer(I4P), parameter :: IRHO  = 1
integer(I4P), parameter :: IU    = 2
integer(I4P), parameter :: IV    = 3
integer(I4P), parameter :: IW    = 4
integer(I4P), parameter :: IYA   = 5
integer(I4P), parameter :: ITEM  = 6
integer(I4P), parameter :: IPRES = 7
integer(I4P), parameter :: IENTA = 8
integer(I4P), parameter :: ISOUN = 9
endmodule adam_nasto_parameters
