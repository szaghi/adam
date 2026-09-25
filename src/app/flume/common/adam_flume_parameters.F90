!< ADAM, FLUME parameters: variables indexes, accepted option values, shared helpers.
module adam_flume_parameters
!< ADAM, FLUME parameters: variables indexes, accepted option values, shared helpers.
!<
!< Every option value FLUME accepts from the INI file is a named constant here, so each parser matches against
!< one list and its error message can name the accepted spellings.

! third party modules
use :: penf, only : I4P

implicit none
private
public :: IQ_R
public :: IQ_RU
public :: IQ_RV
public :: IQ_RW
public :: IQ_RE
public :: IQ_BX
public :: IQ_BY
public :: IQ_BZ
public :: IQ_PSI
public :: IA_R
public :: IA_U
public :: IA_V
public :: IA_W
public :: IA_P
public :: IA_T
public :: IA_H
public :: IA_A
public :: IA_BX
public :: IA_BY
public :: IA_BZ
public :: NV_EULER
public :: NV_AUX
public :: NV_MHD
public :: NV_MHD_GLM
public :: NV_AUX_MHD
public :: S_MAX
public :: MODEL_EULER
public :: MODEL_MHD
public :: MODEL_MHD_GLM
public :: PHYSICAL_MODEL_EULER
public :: PHYSICAL_MODEL_MHD_IDEAL
public :: DIVERGENCE_CONTROL_GLM
public :: DIVERGENCE_CONTROL_NONE
public :: GLM_CH_CHECK_WARNING
public :: GLM_CH_CHECK_ERROR
public :: SCHEME_SPACE_WENO
public :: RECON_CHARACTERISTIC
public :: RECON_CONSERVATIVE
public :: strip_control

! conservative variables
integer(I4P), parameter :: IQ_R   = 1_I4P !< Density.
integer(I4P), parameter :: IQ_RU  = 2_I4P !< Momentum, x component.
integer(I4P), parameter :: IQ_RV  = 3_I4P !< Momentum, y component.
integer(I4P), parameter :: IQ_RW  = 4_I4P !< Momentum, z component.
integer(I4P), parameter :: IQ_RE  = 5_I4P !< Total energy per unit volume (MHD: includes the magnetic energy |B|^2/2).
integer(I4P), parameter :: IQ_BX  = 6_I4P !< Magnetic field, x component (MHD).
integer(I4P), parameter :: IQ_BY  = 7_I4P !< Magnetic field, y component (MHD).
integer(I4P), parameter :: IQ_BZ  = 8_I4P !< Magnetic field, z component (MHD).
integer(I4P), parameter :: IQ_PSI = 9_I4P !< GLM divergence-cleaning scalar (MHD with GLM).
! auxiliary (primitive and derived) variables
integer(I4P), parameter :: IA_R  = 1_I4P  !< Density.
integer(I4P), parameter :: IA_U  = 2_I4P  !< Velocity, x component.
integer(I4P), parameter :: IA_V  = 3_I4P  !< Velocity, y component.
integer(I4P), parameter :: IA_W  = 4_I4P  !< Velocity, z component.
integer(I4P), parameter :: IA_P  = 5_I4P  !< Pressure (thermal).
integer(I4P), parameter :: IA_T  = 6_I4P  !< Temperature.
integer(I4P), parameter :: IA_H  = 7_I4P  !< Total specific enthalpy (MHD: (E + p + |B|^2/2) / rho).
integer(I4P), parameter :: IA_A  = 8_I4P  !< Speed of sound.
integer(I4P), parameter :: IA_BX = 9_I4P  !< Magnetic field, x component (MHD).
integer(I4P), parameter :: IA_BY = 10_I4P !< Magnetic field, y component (MHD).
integer(I4P), parameter :: IA_BZ = 11_I4P !< Magnetic field, z component (MHD).
! dimensions
integer(I4P), parameter :: NV_EULER   = 5_I4P  !< Conservative variables number, Euler model.
integer(I4P), parameter :: NV_AUX     = 8_I4P  !< Auxiliary variables number, Euler model.
integer(I4P), parameter :: NV_MHD     = 8_I4P  !< Conservative variables number, MHD without divergence control.
integer(I4P), parameter :: NV_MHD_GLM = 9_I4P  !< Conservative variables number, MHD with GLM divergence cleaning.
integer(I4P), parameter :: NV_AUX_MHD = 11_I4P !< Auxiliary variables number, MHD (both variants).
integer(I4P), parameter :: S_MAX      = 5_I4P  !< Maximum WENO stencil half-width (weno-u-9).
! physical models: the id the host controllers dispatch on, never inside a kernel (issue #41, section 4)
integer(I4P), parameter :: MODEL_EULER   = 1_I4P !< Compressible Euler (nv = 5).
integer(I4P), parameter :: MODEL_MHD     = 2_I4P !< Ideal compressible MHD without divergence control (nv = 8).
integer(I4P), parameter :: MODEL_MHD_GLM = 3_I4P !< Ideal compressible MHD with GLM divergence cleaning (nv = 9).
! accepted option values
character(len=5),  parameter :: PHYSICAL_MODEL_EULER="euler"          !< [physics].(physical_model): compressible Euler.
character(len=9),  parameter :: PHYSICAL_MODEL_MHD_IDEAL="mhd-ideal"  !< [physics].(physical_model): ideal MHD.
character(len=3),  parameter :: DIVERGENCE_CONTROL_GLM="glm"          !< [mhd].(divergence_control): GLM cleaning.
character(len=4),  parameter :: DIVERGENCE_CONTROL_NONE="none"        !< [mhd].(divergence_control): none.
character(len=7),  parameter :: GLM_CH_CHECK_WARNING="warning"        !< [mhd].(glm_ch_check): warn if c_h is slow.
character(len=5),  parameter :: GLM_CH_CHECK_ERROR="error"            !< [mhd].(glm_ch_check): stop if c_h is slow.
character(len=4),  parameter :: SCHEME_SPACE_WENO="weno"              !< [numerics].(scheme_space): WENO flux splitting.
character(len=14), parameter :: RECON_CHARACTERISTIC="characteristic" !< [numerics].(reconstruction_variables).
character(len=12), parameter :: RECON_CONSERVATIVE="conservative"     !< [numerics].(reconstruction_variables).

contains
   ! public procedures
   pure function strip_control(string) result(stripped)
   !< Return a copy of an INI value with control characters (e.g. a CRLF carriage return) blanked.
   !<
   !< `trim(adjustl())` does not remove a trailing carriage return, so a value read from a CRLF file would miss every
   !< `case` of an option parser and hit its fatal `case default`. Blank them before matching.
   character(len=*), intent(in) :: string   !< Raw value as read from the INI file.
   character(len=len(string))   :: stripped !< Value with control characters blanked.
   integer(I4P)                 :: c        !< Counter.

   stripped = string
   do c=1, len(stripped)
      if (iachar(stripped(c:c)) < 32) stripped(c:c) = ' '
   enddo
   endfunction strip_control
endmodule adam_flume_parameters
