!< ADAM, FLUME parameters: variables indexes, accepted option values, shared helpers.
module adam_flume_parameters
!< ADAM, FLUME parameters: variables indexes, accepted option values, shared helpers.
!<
!< Every option value FLUME accepts from the INI file is a named constant here, so each parser matches against
!< one list and its error message can name the accepted spellings.

! third party modules
use :: penf, only : I4P, R8P

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
public :: IQ_S
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
public :: IA_S
public :: NV_EULER
public :: NV_AUX
public :: NV_MHD
public :: NV_MHD_GLM
public :: NV_MHD_EGLM
public :: NV_AUX_MHD
public :: NV_EULER_DE
public :: NV_AUX_DE
public :: NS_MAX
public :: COOLING_NODES_MAX
public :: NV_EULER_PS
public :: NV_AUX_PS
public :: NV_EULER_DE_PS
public :: NV_AUX_DE_PS
public :: S_MAX
public :: MODEL_EULER
public :: MODEL_MHD
public :: MODEL_MHD_GLM
public :: MODEL_MHD_EGLM
public :: MODEL_EULER_DE
public :: MODEL_EULER_PS
public :: MODEL_EULER_DE_PS
public :: MODEL_EULER_ISO
public :: EOS_IDEAL
public :: EOS_ISOTHERMAL
public :: PHYSICAL_MODEL_EULER
public :: PHYSICAL_MODEL_MHD_IDEAL
public :: PHYSICAL_MODEL_EULER_DE
public :: DUAL_ENERGY_ETA
public :: DIVERGENCE_CONTROL_GLM
public :: DIVERGENCE_CONTROL_EGLM
public :: DIVERGENCE_CONTROL_NONE
public :: GLM_CH_CHECK_WARNING
public :: GLM_CH_CHECK_ERROR
public :: GLM_DAMPING_LENGTH_MIN_CELL
public :: SCHEME_SPACE_WENO
public :: SCHEME_SPACE_WENO_RIEMANN
public :: RECON_CHARACTERISTIC
public :: RECON_CONSERVATIVE
public :: RECON_PRIMITIVE
public :: RIEMANN_SOLVER_LLF
public :: RIEMANN_SOLVER_HLL
public :: RIEMANN_SOLVER_HLLC
public :: RIEMANN_SOLVER_HLLD
public :: FLUX_CORRECTION_6TH
public :: FLUX_CORRECTION_4TH
public :: FLUX_CORRECTION_NONE
public :: FLUX_CORRECTION_SENSOR_WENO
public :: FLUX_CORRECTION_SENSOR_NONE
public :: FLUX_CORRECTION_SENSOR_TAU
public :: POSITIVITY_LIMITER_CELL
public :: POSITIVITY_LIMITER_KAPPA
public :: POSITIVITY_LIMITER_NONE
public :: REFLUX_POLICY_ALWAYS
public :: REFLUX_POLICY_TROUBLED
public :: REFLUX_CADENCE_STEP
public :: REFLUX_CADENCE_STAGE
public :: strip_control
public :: is_hydrodynamic
public :: is_isothermal
public :: has_entropy

! conservative variables
integer(I4P), parameter :: IQ_R   = 1_I4P !< Density.
integer(I4P), parameter :: IQ_RU  = 2_I4P !< Momentum, x component.
integer(I4P), parameter :: IQ_RV  = 3_I4P !< Momentum, y component.
integer(I4P), parameter :: IQ_RW  = 4_I4P !< Momentum, z component.
integer(I4P), parameter :: IQ_RE  = 5_I4P !< Total energy per unit volume (MHD: includes the magnetic energy |B|^2/2).
integer(I4P), parameter :: IQ_BX  = 6_I4P !< Magnetic field, x component (MHD).
integer(I4P), parameter :: IQ_BY  = 7_I4P !< Magnetic field, y component (MHD).
integer(I4P), parameter :: IQ_BZ  = 8_I4P !< Magnetic field, z component (MHD).
integer(I4P), parameter :: IQ_PSI = 9_I4P !< GLM divergence-cleaning scalar (MHD with GLM or EGLM).
integer(I4P), parameter :: IQ_S   = 6_I4P !< Advected entropy rho s, s = p / rho^gamma (Euler with dual energy, #91).
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
integer(I4P), parameter :: IA_S  = 9_I4P  !< Specific entropy s = p / rho^gamma (Euler with dual energy, #91).
! dimensions
integer(I4P), parameter :: NV_EULER   = 5_I4P  !< Conservative variables number, Euler model.
integer(I4P), parameter :: NV_AUX     = 8_I4P  !< Auxiliary variables number, Euler model.
integer(I4P), parameter :: NV_MHD     = 8_I4P  !< Conservative variables number, MHD without divergence control.
integer(I4P), parameter :: NV_MHD_GLM = 9_I4P  !< Conservative variables number, MHD with GLM divergence cleaning.
integer(I4P), parameter :: NV_MHD_EGLM = 9_I4P !< Conservative variables number, MHD with EGLM divergence cleaning.
integer(I4P), parameter :: NV_AUX_MHD = 11_I4P !< Auxiliary variables number, MHD (both variants).
integer(I4P), parameter :: NV_EULER_DE = 6_I4P !< Conservative variables number, Euler with dual energy (#91).
integer(I4P), parameter :: NV_AUX_DE   = 9_I4P !< Auxiliary variables number, Euler with dual energy (#91).
integer(I4P), parameter :: NS_MAX         = 4_I4P  !< Passive scalars capacity of the scalar models (#92, D-M10-1).
integer(I4P), parameter :: COOLING_NODES_MAX = 64_I4P !< Cooling table capacity: the device kernels take fixed-size
                                                     !< tables (constant-shape firstprivate arrays, #93).
integer(I4P), parameter :: NV_EULER_PS    = 9_I4P  !< Conservative variables number, Euler + NS_MAX scalars (#92).
integer(I4P), parameter :: NV_AUX_PS      = 12_I4P !< Auxiliary variables number, Euler + NS_MAX scalars (#92).
integer(I4P), parameter :: NV_EULER_DE_PS = 10_I4P !< Conservative variables number, euler-de + NS_MAX scalars (#92).
integer(I4P), parameter :: NV_AUX_DE_PS   = 13_I4P !< Auxiliary variables number, euler-de + NS_MAX scalars (#92).
integer(I4P), parameter :: S_MAX      = 5_I4P  !< Maximum WENO stencil half-width (weno-u-9).
! physical models: the id the host controllers dispatch on, never inside a kernel (issue #41, section 4)
integer(I4P), parameter :: MODEL_EULER   = 1_I4P !< Compressible Euler (nv = 5).
integer(I4P), parameter :: MODEL_MHD     = 2_I4P !< Ideal compressible MHD without divergence control (nv = 8).
integer(I4P), parameter :: MODEL_MHD_GLM = 3_I4P !< Ideal compressible MHD with GLM divergence cleaning (nv = 9).
integer(I4P), parameter :: MODEL_MHD_EGLM = 4_I4P !< Ideal compressible MHD with EGLM divergence cleaning (nv = 9).
integer(I4P), parameter :: MODEL_EULER_DE = 5_I4P !< Compressible Euler with dual energy, advected entropy (nv = 6, #91).
integer(I4P), parameter :: MODEL_EULER_PS = 6_I4P !< Euler with NS_MAX passive scalar slots (nv = 9, #92).
integer(I4P), parameter :: MODEL_EULER_DE_PS = 7_I4P !< euler-de with NS_MAX passive scalar slots (nv = 10, #92).
integer(I4P), parameter :: MODEL_EULER_ISO = 8_I4P !< Isothermal Euler, `[eos] type = isothermal` (nv = 5, the energy
                                                   !< slot inert: zero flux, zero source, held at 0; #100, D-M13-1).
! accepted option values
character(len=5),  parameter :: PHYSICAL_MODEL_EULER="euler"          !< [physics].(physical_model): compressible Euler.
character(len=9),  parameter :: PHYSICAL_MODEL_MHD_IDEAL="mhd-ideal"  !< [physics].(physical_model): ideal MHD.
character(len=8),  parameter :: PHYSICAL_MODEL_EULER_DE="euler-de"    !< [physics].(physical_model): compressible
                                                                      !< Euler with dual energy (issue #91).
character(len=5),  parameter :: EOS_IDEAL="ideal"                    !< [eos].(type): calorically perfect gas (default).
character(len=10), parameter :: EOS_ISOTHERMAL="isothermal"           !< [eos].(type): P = c0^2 rho (issue #100).
real(R8P),         parameter :: DUAL_ENERGY_ETA=5.e-2_R8P             !< Default [physics].(dual_energy_eta): in a cell
                                                                      !< whose internal energy is below eta times its E
                                                                      !< the entropy is a floor of the internal energy
                                                                      !< (#91 D-M9-1, D-M9-3: Enzo's 1e-3 never acts in
                                                                      !< the #88 collapse, 0.05 does; Cholla 0.035).
character(len=3),  parameter :: DIVERGENCE_CONTROL_GLM="glm"          !< [mhd].(divergence_control): GLM cleaning.
character(len=4),  parameter :: DIVERGENCE_CONTROL_EGLM="eglm"        !< [mhd].(divergence_control): EGLM cleaning.
character(len=4),  parameter :: DIVERGENCE_CONTROL_NONE="none"        !< [mhd].(divergence_control): none.
character(len=7),  parameter :: GLM_CH_CHECK_WARNING="warning"        !< [mhd].(glm_ch_check): warn if c_h is slow.
character(len=5),  parameter :: GLM_CH_CHECK_ERROR="error"            !< [mhd].(glm_ch_check): stop if c_h is slow.
character(len=8),  parameter :: GLM_DAMPING_LENGTH_MIN_CELL="min-cell" !< [mhd].(glm_damping_length): min spacing.
character(len=4),  parameter :: SCHEME_SPACE_WENO="weno"              !< [numerics].(scheme_space): WENO flux splitting.
character(len=14), parameter :: RECON_CHARACTERISTIC="characteristic" !< [numerics].(reconstruction_variables).
character(len=12), parameter :: RECON_CONSERVATIVE="conservative"     !< [numerics].(reconstruction_variables).
character(len=12), parameter :: SCHEME_SPACE_WENO_RIEMANN="weno-riemann" !< [numerics].(scheme_space): WENO
                                                                         !< interpolation + Riemann flux + correction.
character(len=9),  parameter :: RECON_PRIMITIVE="primitive"           !< [numerics].(reconstruction_variables).
character(len=3),  parameter :: RIEMANN_SOLVER_LLF="llf"              !< [numerics].(riemann_solver): local Lax-Friedrichs.
character(len=3),  parameter :: RIEMANN_SOLVER_HLL="hll"              !< [numerics].(riemann_solver): HLL (Einfeldt).
character(len=4),  parameter :: RIEMANN_SOLVER_HLLC="hllc"            !< [numerics].(riemann_solver): HLLC (Euler).
character(len=4),  parameter :: RIEMANN_SOLVER_HLLD="hlld"            !< [numerics].(riemann_solver): HLLD (MHD).
character(len=3),  parameter :: FLUX_CORRECTION_6TH="6th"             !< [numerics].(flux_correction): 6th order.
character(len=3),  parameter :: FLUX_CORRECTION_4TH="4th"             !< [numerics].(flux_correction): 4th order.
character(len=4),  parameter :: FLUX_CORRECTION_NONE="none"           !< [numerics].(flux_correction): none (2nd order).
character(len=4),  parameter :: FLUX_CORRECTION_SENSOR_WENO="weno"    !< [numerics].(flux_correction_sensor): WENO weights.
character(len=4),  parameter :: FLUX_CORRECTION_SENSOR_NONE="none"    !< [numerics].(flux_correction_sensor): always on.
real(R8P),         parameter :: FLUX_CORRECTION_SENSOR_TAU=0.2_R8P    !< Correction off at a face where the WENO weights
                                                                      !< ratio min_k w_k/d_k < TAU (issue #47, M3-P0).
character(len=4),  parameter :: POSITIVITY_LIMITER_NONE="none"        !< [numerics].(positivity_limiter): off.
character(len=6),  parameter :: REFLUX_POLICY_ALWAYS="always"          !< [numerics].(reflux_policy): every face cell.
character(len=8),  parameter :: REFLUX_POLICY_TROUBLED="troubled"      !< [numerics].(reflux_policy): the face cells the
                                                                       !< WENO sensor or the positivity limiter flags.
character(len=4),  parameter :: REFLUX_CADENCE_STEP="step"             !< [numerics].(reflux_cadence): once per step.
character(len=5),  parameter :: REFLUX_CADENCE_STAGE="stage"           !< [numerics].(reflux_cadence): every stage
                                                                       !< (issue #81, P5c).
character(len=4),  parameter :: POSITIVITY_LIMITER_CELL="cell"        !< [numerics].(positivity_limiter): cell-based
                                                                      !< parametrised flux limiter (issue #47, D-9).
real(R8P),         parameter :: POSITIVITY_LIMITER_KAPPA=0.1_R8P      !< Floor of the limiter: rho, the internal
                                                                      !< energy >= KAPPA backbone (issue #47, M3-P5c);
                                                                      !< relative only, no unit (issue #49).

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

   elemental function is_hydrodynamic(model) result(yes)
   !< Return true for the hydrodynamic models (Euler, Euler with dual energy, issue #91), false for the MHD ones: the
   !< predicate behind every "no magnetic field" branch of the host controllers.
   integer(I4P), intent(in) :: model !< Physical model id.
   logical                  :: yes   !< Hydrodynamic model.

   yes = model == MODEL_EULER .or. model == MODEL_EULER_DE .or. model == MODEL_EULER_PS .or. model == MODEL_EULER_DE_PS &
         .or. model == MODEL_EULER_ISO
   endfunction is_hydrodynamic

   elemental function is_isothermal(model) result(yes)
   !< Return true for the isothermal models (issue #100), whose energy slot is inert: no thermodynamics beyond P(rho).
   integer(I4P), intent(in) :: model !< Physical model id.
   logical                  :: yes   !< Isothermal model.

   yes = model == MODEL_EULER_ISO
   endfunction is_isothermal

   elemental function has_entropy(model) result(yes)
   !< Return true for the dual-energy models, whose slot IQ_S is the advected entropy (issue #91; with scalars, #92).
   integer(I4P), intent(in) :: model !< Physical model id.
   logical                  :: yes   !< Dual-energy model.

   yes = model == MODEL_EULER_DE .or. model == MODEL_EULER_DE_PS
   endfunction has_entropy
endmodule adam_flume_parameters
