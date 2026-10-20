!< ADAM, FLUME physics class definition.
module adam_flume_physics_object
!< ADAM, FLUME physics class definition.
!<
!< The physical model is the ONE predicate that decides the state vector width `nv` (and `nv_aux`): the name arrays
!< of the common object are built from the same predicate, so allocation and names cannot disagree. For
!< `mhd-ideal` the predicate is `(physical_model, [mhd].(divergence_control))` (issue #41, D-4): `none` gives
!< `MODEL_MHD` (`nv = 8`), `glm` gives `MODEL_MHD_GLM` (`nv = 9`), `eglm` gives `MODEL_MHD_EGLM` (`nv = 9`).
!< `euler-de` (issue #91) is the Euler model with dual energy: `nv = 6`, the sixth variable the advected entropy
!< `rho s`, `s = p / rho^gamma`; `[physics].(dual_energy_eta)` (optional, default `DUAL_ENERGY_ETA`) is its switch.
!< `[physics].(passive_scalars) = N` (issue #92, optional, 0 by default, 0 <= N <= NS_MAX, the hydrodynamic models only)
!< selects the scalar variant of `euler` or `euler-de`: `NS_MAX` more slots at the end of the state (`nv = 9` or `10`),
!< the first `N` carrying `rho s_k`, the others zero.

! ADAM singleton objects
use :: adam_mpih_global,         only : mpih
! FLUME modules
use :: adam_flume_euler_library, only : primitive_to_conservative
use :: adam_flume_euler_iso_library, only : iso_primitive_to_conservative=>primitive_to_conservative
use :: adam_flume_mhd_library,   only : mhd_primitive_to_conservative
use :: adam_flume_cooling_object, only : flume_cooling_object
use :: adam_flume_dissipation_object, only : flume_dissipation_object
use :: adam_flume_gravity_object, only : flume_gravity_object
use :: adam_flume_mhd_object,    only : flume_mhd_object
use :: adam_flume_parameters,    only : DUAL_ENERGY_ETA, EOS_IDEAL, EOS_ISOTHERMAL, MODEL_EULER_ISO,                       &
                                        has_entropy, IQ_PSI, IQ_S, MODEL_EULER, MODEL_EULER_DE,                         &
                                        MODEL_EULER_DE_PS, MODEL_EULER_PS, MODEL_MHD, MODEL_MHD_EGLM, MODEL_MHD_GLM,     &
                                        NS_MAX, NV_AUX, NV_AUX_DE, NV_AUX_DE_PS, NV_AUX_MHD, NV_AUX_PS, NV_EULER,        &
                                        NV_EULER_DE, NV_EULER_DE_PS, NV_EULER_PS, NV_MHD, NV_MHD_EGLM, NV_MHD_GLM,       &
                                        PHYSICAL_MODEL_EULER, PHYSICAL_MODEL_EULER_DE, PHYSICAL_MODEL_MHD_IDEAL,          &
                                        strip_control
! third party modules
use :: finer,                    only : file_ini
use :: penf,                     only : I4P, R8P, str

implicit none
private
public :: flume_physics_object
public :: primitive_state_to_conservative

character(len=7), parameter :: INI_SECTION_NAME="physics" !< INI (config) file section name containing physics configs.

type :: flume_physics_object
   !< FLUME physics class definition.
   character(:), allocatable :: physical_model !< Physical model.
   integer(I4P)              :: model=0_I4P    !< Physical model id (MODEL_*), the host controllers dispatch on it.
   real(R8P)                 :: cp=0._R8P      !< Specific heat at constant pressure.
   real(R8P)                 :: cv=0._R8P      !< Specific heat at constant volume.
   real(R8P)                 :: gamma=0._R8P   !< Specific heats ratio, cp/cv.
   real(R8P)                 :: R=0._R8P       !< Gas constant, cp-cv.
   character(:), allocatable :: eos            !< `[eos].(type)`: ideal (default) or isothermal (issue #100).
   real(R8P)                 :: sound_speed=0._R8P !< `[eos].(sound_speed)`: the isothermal c0.
   real(R8P)                 :: cs2=0._R8P     !< c0^2, the EOS parameter the isothermal kernels take as `gamma`.
   integer(I4P)              :: nv=0_I4P       !< Conservative variables number.
   integer(I4P)              :: nv_aux=0_I4P   !< Auxiliary variables number.
   integer(I4P)              :: nv_mhd=0_I4P   !< Last slot of the MHD state: NV_EULER (hydrodynamic models) or nv (MHD).
   real(R8P)                 :: dual_energy_eta=DUAL_ENERGY_ETA !< Dual-energy switch (euler-de, issue #91).
   integer(I4P)              :: passive_scalars=0_I4P !< Passive scalars in use, 0..NS_MAX (issue #92).
   integer(I4P)              :: iq_scalar1=0_I4P      !< Slot of the first passive scalar (0 without the scalar slots).
   integer(I4P)              :: nv_entropy=0_I4P      !< Last entropy slot: IQ_S with dual energy, else nv_mhd.
   type(flume_mhd_object)    :: mhd            !< MHD configs (`[mhd]`, loaded with mhd-ideal only).
   type(flume_dissipation_object) :: dissipation !< Dissipative coefficients (`[physics]`, issue #65).
   type(flume_gravity_object) :: gravity !< External gravity (`[gravity]`, issue #73).
   type(flume_cooling_object) :: cooling !< Radiative cooling (`[cooling]`, issue #93).
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize physics.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype flume_physics_object

contains
   ! public methods
   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(flume_physics_object), intent(in) :: self             !< Physics.
   character(len=:), allocatable           :: desc             !< Description.
   character(len=1), parameter             :: NL=new_line('a') !< New line character.

   desc =       mpih%myrankstr//'Physics main data'//NL
   desc = desc//mpih%myrankstr//'  physical_model: '//self%physical_model//NL
   desc = desc//mpih%myrankstr//'  eos:            '//self%eos
   if (self%eos == EOS_ISOTHERMAL) desc = desc//', sound_speed '//trim(str(self%sound_speed))
   desc = desc//NL
   desc = desc//mpih%myrankstr//'  cp, cv:         '//trim(str(self%cp))//', '//trim(str(self%cv))//NL
   desc = desc//mpih%myrankstr//'  gamma, R:       '//trim(str(self%gamma))//', '//trim(str(self%R))//NL
   desc = desc//mpih%myrankstr//'  nv, nv_aux:     '//trim(str(self%nv))//', '//trim(str(self%nv_aux))//NL
   if (has_entropy(self%model)) &
      desc = desc//mpih%myrankstr//'  dual energy:    eta '//trim(str(self%dual_energy_eta))//NL
   if (self%passive_scalars > 0_I4P) &
      desc = desc//mpih%myrankstr//'  passive scalars: '//trim(str(self%passive_scalars, .true.))//' of '// &
                   trim(str(NS_MAX, .true.))//' slots, from slot '//trim(str(self%iq_scalar1, .true.))//NL
   desc = desc//self%dissipation%description()//NL
   desc = desc//self%gravity%description()//NL
   desc = desc//self%cooling%description()
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize physics.
   class(flume_physics_object), intent(inout) :: self            !< Physics.
   type(file_ini),              intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   print '(A)', mpih%myrankstr//'flume_physics_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   ! issue #100, P1: the isothermal gas is an Euler model only; refused before any model loader reads its section
   if (self%eos == EOS_ISOTHERMAL .and. (self%physical_model /= PHYSICAL_MODEL_EULER .or. self%passive_scalars > 0_I4P)) &
      call mpih%error_stop(msg=': [eos].(type)='//EOS_ISOTHERMAL//' needs [physics].(physical_model)='// &
                               PHYSICAL_MODEL_EULER//' without passive_scalars (issue #100, P1)')
   select case(self%physical_model)
   case(PHYSICAL_MODEL_EULER)
      self%model  = MODEL_EULER
      self%nv     = NV_EULER
      self%nv_aux = NV_AUX
      if (self%passive_scalars > 0_I4P) then
         self%model  = MODEL_EULER_PS
         self%nv     = NV_EULER_PS
         self%nv_aux = NV_AUX_PS
      endif
   case(PHYSICAL_MODEL_EULER_DE)
      self%model  = MODEL_EULER_DE
      self%nv     = NV_EULER_DE
      self%nv_aux = NV_AUX_DE
      if (self%passive_scalars > 0_I4P) then
         self%model  = MODEL_EULER_DE_PS
         self%nv     = NV_EULER_DE_PS
         self%nv_aux = NV_AUX_DE_PS
      endif
   case(PHYSICAL_MODEL_MHD_IDEAL)
      call self%mhd%initialize(file_parameters=file_parameters)
      if (self%mhd%has_eglm) then
         self%model = MODEL_MHD_EGLM
         self%nv    = NV_MHD_EGLM
      elseif (self%mhd%has_glm) then
         self%model = MODEL_MHD_GLM
         self%nv    = NV_MHD_GLM
      else
         self%model = MODEL_MHD
         self%nv    = NV_MHD
      endif
      self%nv_aux = NV_AUX_MHD
   endselect
   if (self%eos == EOS_ISOTHERMAL) then ! issue #100, P1: the isothermal Euler model, energy slot inert (D-M13-1)
      if (file_parameters%has_section(section_name='cooling')) &
         call mpih%error_stop(msg=': [cooling] is refused with [eos].(type)='//EOS_ISOTHERMAL// &
                                  ': an isothermal gas has no thermal energy to radiate (issue #100)')
      ! the sinks' virial test and accretion use the internal energy (P4); refused by the section, before any loader
      if (file_parameters%has_section(section_name='sinks')) &
         call mpih%error_stop(msg=': [sinks] is refused with [eos].(type)='//EOS_ISOTHERMAL//' until P4 (issue #100)')
      self%model = MODEL_EULER_ISO
   endif
   self%nv_mhd = merge(self%nv, NV_EULER, self%physical_model == PHYSICAL_MODEL_MHD_IDEAL)
   self%iq_scalar1 = 0_I4P ; if (self%passive_scalars > 0_I4P) self%iq_scalar1 = self%nv - NS_MAX + 1_I4P
   self%nv_entropy = merge(IQ_S, self%nv_mhd, has_entropy(self%model))
   call self%dissipation%load_from_file(file_parameters=file_parameters, model=self%model, cp=self%cp)
   call self%gravity%load_from_file(file_parameters=file_parameters)
   call self%cooling%load_from_file(file_parameters=file_parameters, model=self%model, gamma=self%gamma, R=self%R)
   if (self%model == MODEL_EULER_ISO) then ! issue #100, P1: the energy sources of these terms are not defined yet
      if (self%dissipation%is_active) &
         call mpih%error_stop(msg=': dissipative terms ([physics] viscosity, conductivity, ...) are refused with '// &
                                  '[eos].(type)='//EOS_ISOTHERMAL//' (issue #100, P1)')
      if (self%gravity%is_active) &
         call mpih%error_stop(msg=': [gravity] is refused with [eos].(type)='//EOS_ISOTHERMAL//' until P4 (issue #100)')
   endif
   print '(A)', self%description()
   print '(A)', mpih%myrankstr//'flume_physics_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters)
   !< Load config from file; every key is required and an unknown or inconsistent value is fatal.
   class(flume_physics_object), intent(inout) :: self            !< Physics.
   type(file_ini),              intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   character(999)                             :: buff            !< Option value buffer.
   integer(I4P)                               :: error           !< Error status.

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='physical_model', val=buff, error=error)
   if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(physical_model)')
   self%physical_model = trim(adjustl(strip_control(buff)))
   select case(self%physical_model)
   case(PHYSICAL_MODEL_EULER, PHYSICAL_MODEL_EULER_DE, PHYSICAL_MODEL_MHD_IDEAL)
   case default
      call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(physical_model) "'//self%physical_model// &
                               '"; expected one of '//PHYSICAL_MODEL_EULER//', '//PHYSICAL_MODEL_EULER_DE//', '// &
                               PHYSICAL_MODEL_MHD_IDEAL)
   endselect
   ! passive scalars (issue #92): optional, 0 by default, the hydrodynamic models only
   self%passive_scalars = 0_I4P
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='passive_scalars', val=self%passive_scalars, &
                            error=error)
   if (error > 0) self%passive_scalars = 0_I4P
   if (self%passive_scalars < 0_I4P .or. self%passive_scalars > NS_MAX) &
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(passive_scalars) must lie in 0..'//trim(str(NS_MAX, .true.))// &
                               ', got '//trim(str(self%passive_scalars, .true.)))
   if (self%passive_scalars > 0_I4P .and. self%physical_model == PHYSICAL_MODEL_MHD_IDEAL) &
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(passive_scalars) needs a hydrodynamic physical_model (euler, '// &
                               'euler-de)')
   ! the dual-energy switch (issue #91): optional, read only for euler-de
   self%dual_energy_eta = DUAL_ENERGY_ETA
   if (self%physical_model == PHYSICAL_MODEL_EULER_DE) then
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='dual_energy_eta', val=self%dual_energy_eta, &
                               error=error)
      if (error > 0) self%dual_energy_eta = DUAL_ENERGY_ETA
      if (.not.(self%dual_energy_eta > 0._R8P .and. self%dual_energy_eta < 1._R8P)) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(dual_energy_eta) must lie in (0, 1), got '// &
                                  trim(str(self%dual_energy_eta)))
   endif
   ! the equation of state (issue #100): [eos] optional, type ideal (the gas below) by default, or isothermal
   self%eos = EOS_IDEAL
   if (file_parameters%has_section(section_name='eos')) then
      call file_parameters%get(section_name='eos', option_name='type', val=buff, error=error)
      if (error > 0) call mpih%error_stop(msg=': failed to load [eos].(type)')
      self%eos = trim(adjustl(strip_control(buff)))
      select case(self%eos)
      case(EOS_IDEAL)
      case(EOS_ISOTHERMAL)
         call file_parameters%get(section_name='eos', option_name='sound_speed', val=self%sound_speed, error=error)
         if (error > 0) call mpih%error_stop(msg=': failed to load [eos].(sound_speed) (needed by type = '// &
                                                 EOS_ISOTHERMAL//')')
         if (.not.(self%sound_speed > 0._R8P)) &
            call mpih%error_stop(msg=': [eos].(sound_speed) must be positive, got '//trim(str(self%sound_speed)))
         ! the gas keys have no meaning without an energy equation: refuse them rather than ignore them
         if (has(option_name='gamma') .or. has(option_name='cp') .or. has(option_name='cv')) &
            call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] gamma, cp, cv are not used with [eos].(type) = '// &
                                     EOS_ISOTHERMAL//': remove them')
         self%cs2   = self%sound_speed**2
         self%gamma = 1._R8P ! unused by the isothermal kernels (they take cs2); the gamma-law paths are refused
         self%R     = 1._R8P ! the temperature is then T = c0^2, in code units
         self%cp    = 0._R8P
         self%cv    = 0._R8P
         return
      case default
         call mpih%error_stop(msg=': unknown [eos].(type) "'//self%eos//'"; expected one of '//EOS_IDEAL//', '// &
                                  EOS_ISOTHERMAL)
      endselect
   endif
   ! the gas: cp and cv, or gamma alone (R = 1, the code-unit gas constant of issue #49, ND-2); giving both is fatal
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='gamma', val=self%gamma, error=error)
   if (error <= 0) then
      if (has(option_name='cp') .or. has(option_name='cv')) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] takes either gamma or cp, cv, not both')
      if (.not.(self%gamma > 1._R8P)) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] requires gamma > 1, got gamma='//trim(str(self%gamma)))
      self%R  = 1._R8P
      self%cv = 1._R8P / (self%gamma - 1._R8P)
      self%cp = self%gamma * self%cv
      return
   endif
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='cp', val=self%cp, error=error)
   if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(cp) (or give gamma)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='cv', val=self%cv, error=error)
   if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(cv) (or give gamma)')
   if (.not.(self%cp > self%cv .and. self%cv > 0._R8P)) &
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] requires cp > cv > 0, got cp='//trim(str(self%cp))// &
                               ', cv='//trim(str(self%cv)))
   self%gamma = self%cp / self%cv
   self%R     = self%cp - self%cv
   contains
      function has(option_name) result(is_present)
      !< Return true if the option is in the section.
      character(*), intent(in) :: option_name !< Option name.
      logical                  :: is_present  !< Presence.
      character(999)           :: buff        !< Option value buffer.
      integer(I4P)             :: error       !< Error status.

      call file_parameters%get(section_name=INI_SECTION_NAME, option_name=option_name, val=buff, error=error)
      is_present = error <= 0
      endfunction has
   endsubroutine load_from_file

   ! public procedures
   pure subroutine primitive_state_to_conservative(model, gamma, prim, q)
   !< Convert one primitive state (`r, u, v, w, p` [, `bx, by, bz`]) to the conservative variables of the model; with
   !< dual energy the advected entropy is `rho s = p rho^(1-gamma)` (issue #91).
   !<
   !< Initialisation only (IC regions, BC inflow states): the model branch is outside every kernel and hot loop.
   integer(I4P), intent(in)  :: model   !< Physical model id.
   real(R8P),    intent(in)  :: gamma   !< Specific heats ratio.
   real(R8P),    intent(in)  :: prim(8) !< Primitive state (first 5 used by Euler).
   real(R8P),    intent(out) :: q(:)    !< Conservative variables [nv].

   select case(model)
   case(MODEL_EULER)
      call primitive_to_conservative(gamma=gamma, r=prim(1), u=prim(2), v=prim(3), w=prim(4), p=prim(5), q=q)
   case(MODEL_EULER_DE)
      call primitive_to_conservative(gamma=gamma, r=prim(1), u=prim(2), v=prim(3), w=prim(4), p=prim(5), q=q(1:NV_EULER))
      q(IQ_S) = prim(5) * prim(1)**(1._R8P - gamma)
   case(MODEL_EULER_PS) ! the scalars are the caller's (the IC regions, the inflow states): zero here
      call primitive_to_conservative(gamma=gamma, r=prim(1), u=prim(2), v=prim(3), w=prim(4), p=prim(5), q=q(1:NV_EULER))
      q(NV_EULER+1:) = 0._R8P
   case(MODEL_EULER_DE_PS)
      call primitive_to_conservative(gamma=gamma, r=prim(1), u=prim(2), v=prim(3), w=prim(4), p=prim(5), q=q(1:NV_EULER))
      q(IQ_S) = prim(5) * prim(1)**(1._R8P - gamma)
      q(IQ_S+1:) = 0._R8P
   case(MODEL_EULER_ISO) ! issue #100: the pressure is not a state variable, prim(5) is ignored, the inert slot 0
      call iso_primitive_to_conservative(gamma=gamma, r=prim(1), u=prim(2), v=prim(3), w=prim(4), p=prim(5), q=q)
   case(MODEL_MHD)
      call mhd_primitive_to_conservative(gamma=gamma, r=prim(1), u=prim(2), v=prim(3), w=prim(4), p=prim(5), &
                                         bx=prim(6), by=prim(7), bz=prim(8), q=q)
   case(MODEL_MHD_GLM, MODEL_MHD_EGLM) ! psi = 0: the EGLM energy psi^2 / 2 vanishes
      call mhd_primitive_to_conservative(gamma=gamma, r=prim(1), u=prim(2), v=prim(3), w=prim(4), p=prim(5), &
                                         bx=prim(6), by=prim(7), bz=prim(8), q=q)
      q(IQ_PSI) = 0._R8P
   endselect
   endsubroutine primitive_state_to_conservative
endmodule adam_flume_physics_object
