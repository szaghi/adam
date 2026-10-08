!< ADAM, FLUME dissipative coefficients: viscosity, heat conduction, Ohmic resistivity (issue #65, M4).
module adam_flume_dissipation_object
!< ADAM, FLUME dissipative coefficients: viscosity, heat conduction, Ohmic resistivity (issue #65, M4).
!<
!< The coefficients are read from `[physics]`, each either as a coefficient or as the dimensionless number it stands
!< for (issue #49, N4, section 4.5); giving two keys of a term is fatal, giving none leaves the term off (ideal):
!<
!< | term | coefficient (code units) | number | coefficient from the number |
!< |---|---|---|---|
!< | viscosity | `viscosity` (dynamic, mu) | `reynolds` | `mu = 1 / Re` |
!< | heat conduction | `conductivity` (k) | `prandtl` | `k = mu cp / Pr` |
!< | Ohmic resistivity (MHD) | `resistivity` (magnetic diffusivity eta) | `magnetic_reynolds`, `lundquist` | `1 / Rm`, `1 / S` |
!<
!< In code units the references are 1, so a number is the reciprocal coefficient. The `[reference]` layer converts the
!< dimensional coefficients (`mu`: length x velocity x density; `k` the same divided by the gas constant; `eta`: length
!< x velocity) and leaves the numbers as they are; `lundquist` needs its Alfvenic preset (the velocity unit is then the
!< Alfven speed and `S` is `Rm`). The heat flux is `-k grad T` with `T = p / (rho R)`.
!<
!< `viscosity_law = power-law` makes the viscosity depend on the temperature, `mu (T / reference_temperature) **
!< viscosity_exponent`; with `prandtl` the conductivity follows it, `k = mu(T) cp / Pr`.
!<
!< The viscous and heat-conduction kernels landed in issue #65 P2 (Euler), the MHD ones land in P3: until then a
!< non-zero coefficient with an MHD model is refused at initialisation by the common object, so no coefficient is ever
!< ignored silently.

! ADAM singleton objects
use :: adam_mpih_global,      only : mpih
! FLUME modules
use :: adam_flume_parameters, only : MODEL_EULER, strip_control
! third party modules
use :: finer,                 only : file_ini
use :: penf,                  only : I4P, R8P, str

implicit none
private
public :: flume_dissipation_object
public :: VISCOSITY_LAW_CONSTANT
public :: VISCOSITY_LAW_POWER

character(len=7), parameter :: INI_SECTION_NAME="physics"        !< INI section of the coefficients.
character(len=8), parameter :: VISCOSITY_LAW_CONSTANT="constant" !< [physics].(viscosity_law): constant mu.
character(len=9), parameter :: VISCOSITY_LAW_POWER="power-law"   !< [physics].(viscosity_law): mu (T / T_ref)^omega.

type :: flume_dissipation_object
   !< FLUME dissipative coefficients, code units.
   real(R8P)                 :: mu=0._R8P                    !< Dynamic viscosity (at the reference temperature).
   real(R8P)                 :: k=0._R8P                     !< Thermal conductivity (with Prandtl: at mu).
   real(R8P)                 :: eta=0._R8P                   !< Magnetic diffusivity (Ohmic resistivity).
   real(R8P)                 :: prandtl=0._R8P               !< Prandtl number, 0 when k is given directly.
   character(:), allocatable :: viscosity_law                !< Viscosity law: constant or power-law.
   real(R8P)                 :: viscosity_exponent=0._R8P    !< Power-law exponent omega.
   real(R8P)                 :: reference_temperature=0._R8P !< Power-law reference temperature.
   logical                   :: has_viscosity=.false.        !< mu > 0.
   logical                   :: has_conduction=.false.       !< k > 0.
   logical                   :: has_resistivity=.false.      !< eta > 0.
   logical                   :: is_active=.false.            !< Any term on.
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: laws           !< Return the temperature laws of the kernels.
      procedure, pass(self) :: load_from_file !< Load the coefficients from file.
endtype flume_dissipation_object

contains
   ! public methods
   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(flume_dissipation_object), intent(in) :: self             !< Dissipation.
   character(len=:), allocatable               :: desc             !< Description.
   character(len=1), parameter                 :: NL=new_line('a') !< New line character.

   desc = mpih%myrankstr//'  viscosity mu:    '//trim(str(self%mu))//' ('//self%viscosity_law//')'//NL
   if (self%viscosity_law == VISCOSITY_LAW_POWER) &
      desc = desc//mpih%myrankstr//'  power law:       exponent '//trim(str(self%viscosity_exponent))//            &
                                   ', reference temperature '//trim(str(self%reference_temperature))//NL
   if (self%prandtl > 0._R8P) then
      desc = desc//mpih%myrankstr//'  conductivity k:  '//trim(str(self%k))//' (mu cp / Pr, Pr '//               &
                                   trim(str(self%prandtl))//')'//NL
   else
      desc = desc//mpih%myrankstr//'  conductivity k:  '//trim(str(self%k))//NL
   endif
   desc = desc//mpih%myrankstr//'  resistivity eta: '//trim(str(self%eta))
   endfunction description

   pure subroutine laws(self, tref, omega_mu, omega_k)
   !< Return the temperature laws as the kernels take them, `mu = mu0 (T / tref)^omega_mu`, `k = k0 (T / tref)^omega_k`:
   !< the constant law is `tref = 1` with zero exponents; under the power law the conductivity follows the viscosity
   !< only when it comes from the Prandtl number.
   class(flume_dissipation_object), intent(in)  :: self              !< Dissipation.
   real(R8P),                       intent(out) :: tref              !< Reference temperature.
   real(R8P),                       intent(out) :: omega_mu, omega_k !< Exponents.

   tref = 1._R8P ; omega_mu = 0._R8P ; omega_k = 0._R8P
   if (self%viscosity_law == VISCOSITY_LAW_POWER) then
      tref     = self%reference_temperature
      omega_mu = self%viscosity_exponent
      if (self%prandtl > 0._R8P) omega_k = self%viscosity_exponent
   endif
   endsubroutine laws

   subroutine load_from_file(self, file_parameters, model, cp)
   !< Load the coefficients: a coefficient or its number per term (two keys of a term are fatal, none is the ideal 0).
   class(flume_dissipation_object), intent(inout) :: self            !< Dissipation.
   type(file_ini),                  intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   integer(I4P),                    intent(in)    :: model           !< Physical model id.
   real(R8P),                       intent(in)    :: cp              !< Specific heat at constant pressure, code units.
   character(999)                                 :: buff            !< Option value buffer.
   integer(I4P)                                   :: error           !< Error status.

   ! viscosity: viscosity or reynolds
   call pair(coefficient='viscosity', numbers=['reynolds'], val=self%mu)
   ! heat conduction: conductivity or prandtl (the latter scales the viscosity)
   call pair(coefficient='conductivity', numbers=['prandtl'], val=self%k, number=self%prandtl)
   if (self%prandtl > 0._R8P) then
      if (.not.(self%mu > 0._R8P)) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(prandtl) needs a viscosity '// &
                                                            '(viscosity or reynolds)')
      self%k = self%mu * cp / self%prandtl
   endif
   ! Ohmic resistivity: resistivity, magnetic_reynolds or lundquist (MHD only)
   call pair(coefficient='resistivity', numbers=['magnetic_reynolds', 'lundquist        '], val=self%eta)
   if (self%eta > 0._R8P .and. model == MODEL_EULER) &
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] resistivity (resistivity, magnetic_reynolds, lundquist) '// &
                               'needs [physics].(physical_model)=mhd-ideal')
   ! viscosity law
   self%viscosity_law = VISCOSITY_LAW_CONSTANT
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='viscosity_law', val=buff, error=error)
   if (error <= 0) self%viscosity_law = trim(adjustl(strip_control(buff)))
   select case(self%viscosity_law)
   case(VISCOSITY_LAW_CONSTANT)
      if (has('viscosity_exponent') .or. has('reference_temperature')) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(viscosity_exponent, reference_temperature) need '// &
                                  '[physics].(viscosity_law)='//VISCOSITY_LAW_POWER)
   case(VISCOSITY_LAW_POWER)
      if (.not.(self%mu > 0._R8P)) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(viscosity_law)='// &
                                                            VISCOSITY_LAW_POWER//' needs a viscosity')
      call required(key='viscosity_exponent', val=self%viscosity_exponent)
      ! the diffusive time step bound evaluates the laws at the largest temperature a face can see (issue #65, P2)
      if (self%viscosity_exponent < 0._R8P) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(viscosity_exponent) must not be negative, got '// &
                                  trim(str(self%viscosity_exponent)))
      call required(key='reference_temperature', val=self%reference_temperature)
      if (.not.(self%reference_temperature > 0._R8P)) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(reference_temperature) must be positive')
   case default
      call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(viscosity_law) "'//self%viscosity_law// &
                               '"; expected one of '//VISCOSITY_LAW_CONSTANT//', '//VISCOSITY_LAW_POWER)
   endselect
   self%has_viscosity   = self%mu > 0._R8P
   self%has_conduction  = self%k > 0._R8P
   self%has_resistivity = self%eta > 0._R8P
   self%is_active       = self%has_viscosity .or. self%has_conduction .or. self%has_resistivity
   contains
      subroutine pair(coefficient, numbers, val, number)
      !< Load a term given as `coefficient` (>= 0) or as one of `numbers` (> 0, the coefficient is its reciprocal):
      !< more than one key of the term is fatal, none leaves the coefficient 0.
      character(*),        intent(in)  :: coefficient !< Coefficient key.
      character(*),        intent(in)  :: numbers(:)  !< Number keys.
      real(R8P),           intent(out) :: val         !< Coefficient.
      real(R8P), optional, intent(out) :: number      !< The number given, 0 if none.
      character(:), allocatable        :: given       !< Keys given, for the message.
      real(R8P)                        :: x           !< Value read.
      integer(I4P)                     :: n, i        !< Keys given, counter.

      val = 0._R8P ; n = 0_I4P ; given = ''
      if (present(number)) number = 0._R8P
      if (has(coefficient)) then
         call required(key=coefficient, val=val)
         if (val < 0._R8P) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].('//coefficient//') must not be '// &
                                                    'negative, got '//trim(str(val)))
         n = n + 1_I4P ; given = coefficient
      endif
      do i=1, size(numbers)
         if (.not.has(trim(numbers(i)))) cycle
         call required(key=trim(numbers(i)), val=x)
         if (.not.(x > 0._R8P)) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].('//trim(numbers(i))// &
                                                         ') must be positive, got '//trim(str(x)))
         val = 1._R8P / x
         if (present(number)) number = x
         n = n + 1_I4P ; given = given//' '//trim(numbers(i))
      enddo
      if (n > 1_I4P) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] gives the same term twice ('// &
                                              trim(adjustl(given))//'): give the coefficient or the number, not both')
      endsubroutine pair

      subroutine required(key, val)
      !< Load a required real key.
      character(*), intent(in)  :: key   !< Key.
      real(R8P),    intent(out) :: val   !< Value.
      integer(I4P)              :: error !< Error status.

      call file_parameters%get(section_name=INI_SECTION_NAME, option_name=key, val=val, error=error)
      if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].('//key//')')
      endsubroutine required

      function has(key) result(is_present)
      !< Return true if the key is in the section.
      character(*), intent(in) :: key        !< Key.
      logical                  :: is_present !< Presence.
      character(999)           :: buff       !< Option value buffer.
      integer(I4P)             :: error      !< Error status.

      call file_parameters%get(section_name=INI_SECTION_NAME, option_name=key, val=buff, error=error)
      is_present = error <= 0
      endfunction has
   endsubroutine load_from_file
endmodule adam_flume_dissipation_object
