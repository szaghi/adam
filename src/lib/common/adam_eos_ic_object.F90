!< ADAM, Equations Of State (EOS) and physics for Ideal Compressible fluids, class definition, common to all backends.
module adam_eos_ic_object
!< ADAM, Equations Of State (EOS) and physics for Ideal Compressible fluids, class definition, common to all backends.

! ADAM modules
use :: adam_global_mpih, only: mpih
! Third party modules
use :: finer, only : file_ini
use :: penf, only : I4P, R8P, str

implicit none
private
public :: eos_ic_object
public :: ic_density
public :: ic_internal_energy
public :: ic_pressure
public :: ic_speed_of_sound
public :: ic_temperature
public :: ic_total_energy
public :: ic_total_entalpy
public :: ic_conservative2primitive
public :: ic_primitive2conservative

character(len=15), parameter :: INI_SECTION_NAME='physics_specie_' !< INI file section name containing fluid eos.

type :: eos_ic_object
   !< Equations Of State (EOS, ideal, compressible fluid) class definition.
   integer(I4P)      :: id    = 0_I4P  !< Fluid specie unique ID.
   real(R8P)         :: cp    = 0._R8P !< Specific heat at constant pressure `cp`.
   real(R8P)         :: cv    = 0._R8P !< Specific heat at constant volume `cv`.
   real(R8P)         :: g     = 0._R8P !< Specific heats ratio `gamma = cp / cv`.
   real(R8P)         :: R     = 0._R8P !< Fluid constant `R = cp - cv`.
   real(R8P)         :: gm1   = 0._R8P !< `gamma - 1`.
   real(R8P)         :: gp1   = 0._R8P !< `gamma + 1`.
   real(R8P)         :: delta = 0._R8P !< `(gamma - 1) / 2`.
   real(R8P)         :: eta   = 0._R8P !< `2 * gamma / (gamma - 1)`.
   real(R8P)         :: mu    = 0._R8P !< Dynamic viscosity.
   real(R8P)         :: kd    = 0._R8P !< Thermal diffusivity.
   real(R8P)         :: dha   = 0._R8P !< Entalpy formation.
   contains
      ! public methods
      procedure, pass(self) :: compute_derivate       !< Compute derivate quantities (from `cp` and `cv`).
      ! procedure, pass(self) :: conservative2primitive !< Return primitive variables from conservative ones.
      procedure, pass(self) :: description            !< Return pretty-printed object description.
      procedure, pass(self) :: destroy                !< Destroy physics.
      procedure, pass(self) :: initialize             !< initialize physics.
      procedure, pass(self) :: load_from_file         !< Load from finer file.
      procedure, pass(self) :: primitive2conservative !< Return conservative variables from primitive ones.
      procedure, pass(self) :: save_into_file         !< Save into finer file.
      ! states functions
      procedure, pass(self) :: density         !< Return density.
      procedure, pass(self) :: internal_energy !< Return specific internal energy.
      procedure, pass(self) :: pressure        !< Return pressure.
      procedure, pass(self) :: speed_of_sound  !< Return speed of sound.
      procedure, pass(self) :: temperature     !< Return temperature.
      procedure, pass(self) :: total_energy    !< Return total specific energy.
      procedure, pass(self) :: total_entalpy   !< Return total specific entalpy.
      ! operators
      generic :: assignment(=) => eos_assign_eos !< Overload `=`.
      ! private methods
      procedure, pass(lhs), private :: eos_assign_eos !< Operator `=`.
endtype eos_ic_object

interface eos_ic_object
   !< Overload [[eos_ic_object]] name with its constructor.
   module procedure eos_ic_instance
endinterface
contains
   ! public methods
   elemental subroutine compute_derivate(self)
   !< Compute derivate quantities (from `cp` and `cv`).
   class(eos_ic_object), intent(inout) :: self !< Fluid physics.

   self%g     = self%cp / self%cv
   self%R     = self%cp - self%cv
   self%gm1   = self%g - 1._R8P
   self%gp1   = self%g + 1._R8P
   self%delta = (self%g - 1._R8P) * 0.5_R8P
   self%eta   = 2._R8P * self%g / (self%g - 1._R8P)
   endsubroutine compute_derivate

   pure function conservative2primitive(self, conservative) result(primitive)
   !< Return primitive variables (r, u, v, w, p) from conservative variables (r, ru, rv, rw, rE).
   class(eos_ic_object), intent(in) :: self            !< Fluid physics.
   real(R8P),            intent(in) :: conservative(5) !< Conservative variables
   real(R8P)                        :: primitive(5)    !< Primitive variables

   primitive(1) =  conservative(1)
   primitive(2) =  conservative(2) / primitive(1)
   primitive(3) =  conservative(3) / primitive(1)
   primitive(4) =  conservative(4) / primitive(1)
   primitive(5) = (conservative(5) - 0.5_R8P * primitive(1) * (primitive(2)**2+primitive(3)**2+primitive(4)**2)) * self%gm1
   endfunction conservative2primitive

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(eos_ic_object), intent(in) :: self             !< Fluid physics.
   character(len=:), allocatable    :: desc             !< Description.
   character(len=1), parameter      :: NL=new_line('a') !< New line character.

   desc =       mpih%myrankstr//'EOS, ideal, compressible fluid properities of specie "'//trim(str(self%id,.true.))//'"'//NL
   desc = desc//mpih%myrankstr//'  cp    (specific heat constant pressure) = '          //trim(str(self%cp       ))     //NL
   desc = desc//mpih%myrankstr//'  cv    (specific heat constant volume)   = '          //trim(str(self%cv       ))     //NL
   desc = desc//mpih%myrankstr//'  g     (cp/cv)                           = '          //trim(str(self%g        ))     //NL
   desc = desc//mpih%myrankstr//'  R     (cp-cv)                           = '          //trim(str(self%R        ))     //NL
   desc = desc//mpih%myrankstr//'  gm1   (g-1)                             = '          //trim(str(self%gm1      ))     //NL
   desc = desc//mpih%myrankstr//'  gp1   (g+1)                             = '          //trim(str(self%gp1      ))     //NL
   desc = desc//mpih%myrankstr//'  delta (gm1/2)                           = '          //trim(str(self%delta    ))     //NL
   desc = desc//mpih%myrankstr//'  eta   (2*g/gm1)                         = '          //trim(str(self%eta      ))     //NL
   desc = desc//mpih%myrankstr//'  mu    (dinamic viscosity)               = '          //trim(str(self%mu       ))     //NL
   desc = desc//mpih%myrankstr//'  kd    (thermal diffusivity)             = '          //trim(str(self%kd       ))     //NL
   desc = desc//mpih%myrankstr//'  dha   (entalphy diffusivity)            = '          //trim(str(self%dha      ))
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy physics.
   class(eos_ic_object), intent(inout) :: self  !< Fluid physics.
   type(eos_ic_object)                 :: fresh !< Fresh fluid physics.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, file_parameters, s, go_on_fail, physics, mu, kd, dha, cp, cv, gam, R, verbose)
   !< Initialize fluid physics.
   class(eos_ic_object), intent(inout)        :: self            !< Fluid physics.
   type(file_ini),       intent(in), optional :: file_parameters !< Simulation parameters ini file handler.
   integer(I4P),         intent(in), optional :: s               !< Fluid specie number.
   logical,              intent(in), optional :: go_on_fail      !< Go on if load fails.
   type(eos_ic_object),  intent(in), optional :: physics         !< Fluid physics.
   real(R8P),            intent(in), optional :: mu              !< Dynamic viscosity.
   real(R8P),            intent(in), optional :: kd              !< Thermal diffusivity.
   real(R8P),            intent(in), optional :: dha             !< Entalpy formation.
   real(R8P),            intent(in), optional :: cp              !< Specific heat at constant pressure `cp` value.
   real(R8P),            intent(in), optional :: cv              !< Specific heat at constant volume `cv` value.
   real(R8P),            intent(in), optional :: gam             !< Specific heats ratio `gamma=cp/cv` value.
   real(R8P),            intent(in), optional :: R               !< Fluid constant `R=cp-cv` value.
   logical,              intent(in), optional :: verbose         !< Flag to activate verbose output.

   self%id = 1_I4P ; if (present(s)) self%id = s
   if     (present(file_parameters)) then
      call self%load_from_file(file_parameters=file_parameters, s=self%id, go_on_fail=go_on_fail)
   elseif (present(physics)) then
      self = physics
   else
      self = eos_ic_object(mu=mu, kd=kd, dha=dha, cp=cp, cv=cv, gam=gam, R=R)
   endif
   if (present(verbose)) then
      if (verbose) print '(A)', self%description()
   endif
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, s, go_on_fail)
   !< Load from file.
   class(eos_ic_object), intent(inout)        :: self            !< Fluid physics.
   type(file_ini),       intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   integer(I4P),         intent(in)           :: s               !< Fluid specie number.
   logical,              intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                    :: go_on_fail_     !< Go on if load fails.
   integer(I4P)                               :: error           !< Error status.
   character(:), allocatable                  :: sname           !< Section name.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
   sname = INI_SECTION_NAME//trim(str(s,.true.))

   call self%destroy
   call file_parameters%get(section_name=sname, option_name='cp', val=self%cp, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(cp)')
   call file_parameters%get(section_name=sname, option_name='cv', val=self%cv, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(cv)')
   call self%compute_derivate
   call file_parameters%get(section_name=sname, option_name='mu', val=self%mu, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(mu)')
   call file_parameters%get(section_name=sname, option_name='kd', val=self%kd, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(kd)')
   call file_parameters%get(section_name=sname, option_name='dha', val=self%dha, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(dha)')
   endsubroutine load_from_file

   pure function primitive2conservative(self, primitive) result(conservative)
   !< Return conservative variables (r, ru, rv, rw, rE) from primitive variables (r, u, v, w, p).
   class(eos_ic_object), intent(in) :: self            !< Fluid physics.
   real(R8P),            intent(in) :: primitive(5)    !< Primitive variables
   real(R8P)                        :: conservative(5) !< Conservative variables

   conservative(1) = primitive(1)
   conservative(2) = primitive(1) * primitive(2)
   conservative(3) = primitive(1) * primitive(3)
   conservative(4) = primitive(1) * primitive(4)
   conservative(5) = primitive(1) * self%total_energy(density          = primitive(1), &
                                                      pressure         = primitive(5), &
                                                      velocity_sq_norm = primitive(2)**2+primitive(3)**2+primitive(4)**2)
   endfunction primitive2conservative

   subroutine save_into_file(self, file_parameters, s)
   !< Save into file.
   class(eos_ic_object), intent(inout) :: self            !< Free conditions.
   type(file_ini),       intent(inout) :: file_parameters !< Simulation parameters ini file handler.
   integer(I4P),         intent(in)    :: s               !< Fluid specie number.
   character(:), allocatable           :: sname           !< Section name.

   sname = INI_SECTION_NAME//trim(str(s,.true.))

   call file_parameters%add(section_name=sname, option_name='cp',  val=self%cp )
   call file_parameters%add(section_name=sname, option_name='cv',  val=self%cv )
   call file_parameters%add(section_name=sname, option_name='mu',  val=self%mu )
   call file_parameters%add(section_name=sname, option_name='kd',  val=self%kd )
   call file_parameters%add(section_name=sname, option_name='dha', val=self%dha)
   endsubroutine save_into_file

   ! states functions
   elemental function density(self, pressure, speed_of_sound) result(density_)
   !< Return density.
   class(eos_ic_object), intent(in) :: self           !< Fluid physics.
   real(R8P),            intent(in) :: pressure       !< Pressure value.
   real(R8P),            intent(in) :: speed_of_sound !< Speed of sound value.
   real(R8P)                        :: density_       !< Density value.

   density_ = ic_density(g=self%g, pressure=pressure, speed_of_sound=speed_of_sound)
   endfunction density

   elemental function internal_energy(self, density, pressure) result(energy_)
   !< Return specific internal energy.
   class(eos_ic_object), intent(in) :: self     !< Fluid physics.
   real(R8P),            intent(in) :: density  !< Density value.
   real(R8P),            intent(in) :: pressure !< Pressure value.
   real(R8P)                        :: energy_  !< Energy value.

   energy_ = ic_internal_energy(gm1=self%gm1, density=density, pressure=pressure)
   endfunction internal_energy

   elemental function pressure(self, density, energy) result(pressure_)
   !< Return pressure.
   class(eos_ic_object), intent(in) :: self      !< Fluid physics.
   real(R8P),            intent(in) :: density   !< Density value.
   real(R8P),            intent(in) :: energy    !< Specific internal energy value.
   real(R8P)                        :: pressure_ !< Pressure value.

   pressure_ = ic_pressure(gm1=self%gm1, density=density, energy=energy)
   endfunction pressure

   elemental function speed_of_sound(self, density, pressure) result(speed_of_sound_)
   !< Return speed of sound.
   class(eos_ic_object), intent(in) :: self            !< Fluid physics.
   real(R8P),            intent(in) :: density         !< Density value.
   real(R8P),            intent(in) :: pressure        !< Pressure value.
   real(R8P)                        :: speed_of_sound_ !< Speed of sound value.

   speed_of_sound_ = ic_speed_of_sound(g=self%g, density=density, pressure=pressure)
   endfunction speed_of_sound

   elemental function temperature(self, density, pressure) result(temperature_)
   !< Return temperature.
   class(eos_ic_object), intent(in) :: self         !< Fluid physics.
   real(R8P),            intent(in) :: density      !< Density value.
   real(R8P),            intent(in) :: pressure     !< Pressure value.
   real(R8P)                        :: temperature_ !< Temperature value.

   temperature_ = ic_temperature(R=self%R, density=density, pressure=pressure)
   endfunction temperature

   elemental function total_energy(self, density, pressure, velocity_sq_norm) result(energy_)
   !< Return total specific energy.
   class(eos_ic_object), intent(in) :: self             !< Fluid physics.
   real(R8P),            intent(in) :: density          !< Density value.
   real(R8P),            intent(in) :: pressure         !< Pressure value.
   real(R8P),            intent(in) :: velocity_sq_norm !< Velocity vector square norm `||velocity||^2`.
   real(R8P)                        :: energy_          !< Total specific energy (per unit of mass).

   energy_ = ic_total_energy(gm1=self%gm1, density=density, pressure=pressure, velocity_sq_norm=velocity_sq_norm)
   endfunction total_energy

   elemental function total_entalpy(self, density, pressure, velocity_sq_norm) result(entalpy_)
   !< Return total specific entalpy.
   class(eos_ic_object), intent(in) :: self             !< Fluid physics.
   real(R8P),            intent(in) :: density          !< Density value.
   real(R8P),            intent(in) :: pressure         !< Pressure value.
   real(R8P),            intent(in) :: velocity_sq_norm !< Velocity vector square norm `||velocity||^2`.
   real(R8P)                        :: entalpy_         !< Total specific entalpy (per unit of mass).

   entalpy_ = ic_total_entalpy(g=self%g, density=density, pressure=pressure, velocity_sq_norm=velocity_sq_norm)
   endfunction total_entalpy

   ! operators
   pure subroutine eos_assign_eos(lhs, rhs)
   !< Operator `=`.
   class(eos_ic_object), intent(inout) :: lhs !< Left hand side.
   type(eos_ic_object),  intent(in)    :: rhs !< Right hand side.

   lhs%cp    = rhs%cp
   lhs%cv    = rhs%cv
   lhs%g     = rhs%g
   lhs%R     = rhs%R
   lhs%delta = rhs%delta
   lhs%eta   = rhs%eta
   lhs%gm1   = rhs%gm1
   lhs%gp1   = rhs%gp1
   lhs%mu    = rhs%mu
   lhs%kd    = rhs%kd
   lhs%dha   = rhs%dha
   endsubroutine eos_assign_eos

   ! public non TBP
   elemental function ic_density(g, pressure, speed_of_sound) result(density_)
   !< Return density for ideal compressible fluid.
   real(R8P), intent(in) :: g              !< Specific heats ratio.
   real(R8P), intent(in) :: pressure       !< Pressure value.
   real(R8P), intent(in) :: speed_of_sound !< Speed of sound value.
   real(R8P)             :: density_       !< Density value.

   density_ = g * pressure / (speed_of_sound * speed_of_sound)
   endfunction ic_density

   elemental function ic_internal_energy(gm1, density, pressure) result(energy_)
   !< Return specific internal energy for ideal compressible fluid.
   real(R8P), intent(in) :: gm1      !< Specific heats ratio minus 1.
   real(R8P), intent(in) :: density  !< Density value.
   real(R8P), intent(in) :: pressure !< Pressure value.
   real(R8P)             :: energy_  !< Energy value.

   energy_ = pressure / (gm1 * density)
   endfunction ic_internal_energy

   elemental function ic_pressure(gm1, density, energy) result(pressure_)
   !< Return pressure for ideal compressible fluid.
   real(R8P), intent(in) :: gm1       !< Specific heats ratio minus 1.
   real(R8P), intent(in) :: density   !< Density value.
   real(R8P), intent(in) :: energy    !< Specific internal energy value.
   real(R8P)             :: pressure_ !< Pressure value.

   pressure_ = density * gm1 * energy
   endfunction ic_pressure

   elemental function ic_speed_of_sound(g, density, pressure) result(speed_of_sound_)
   !< Return speed of sound for ideal compressible fluid.
   real(R8P), intent(in) :: g               !< Specific heats ratio.
   real(R8P), intent(in) :: density         !< Density value.
   real(R8P), intent(in) :: pressure        !< Pressure value.
   real(R8P)             :: speed_of_sound_ !< Speed of sound value.

   speed_of_sound_ = sqrt(g * pressure / density)
   endfunction ic_speed_of_sound

   elemental function ic_temperature(R, density, pressure) result(temperature_)
   !< Return temperature for ideal compressible fluid.
   real(R8P), intent(in) :: R            !< Fluid constant.
   real(R8P), intent(in) :: density      !< Density value.
   real(R8P), intent(in) :: pressure     !< Pressure value.
   real(R8P)             :: temperature_ !< Temperature value.

   temperature_ = pressure / (R * density)
   endfunction ic_temperature

   elemental function ic_total_energy(gm1, density, pressure, velocity_sq_norm) result(energy_)
   !< Return total specific energy for ideal compressible fluid.
   real(R8P), intent(in) :: gm1              !< Specific heats ratio minus 1.
   real(R8P), intent(in) :: density          !< Density value.
   real(R8P), intent(in) :: pressure         !< Pressure value.
   real(R8P), intent(in) :: velocity_sq_norm !< Velocity vector square norm `||velocity||^2`.
   real(R8P)             :: energy_          !< Total specific energy (per unit of mass).

   energy_ = ic_internal_energy(gm1=gm1, density=density, pressure=pressure) + 0.5_R8P * velocity_sq_norm
   endfunction ic_total_energy

   elemental function ic_total_entalpy(g, density, pressure, velocity_sq_norm) result(entalpy_)
   !< Return total specific entalpy for ideal compressible fluid.
   real(R8P), intent(in) :: g                !< Specific heats ratio.
   real(R8P), intent(in) :: density          !< Density value.
   real(R8P), intent(in) :: pressure         !< Pressure value.
   real(R8P), intent(in) :: velocity_sq_norm !< Velocity vector square norm `||velocity||^2`.
   real(R8P)             :: entalpy_         !< Total specific entalpy (per unit of mass).

   entalpy_ = g * pressure/((g-1._R8P) * density) + 0.5_R8P * velocity_sq_norm
   endfunction ic_total_entalpy

   pure function ic_conservative2primitive(gm1, conservative) result(primitive)
   !< Return primitive variables (r, u, v, w, p) from conservative variables (r, ru, rv, rw, rE).
   real(R8P), intent(in) :: gm1             !< Specific heats ratio minus 1.
   real(R8P), intent(in) :: conservative(5) !< Conservative variables
   real(R8P)             :: primitive(5)    !< Primitive variables

   primitive(1) =  conservative(1)
   primitive(2) =  conservative(2) / primitive(1)
   primitive(3) =  conservative(3) / primitive(1)
   primitive(4) =  conservative(4) / primitive(1)
   primitive(5) = (conservative(5) - 0.5_R8P * primitive(1) * &
                  (primitive(2)*primitive(2)+primitive(3)*primitive(3)+primitive(4)*primitive(4))) * gm1
   endfunction ic_conservative2primitive

   pure function ic_primitive2conservative(gm1, primitive) result(conservative)
   !< Return conservative variables (r, ru, rv, rw, rE) from primitive variables (r, u, v, w, p).
   real(R8P), intent(in) :: gm1             !< Specific heats ratio minus 1.
   real(R8P), intent(in) :: primitive(5)    !< Primitive variables
   real(R8P)             :: conservative(5) !< Conservative variables

   conservative(1) = primitive(1)
   conservative(2) = primitive(1) * primitive(2)
   conservative(3) = primitive(1) * primitive(3)
   conservative(4) = primitive(1) * primitive(4)
   conservative(5) = primitive(1) * ic_total_energy(gm1              = gm1,                      &
                                                    density          = primitive(1),             &
                                                    pressure         = primitive(5),             &
                                                    velocity_sq_norm = primitive(2)*primitive(2)+&
                                                                       primitive(3)*primitive(3)+&
                                                                       primitive(4)*primitive(4))
   endfunction ic_primitive2conservative

   ! private non TBP
   elemental function eos_ic_instance(mu, kd, dha, cp, cv, gam, R) result(instance)
   !< Return and instance of [[eos_ic_object]].
   !<
   !< @note This procedure is used for overloading [[eos_ic_object]] name.
   real(R8P), intent(in)           :: mu       !< Dynamic viscosity.
   real(R8P), intent(in)           :: kd       !< Thermal diffusivity.
   real(R8P), intent(in)           :: dha      !< Entalpy formation.
   real(R8P), intent(in), optional :: cp       !< Specific heat at constant pressure `cp` value.
   real(R8P), intent(in), optional :: cv       !< Specific heat at constant volume `cv` value.
   real(R8P), intent(in), optional :: gam      !< Specific heats ratio `gamma=cp/cv` value.
   real(R8P), intent(in), optional :: R        !< Fluid constant `R=cp-cv` value.
   type(eos_ic_object)             :: instance !< Instance of [[eos_ic_object]].

   instance%mu  = mu
   instance%kd  = kd
   instance%dha = dha
   if (present(cp).and.present(cv)) then
      instance%cp = cp
      instance%cv = cv
   elseif (present(gam).and.present(R)) then
      instance%cv = R/(gam - 1._R8P)
      instance%cp = gam * instance%cv
   elseif (present(gam).and.present(cp)) then
      instance%cp = cp
      instance%cv = cp / gam
   elseif (present(gam).and.present(cv)) then
      instance%cp = gam * cv
      instance%cv = cv
   elseif (present(R).and.present(cp)) then
      instance%cp = cp
      instance%cv = cp - R
   elseif (present(R).and.present(cv)) then
      instance%cp = cv + R
      instance%cv = cv
   endif
   call instance%compute_derivate
   endfunction eos_ic_instance
endmodule adam_eos_ic_object
