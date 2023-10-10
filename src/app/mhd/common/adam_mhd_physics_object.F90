!< ADAM, MHD fluid physics class definition, CPU backend.
module adam_mhd_physics_object
!< ADAM, MHD fluid physics class definition, CPU backend.

use adam_mpih_object, only : mpih_object
use adam_mhd_eos_object, only : mhd_eos_object
use finer, only : file_ini
use penf, only : I4P, R8P, str

implicit none
private
public :: mhd_physics_object

character(len=7), parameter :: INI_SECTION_NAME='physics' !< INI file section name containing fluid physics.

! Named indexes of q_aux variables, to be initialized when ns is known.
integer(I4P) :: IR = 2_I4P
integer(I4P) :: IU = 3_I4P
integer(I4P) :: IV = 4_I4P
integer(I4P) :: IW = 5_I4P
integer(I4P) :: IG = 6_I4P
integer(I4P) :: IP = 7_I4P

! conservative variables are arranged as follows:
! q(1): rho
! q(2): rho * u
! q(3): rho * v
! q(4): rho * w
! q(5): Bx
! q(6): By
! q(7): Bz
! q(8): rho * E
! q(9): psi
! auxiliary variables are arranged as follows:
! q_aux(1 ): density
! q_aux(2 ): u
! q_aux(3 ): v
! q_aux(4 ): w
! q_aux(5 ): pressure
! q_aux(6 ): ca Alfven speed
! q_aux(7 ): cf fast magneto sonic speed
! q_aux(8 ): cs slow magneto sonic speed
! q_aux(9 ): Bx
! q_aux(10): By
! q_aux(11): Bz
! q_aux(12): psi gauss unphysical potential for divergence-free of B

type :: mhd_physics_object
   !< mhd fluid physics class definition.
   type(mpih_object)    :: mpih          !< MPI handler.
   integer(I4P)         :: nv=9_I4P      !< Number of variables.
   integer(I4P)         :: nv_aux=12_I4P !< Number of auxiliary variables.
   integer(I4P)         :: np=6_I4P      !< Number of 1D primitive variables.
   type(mhd_eos_object) :: eos           !< Equations of state.
   contains
      ! public methods
      procedure, pass(self) :: conservative2primitive !< Return primitive variables from conservative ones.
      procedure, pass(self) :: description            !< Return pretty-printed object description.
      procedure, pass(self) :: initialize             !< Initialize physics.
      procedure, pass(self) :: load_from_file         !< Load config from file.
      procedure, pass(self) :: primitive2conservative !< Return conservative variables from primitive ones.
endtype mhd_physics_object

contains
   ! public methods
   pure function conservative2primitive(self, conservative) result(primitive)
   !< Return primitive variables (rs, u, v, w, r, p, g) from conservative variables (rs, ru, rv, rw, rE).
   class(mhd_physics_object), intent(in) :: self                   !< Equation of state.
   real(R8P),                   intent(in) :: conservative(self%nv)  !< Conservative variables
   real(R8P)                               :: primitive(self%nv_aux) !< Primitive variables
   real(R8P)                               :: c(1:self%ns)           !< Species concentration.

   associate(ns=>self%ns, eos=>self%eos)
   primitive(1:ns) =  conservative(1:ns)
   primitive(IR)   =  sum(conservative(1:ns)) ; c = primitive(1:ns) / primitive(IR)
   primitive(IU)   =  conservative(ns+1) / primitive(IR)
   primitive(IV)   =  conservative(ns+2) / primitive(IR)
   primitive(IW)   =  conservative(ns+3) / primitive(IR)
   primitive(IG)   =  dot_product(c, eos%cp) / dot_product(c, eos%cv)
   primitive(IP)   = (conservative(ns+4) - 0.5_R8P * primitive(IR) * (primitive(IU)**2+primitive(IV)**2+primitive(IW)**2)) * &
                     (primitive(IG) - 1._R8P)
   endassociate
   endfunction conservative2primitive

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(mhd_physics_object), intent(in) :: self             !< Physics.
   character(len=:), allocatable           :: desc             !< Description.
   character(len=1), parameter             :: NL=new_line('a') !< New line character.
   integer(I4P)                            :: s                !< Counter.

   desc =       self%mpih%myrankstr//'Physics main data:'                                     //NL
   desc = desc//self%mpih%myrankstr//'  mhd_physics_object ns:     '//trim(str(self%ns    ))//NL
   desc = desc//self%mpih%myrankstr//'  mhd_physics_object nv:     '//trim(str(self%nv    ))//NL
   desc = desc//self%mpih%myrankstr//'  mhd_physics_object nv_aux: '//trim(str(self%nv_aux))//NL
   desc = desc//self%mpih%myrankstr//'  mhd_physics_object np:     '//trim(str(self%np    ))//NL
   desc = desc//self%mpih%myrankstr//'  mhd_physics_object IR:     '//trim(str(IR         ))//NL
   desc = desc//self%mpih%myrankstr//'  mhd_physics_object IU:     '//trim(str(IU         ))//NL
   desc = desc//self%mpih%myrankstr//'  mhd_physics_object IV:     '//trim(str(IV         ))//NL
   desc = desc//self%mpih%myrankstr//'  mhd_physics_object IW:     '//trim(str(IW         ))//NL
   desc = desc//self%mpih%myrankstr//'  mhd_physics_object IG:     '//trim(str(IG         ))//NL
   desc = desc//self%mpih%myrankstr//'  mhd_physics_object IP:     '//trim(str(IP         ))
   do s=1, self%ns
      desc = desc//NL//self%mpih%myrankstr//'  mhd_physics_object cp('//trim(str(s,.true.))//'):  '//trim(str(self%eos(s)%cp))
      desc = desc//NL//self%mpih%myrankstr//'  mhd_physics_object cv('//trim(str(s,.true.))//'):  '//trim(str(self%eos(s)%cv))
      desc = desc//NL//self%mpih%myrankstr//'  mhd_physics_object  g('//trim(str(s,.true.))//'):  '//trim(str(self%eos(s)%g ))
      desc = desc//NL//self%mpih%myrankstr//'  mhd_physics_object  R('//trim(str(s,.true.))//'):  '//trim(str(self%eos(s)%R ))
   enddo
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize the equation.
   class(mhd_physics_object), intent(inout) :: self            !< Physics.
   type(file_ini),              intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'mhd_physics_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   self%nv     = self%ns + 4
   self%nv_aux = self%ns + 8
   self%np     = self%ns + 4
   ! initialize named index of q_aux array
   IR = self%ns + 1
   IU = self%ns + 2
   IV = self%ns + 3
   IW = self%ns + 4
   IG = self%ns + 5
   IP = self%ns + 6
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'mhd_physics_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(mhd_physics_object), intent(inout)        :: self            !< Physics.
   type(file_ini),              intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                     intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                           :: go_on_fail_     !< Go on if load fails.
   character(:), allocatable                         :: sname           !< Section name.
   real(R8P)                                         :: cp, cv          !< Constant specific heats.
   integer(I4P)                                      :: s               !< Counter.
   integer(I4P)                                      :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='ns', val=self%ns, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(ns)')
   allocate(self%eos(1:self%ns))
   do s=1, self%ns
      call self%eos(s)%load_from_file(file_parameters=file_parameters, s=s)
   enddo
   endsubroutine load_from_file

   function primitive2conservative(self, primitive) result(conservative)
   !< Return conservative variables (rs, ru, rv, rw, rE) from primitive variables (rs, u, v, w, r, p, g).
   class(mhd_physics_object), intent(in) :: self                   !< Equation of state.
   real(R8P),                   intent(in) :: primitive(self%nv_aux) !< Primitive variables
   real(R8P)                               :: conservative(self%nv)  !< Conservative variables

   associate(ns=>self%ns)
   conservative(1:ns) = primitive(1:ns)
   conservative(ns+1) = primitive(IR) * primitive(IU)
   conservative(ns+2) = primitive(IR) * primitive(IV)
   conservative(ns+3) = primitive(IR) * primitive(IW)
   conservative(ns+4) = primitive(IR) * total_energy(g                = primitive(IG), &
                                                     density          = primitive(IR), &
                                                     pressure         = primitive(IP), &
                                                     velocity_sq_norm = primitive(IU)**2+primitive(IV)**2+primitive(IW)**2)
   endassociate
   endfunction primitive2conservative

   ! non TBP
   elemental function internal_energy(g, density, pressure) result(energy_)
   !< Return specific internal energy.
   real(R8P), intent(in) :: g        !< Specific heats ratio.
   real(R8P), intent(in) :: density  !< Density value.
   real(R8P), intent(in) :: pressure !< Pressure value.
   real(R8P)             :: energy_  !< Energy value.

   energy_ = pressure / ((g - 1._R8P) * density)
   endfunction internal_energy

   elemental function total_energy(g, density, pressure, velocity_sq_norm) result(energy_)
   !< Return total specific energy.
   real(R8P), intent(in) :: g                !< Specific heats ratio.
   real(R8P), intent(in) :: density          !< Density value.
   real(R8P), intent(in) :: pressure         !< Pressure value.
   real(R8P), intent(in) :: velocity_sq_norm !< Velocity vector square norm `||velocity||^2`.
   real(R8P)             :: energy_          !< Total specific energy (per unit of mass).

   energy_ = internal_energy(g=g, density=density, pressure=pressure) + 0.5_R8P * velocity_sq_norm
   endfunction total_energy
endmodule adam_mhd_physics_object
