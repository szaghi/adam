!< ADAM, PRISM PIC particle boundary conditions (particle closure at the domain faces).
module adam_prism_pic_bc_object
!< ADAM, PRISM PIC particle boundary conditions (particle closure at the domain faces).
!<
!< The particle closure of a face is derived from the field BC of that face (`[bc_*] type`, no input key of its own), by
!< a total map: every field BC kind maps to a particle closure, an unknown kind is fatal, so a new field BC has to state
!< its particle closure here too.
!<
!<+ `periodic` -> PIC_BC_PERIODIC: a particle leaving the domain across the face re-enters from the opposite face. The
!<  deposits and gathers of a particle near the face need no special treatment: the periodic ghosts are images of the
!<  opposite interior (library periodicity), so a particle beyond the face is located in the ghost halo of the block on
!<  that face, exactly as a particle beyond the blocks of its rank. Its position is wrapped at the end of a time step
!<  (never inside the RK stages, whose positions are combined linearly).
!<+ every other kind -> PIC_BC_LEGACY, with a warning: no particle closure, a particle must stay far from the face (a
!<  particle outside the domain stops the run, a conserving current whose support reaches the face stops the run, the
!<  standard deposit drops what spills into the boundary ghosts).

! ADAM modules
use :: adam_grid_object,      only : grid_object
use :: adam_parameters,       only : BC_PERIODIC
! ADAM singleton objects
use :: adam_mpih_global,      only : mpih
! PRISM modules
use :: adam_prism_bc_object,  only : BC_EXTRAPOLATION, BC_NEUMANN, BC_DIRICHLET, BC_Silver_Muller, BC_radiative, BC_PEC
! third party modules
use :: penf,                  only : I4P, R8P, str

implicit none
private
public :: prism_pic_bc_object
public :: PIC_BC_LEGACY
public :: PIC_BC_PERIODIC

integer(I4P), parameter :: PIC_BC_LEGACY   = 0_I4P !< No particle closure: particles must stay far from the face.
integer(I4P), parameter :: PIC_BC_PERIODIC = 1_I4P !< Periodic particle closure.

character(len=8), parameter :: FACE_NAMES(6) = ["bc_x_min", "bc_x_max", "bc_y_min", "bc_y_max", "bc_z_min", "bc_z_max"]

type :: prism_pic_bc_object
   !< PIC particle boundary conditions.
   integer(I4P) :: pic_bc_type(6) = PIC_BC_LEGACY !< Particle closure of each face.
   logical      :: is_periodic(3) = .false.       !< Periodic particle closure along each axis.
   real(R8P)    :: domain_emin(3) = 0._R8P        !< Domain minimum coordinates.
   real(R8P)    :: domain_emax(3) = 0._R8P        !< Domain maximum coordinates.
   real(R8P)    :: domain_length(3) = 0._R8P      !< Domain lengths.
   contains
      ! public methods
      procedure, pass(self) :: description             !< Return pretty-printed object description.
      procedure, pass(self) :: initialize              !< Initialize the particle closures from the field BCs.
      procedure, pass(self) :: is_outside              !< True if a position is outside the domain along a non-periodic axis.
      procedure, pass(self) :: wrap_periodic_positions !< Bring the particles that crossed a periodic face back in the domain.
endtype prism_pic_bc_object

contains
   ! public methods
   pure function description(self, prefix) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_pic_bc_object), intent(in)           :: self             !< PIC BC object.
   character(*),               intent(in), optional :: prefix           !< Prefixing string.
   character(len=:), allocatable                    :: prefix_          !< Prefixing string, local variable.
   character(len=:), allocatable                    :: desc             !< Description.
   character(len=1), parameter                      :: NL=new_line('a') !< New line character.
   integer(I4P)                                     :: f                !< Counter.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = prefix_//'PIC particle boundary conditions:'
   do f=1, 6
      select case(self%pic_bc_type(f))
      case(PIC_BC_PERIODIC)
         desc = desc//NL//prefix_//'  '//FACE_NAMES(f)//': periodic'
      case default
         desc = desc//NL//prefix_//'  '//FACE_NAMES(f)//': legacy (no particle closure)'
      endselect
   enddo
   endfunction description

   subroutine initialize(self, grid)
   !< Initialize the particle closures from the field BCs of the grid (grid%set_bc_type already called).
   class(prism_pic_bc_object), intent(inout) :: self !< PIC BC object.
   type(grid_object),          intent(in)    :: grid !< Grid (sibling realm component, threaded in).
   integer(I4P)                              :: f    !< Counter.

   do f=1, 6
      select case(grid%bc_type(f))
      case(BC_PERIODIC)
         self%pic_bc_type(f) = PIC_BC_PERIODIC
      case(BC_EXTRAPOLATION, BC_NEUMANN, BC_DIRICHLET, BC_Silver_Muller, BC_radiative, BC_PEC)
         self%pic_bc_type(f) = PIC_BC_LEGACY
         if (mpih%myrank == 0) &
            call mpih%print_message('WARNING: ['//FACE_NAMES(f)//'] has no particle closure yet: particles are assumed to '// &
                                    'stay far from the face (legacy closure)')
      case default
         call mpih%error_stop(msg=': no particle closure for the field BC '//trim(str(grid%bc_type(f)))//' of ['// &
                                  FACE_NAMES(f)//']')
      endselect
   enddo
   do f=1, 3
      self%is_periodic(f) = self%pic_bc_type(2*f-1) == PIC_BC_PERIODIC
   enddo
   self%domain_emin = grid%domain_emin
   self%domain_emax = grid%domain_emax
   self%domain_length = grid%domain_emax - grid%domain_emin
   endsubroutine initialize

   pure function is_outside(self, x) result(outside)
   !< True if a position is outside the domain (emin, emax] along an axis without periodic closure.
   class(prism_pic_bc_object), intent(in) :: self    !< PIC BC object.
   real(R8P),                  intent(in) :: x(3)    !< Position.
   logical                                :: outside !< Outside the domain along a non-periodic axis.

   outside = any((.not.self%is_periodic) .and. (x <= self%domain_emin .or. x > self%domain_emax))
   endfunction is_outside

   subroutine wrap_periodic_positions(self, q_pic, particle_number)
   !< Bring the particles that crossed a periodic face back in the domain (emin, emax], by one domain length (a particle
   !< moves less than a cell in a step). Particles inside the domain are not touched (bitwise).
   class(prism_pic_bc_object), intent(in)    :: self            !< PIC BC object.
   real(R8P),                  intent(inout) :: q_pic(1:,1:)    !< PIC variables [variable, particle].
   integer(I4P),               intent(in)    :: particle_number !< Number of particles.
   integer(I4P)                              :: n, d            !< Counters.

   if (.not.any(self%is_periodic)) return
   do n=1, particle_number
      do d=1, 3
         if (.not.self%is_periodic(d)) cycle
         if (q_pic(d,n) > self%domain_emax(d)) then
            q_pic(d,n) = q_pic(d,n) - self%domain_length(d)
         elseif (q_pic(d,n) <= self%domain_emin(d)) then
            q_pic(d,n) = q_pic(d,n) + self%domain_length(d)
         endif
         if (q_pic(d,n) > self%domain_emax(d) .or. q_pic(d,n) <= self%domain_emin(d)) &
            call mpih%error_stop(msg=': PIC particle '//trim(str(n,.true.))//' moved more than one domain length in a step')
      enddo
   enddo
   endsubroutine wrap_periodic_positions
endmodule adam_prism_pic_bc_object
