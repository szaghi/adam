!< ADAM, FLUME boundary conditions class definition.
module adam_flume_bc_object
!< ADAM, FLUME boundary conditions class definition.
!<
!< `periodic` is the library `BC_PERIODIC`, so the tree builds true periodic neighbors and the ghost exchange fills
!< periodic ghosts across blocks and ranks (verified by `src/tests/amr/test_periodic_ghost`); the other kinds are
!< filled by the backends on the boundary crown maps.

! ADAM classes, libraries, parameters
use :: adam_parameters,           only : BC_PERIODIC
! ADAM singleton objects
use :: adam_mpih_global,          only : mpih
! FLUME modules
use :: adam_flume_parameters,     only : IQ_BX, IQ_RU, MODEL_EULER, MODEL_MHD, MODEL_MHD_GLM, strip_control
use :: adam_flume_physics_object, only : flume_physics_object, primitive_state_to_conservative
! third party modules
use :: finer,                     only : file_ini
use :: penf,                      only : I4P, R8P, str

implicit none
private
public :: flume_bc_object
public :: BC_EXTRAPOLATION
public :: BC_INFLOW
public :: BC_WALL_INVISCID
public :: BC_PERIODIC

integer(I4P), parameter :: BC_EXTRAPOLATION = 1_I4P !< Zeroth-order extrapolation.
integer(I4P), parameter :: BC_INFLOW        = 2_I4P !< Prescribed state.
integer(I4P), parameter :: BC_WALL_INVISCID = 3_I4P !< Inviscid (slip) wall: mirror with normal momentum negated.

character(len=13), parameter :: BC_EXTRAPOLATION_STR="extrapolation"          !< Accepted spelling of BC_EXTRAPOLATION.
character(len=6),  parameter :: BC_INFLOW_STR="inflow"                        !< Accepted spelling of BC_INFLOW.
character(len=13), parameter :: BC_WALL_INVISCID_STR="wall-inviscid"          !< Accepted spelling of BC_WALL_INVISCID.
character(len=8),  parameter :: BC_PERIODIC_STR="periodic"                    !< Accepted spelling of BC_PERIODIC.
character(len=2),  parameter :: INFLOW_KEY(8)=['r ', 'u ', 'v ', 'w ', &
                                               'p ', 'bx', 'by', 'bz']       !< Inflow primitive state keys (MHD: all 8).
character(len=8),  parameter :: SECTION_NAME(6)=['bc_x_min', 'bc_x_max', &
                                                 'bc_y_min', 'bc_y_max', &
                                                 'bc_z_min', 'bc_z_max']      !< INI section names of the 6 faces.

type :: flume_bc_object
   !< FLUME boundary conditions class definition.
   integer(I4P)           :: bc_type(6)=0_I4P !< Boundary condition type of each face.
   real(R8P), allocatable :: q_inflow(:,:)    !< Conservative inflow state of each face [nv, 6].
   real(R8P), allocatable :: wall_sign(:,:)   !< Wall mirror sign of each variable per direction [nv, 3] (+1 or -1).
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize boundary conditions.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype flume_bc_object

contains
   ! public methods
   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(flume_bc_object), intent(in) :: self             !< Boundary conditions.
   character(len=:), allocatable      :: desc             !< Description.
   character(len=1), parameter        :: NL=new_line('a') !< New line character.
   integer(I4P)                       :: f                !< Counter.

   desc = mpih%myrankstr//'Boundary conditions main data'
   do f=1, 6
      desc = desc//NL//mpih%myrankstr//'  '//SECTION_NAME(f)//': '//trim(str(self%bc_type(f)))
   enddo
   endfunction description

   subroutine initialize(self, file_parameters, physics)
   !< Initialize boundary conditions.
   class(flume_bc_object),     intent(inout) :: self            !< Boundary conditions.
   type(file_ini),             intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   type(flume_physics_object), intent(in)    :: physics         !< Physics (for the inflow state conversion).

   print '(A)', mpih%myrankstr//'flume_bc_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters, physics=physics)
   print '(A)', self%description()
   print '(A)', mpih%myrankstr//'flume_bc_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, physics)
   !< Load config from file; `type` is required on every face and an unknown value is fatal.
   class(flume_bc_object),     intent(inout) :: self            !< Boundary conditions.
   type(file_ini),             intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   type(flume_physics_object), intent(in)    :: physics         !< Physics (for the inflow state conversion).
   character(999)                            :: buff            !< Option value buffer.
   character(:), allocatable                 :: bc_str          !< BC type string.
   real(R8P)                                 :: prim(8)         !< Inflow primitive state (first nprim used).
   integer(I4P)                              :: nprim           !< Primitive keys number: 5 (Euler) or 8 (MHD).
   integer(I4P)                              :: error           !< Error status.
   integer(I4P)                              :: d, f, k         !< Counters.

   if (allocated(self%q_inflow)) deallocate(self%q_inflow)
   if (allocated(self%wall_sign)) deallocate(self%wall_sign)
   allocate(self%q_inflow(physics%nv,6), self%wall_sign(physics%nv,3))
   self%q_inflow = 0._R8P
   ! wall rule (issue #41, section 3.6): the mirror state negates the wall-normal momentum; MHD (perfectly conducting
   ! reflecting wall, D-12) also negates the wall-normal magnetic field, psi stays even
   self%wall_sign = 1._R8P
   select case(physics%model)
   case(MODEL_EULER)
      nprim = 5_I4P
      do d=1, 3
         self%wall_sign(IQ_RU+d-1,d) = -1._R8P
      enddo
   case(MODEL_MHD, MODEL_MHD_GLM)
      nprim = 8_I4P
      do d=1, 3
         self%wall_sign(IQ_RU+d-1,d) = -1._R8P
         self%wall_sign(IQ_BX+d-1,d) = -1._R8P
      enddo
   case default
      nprim = 0_I4P
      call mpih%error_stop(msg=': no boundary conditions for physical model "'//physics%physical_model//'"')
   endselect
   prim = 0._R8P
   do f=1, 6
      call file_parameters%get(section_name=SECTION_NAME(f), option_name='type', val=buff, error=error)
      if (error > 0) call mpih%error_stop(msg=': failed to load ['//SECTION_NAME(f)//'].(type)')
      bc_str = trim(adjustl(strip_control(buff)))
      select case(bc_str)
      case(BC_EXTRAPOLATION_STR)
         self%bc_type(f) = BC_EXTRAPOLATION
      case(BC_INFLOW_STR)
         self%bc_type(f) = BC_INFLOW
         do k=1, nprim
            call file_parameters%get(section_name=SECTION_NAME(f), option_name=trim(INFLOW_KEY(k)), val=prim(k), &
                                     error=error)
            if (error > 0) &
               call mpih%error_stop(msg=': failed to load ['//SECTION_NAME(f)//'].('//trim(INFLOW_KEY(k))//')')
         enddo
         call primitive_state_to_conservative(model=physics%model, gamma=physics%gamma, prim=prim, q=self%q_inflow(:,f))
      case(BC_WALL_INVISCID_STR)
         self%bc_type(f) = BC_WALL_INVISCID
      case(BC_PERIODIC_STR)
         self%bc_type(f) = BC_PERIODIC
      case default
         call mpih%error_stop(msg=': unknown ['//SECTION_NAME(f)//'].(type) "'//bc_str//'"; expected one of '// &
                                  BC_EXTRAPOLATION_STR//', '//BC_INFLOW_STR//', '//BC_WALL_INVISCID_STR//', '// &
                                  BC_PERIODIC_STR)
      endselect
   enddo
   do f=1, 5, 2
      if ((self%bc_type(f) == BC_PERIODIC) .neqv. (self%bc_type(f+1) == BC_PERIODIC)) &
         call mpih%error_stop(msg=': ['//SECTION_NAME(f)//'] and ['//SECTION_NAME(f+1)//'] must be both periodic or neither')
   enddo
   endsubroutine load_from_file
endmodule adam_flume_bc_object
