!< ADAM, PRISM Boundary Conditions class definition, CPU backend.
module adam_prism_bc_object
!< ADAM, PRISM Boundary Conditions class definition, CPU backend.

! ADAM modules
use adam_mpih_object
! third party modules
use finer
use penf

implicit none
private
public :: prism_bc_object
public :: BC_EXTRAPOLATION
public :: BC_fWLayer
public :: BC_Silver_Muller
public :: BC_EXTRAP_DIRICHLET
public :: BC_PERIOD

character(len=8), parameter :: INI_SECTION_NAMES(6)=["bc_x_min", "bc_x_max", &
                                                     "bc_y_min", "bc_y_max", &
                                                     "bc_z_min", "bc_z_max"] !< INI (config) file section name containing BC configs.

integer(I4P), parameter :: BC_EXTRAPOLATION   = 1_I4P !< Extrapolation.
integer(I4P), parameter :: BC_fWLayer         = 2_I4P !< fWLayer BC
integer(I4P), parameter :: BC_Silver_Muller   = 3_I4P !< Silver-Muller BC.
integer(I4P), parameter :: BC_EXTRAP_DIRICHLET= 4_I4P !< Prova
integer(I4P), parameter :: BC_PERIOD          = 5_I4P !< Periodic BC.

type :: prism_bc_object
   !< Boundary Conditions class definition, CPU backend.
   type(mpih_object)      :: mpih       !< MPI handler.
   integer(I4P)           :: bc_type(6) !< Boundary condition type.
   real(R8P), allocatable :: q(:,:)     !< Primitive variables (Dx,Dy,Dz,Bx,By,Bz,Jx,Jy,Jz) at BC.
   contains
      ! public methods
      procedure, pass(self) :: initialize     !< Initialize BC.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype prism_bc_object

contains
   ! public methods
   subroutine initialize(self, file_parameters)
   !< Initialize the equation.
   class(prism_bc_object), intent(inout) :: self            !< BC.
   type(file_ini),         intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'prism_bc_object%initialize start'
   allocate(self%q(9,6))
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%mpih%myrankstr//'prism_bc_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(prism_bc_object), intent(inout)        :: self            !< BC.
   type(file_ini),         intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                      :: go_on_fail_     !< Go on if load fails.
   character(8)                                 :: sname           !< Section name.
   character(99)                                :: buff_c          !< Character buffer.
   integer(I4P)                                 :: b               !< Counter.
   integer(I4P)                                 :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
   do b=1, 6
      sname = INI_SECTION_NAMES(b)
      call file_parameters%get(section_name=sname, option_name='type', val=buff_c, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(type)')
      select case(trim(adjustl(buff_c)))
      case('extrapolation')
         self%bc_type(b) = BC_EXTRAPOLATION
      case('fWLayer')
         self%bc_type(b) = BC_fWLayer
      case('Silver_Muller')
         self%bc_type(b) = BC_Silver_Muller
      case('Extr_Diric')
         self%bc_type(b) = BC_EXTRAP_DIRICHLET
      case('periodic')
         self%bc_type(b) = BC_PERIOD
      endselect
   enddo
   endsubroutine load_from_file
endmodule adam_prism_bc_object
