!< ADAM, AMR markers class definition, CPU backend.
module adam_amr_cpu_object
!< ADAM, AMR markers class definition, CPU backend.

use adam_grid_object, only : grid_object
use adam_mpih_object, only : mpih_object
use finer
use penf

implicit none
private
public :: amr_cpu_object
public :: amr_marker_cpu_object
public :: AMR_GEO
public :: AMR_GRAD

character(len=3), parameter :: INI_SECTION_NAME="amr" !< INI (config) file section name containing AMR markers configs.

integer(I4P), parameter :: AMR_GEO  = 1_I4P !< Geometrical marker.
integer(I4P), parameter :: AMR_GRAD = 2_I4P !< Field gradient marker.

type :: amr_marker_cpu_object
   !< AMR marker object.
   integer(I4P) :: mode=AMR_GRAD !< Marker mode.
   real(R8P)    :: delta_fine    !< Fine cell space step.
   real(R8P)    :: delta_coarse  !< Coarse cell space step.
   integer(I4P) :: field=1_I4P   !< Field array containing the marker variable, 1=q, 2=q_aux.
   integer(I4P) :: ivar=1_I4P    !< ivar.
   real(R8P)    :: tol=0.5_R8P   !< Tolerance.
   integer(I4P) :: solid=1_I4P   !< Solid number.
endtype amr_marker_cpu_object

type :: amr_cpu_object
   !< AMR markers class definition, CPU backend.
   type(mpih_object)                        :: mpih                 !< MPI handler.
   integer(I4P)                             :: iters=5_I4P          !< AMR updates iterations number.
   integer(I4P)                             :: frequency=100_I4P    !< AMR update time step frequency.
   integer(I4P)                             :: markers_number=1_I4P !< AMR number of markers.
   type(amr_marker_cpu_object), allocatable :: markers(:)           !< AMR array of marker objects.
   contains
      ! public methods
      procedure, pass(self) :: initialize     !< Initialize IC.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype amr_cpu_object

contains
   ! public methods
   subroutine initialize(self, file_parameters)
   !< Initialize the equation.
   class(amr_cpu_object), intent(inout) :: self            !< AMR.
   type(file_ini),        intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'amr_cpu_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%mpih%myrankstr//'amr_cpu_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(amr_cpu_object), intent(inout)        :: self            !< AMR.
   type(file_ini),        intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,               intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                     :: go_on_fail_     !< Go on if load fails.
   character(:), allocatable                   :: sname           !< Section name.
   integer(I4P)                                :: i_marker        !< Counter.
   integer(I4P)                                :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='frequency', val=self%frequency, error=error)
   ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//INI_SECTION_NAME//'].(frequency)'
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='iters', val=self%iters, error=error)
   ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//INI_SECTION_NAME//'].(iters)'
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='markers_number', val=self%markers_number, error=error)
   ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//&
                                                ! 'error: failed to load ['//INI_SECTION_NAME//'].(markers_number)'

   allocate(self%markers(self%markers_number))
   do i_marker=1, self%markers_number
      sname = INI_SECTION_NAME//'_marker_'//trim(str(i_marker,.true.))
      call file_parameters%get(section_name=sname, option_name='mode', val=self%markers(i_marker)%mode, error=error)
      ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//sname//'].(mode)'
      call file_parameters%get(section_name=sname, option_name='delta_fine', val=self%markers(i_marker)%delta_fine, error=error)
      ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//sname//'].(delta_fine)'
      call file_parameters%get(section_name=sname, option_name='delta_coarse', val=self%markers(i_marker)%delta_coarse, error=error)
      ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//sname//'].(delta_coarse)'
      select case(self%markers(i_marker)%mode)
      case(AMR_GEO)
         call file_parameters%get(section_name=sname, option_name='solid', val=self%markers(i_marker)%solid, error=error)
         ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//sname//'].(solid)'
      case(AMR_GRAD)
         call file_parameters%get(section_name=sname, option_name='field', val=self%markers(i_marker)%field, error=error)
         ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//sname//'].(field)'
         call file_parameters%get(section_name=sname, option_name='var', val=self%markers(i_marker)%ivar, error=error)
         ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//sname//'].(var)'
         call file_parameters%get(section_name=sname, option_name='tol', val=self%markers(i_marker)%tol, error=error)
         ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//sname//'].(tol)'
      endselect
   enddo
   endsubroutine load_from_file
endmodule adam_amr_cpu_object
