!< ADAM, PRISM input/output options object.
module adam_prism_io_object
!< ADAM, PRISM input/output options object.
!<
!< The `[IO]` options only PRISM consumes (issue #80 P5): the library `io_object` reads the keys every application shares,
!< PRISM reads these, so an input of another application that carries them reports them unread.

use :: adam_ini_library, only : ini_check, ini_get
use :: finer, only : file_ini
use :: penf

implicit none
private
public :: prism_io_object

character(len=2), parameter :: INI_SECTION_NAME='IO' !< INI file section name, shared with the library `io_object`.

type :: prism_io_object
   !< PRISM `[IO]` options.
   logical   :: save_memory_status=.false.     !< Save the memory status every iteration.
   logical   :: save_curl_fields=.false.       !< Save the curl fields.
   logical   :: save_divergence_fields=.false. !< Save the divergence fields.
   ! seam div(B) guard-rail (issue #29, accept-truncation): the 2:1 AMR seam injects an O(h^p) div(B) source that is
   ! refinement-convergent but unbounded in t at fixed h; no local collocated-FD fix exists (E1/E2/E3/E/A all ruled out).
   ! This monitor turns the otherwise silent growth into a visible signal on AMR runs. Disabled by default.
   real(R8P) :: seam_divB_tol=-1._R8P          !< Seam max|div(B)| tolerance; <=0 disables the monitor.
   logical   :: seam_divB_error=.false.        !< If .true., exceeding seam_divB_tol is fatal; else warn-only.
   contains
      procedure, pass(self) :: load_from_file !< Load the options from file.
endtype prism_io_object

contains
   subroutine load_from_file(self, file_parameters)
   !< Load the options from file; every one is optional.
   class(prism_io_object), intent(inout) :: self            !< PRISM IO options.
   type(file_ini),         intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   associate(ISN=>INI_SECTION_NAME)
      call ini_get(file_parameters, ISN, 'save_memory_status', self%save_memory_status, default=.false.)
      call ini_get(file_parameters, ISN, 'save_curl_fields', self%save_curl_fields, default=.false.)
      call ini_get(file_parameters, ISN, 'save_divergence_fields', self%save_divergence_fields, default=.false.)
      call ini_get(file_parameters, ISN, 'seam_divB_tol', self%seam_divB_tol, default=-1._R8P)
      call ini_get(file_parameters, ISN, 'seam_divB_error', self%seam_divB_error, default=.false.)
   endassociate
   call ini_check()
   endsubroutine load_from_file
endmodule adam_prism_io_object
