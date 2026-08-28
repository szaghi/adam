!< ADAM, PRISM magnetic-field-at-domain-center diagnostic object.
module adam_prism_magnetic_field_at_center_domain_object
!< ADAM, PRISM magnetic-field-at-domain-center diagnostic object.

use :: adam_mpih_global, only : mpih
use :: finer
use :: penf

implicit none
private
public :: prism_magnetic_field_at_center_domain_object

character(len=31), parameter :: MFCD_SECTION_NAME="magnetic_field_at_center_domain"

type :: prism_magnetic_field_at_center_domain_object
   !< PRISM domain-center magnetic-field diagnostic configuration and history writer.
   logical   :: do_save_history = .false.             !< Enable history output.
   integer(I4P) :: history_save = 10_I4P              !< History output cadence.
   integer(I4P) :: history_unit = -1_I4P              !< History file unit.
   character(:), allocatable :: output_basename       !< Basename of output files.
   real(R8P) :: center(3) = [0._R8P, 0._R8P, 0._R8P]  !< Geometrical domain center.
   real(R8P) :: sample_point(3) = [0._R8P, 0._R8P, 0._R8P] !< Cell center used for sampling.
   real(R8P) :: magnetic_field(3) = [0._R8P, 0._R8P, 0._R8P] !< B at sample point.
   real(R8P) :: distance = huge(1._R8P)               !< Distance from sample point to center.
   contains
      procedure, pass(self) :: initialize
      procedure, pass(self) :: load_from_file
      procedure, pass(self) :: save_history
endtype prism_magnetic_field_at_center_domain_object

contains

   subroutine initialize(self, file_parameters, output_basename, verbose)
   !< Initialize the domain-center magnetic-field diagnostic.
   class(prism_magnetic_field_at_center_domain_object), intent(inout)        :: self
   type(file_ini),                                      intent(in)           :: file_parameters
   character(*),                                        intent(in)           :: output_basename
   logical,                                             intent(in), optional :: verbose
   logical                                                                   :: verbose_

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   self%output_basename = trim(output_basename)
   call self%load_from_file(file_parameters=file_parameters)
   if (verbose_ .and. self%do_save_history) then
      print '(A)', mpih%myrankstr//'prism_magnetic_field_at_center_domain_object initialized'
   endif
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters)
   !< Load domain-center magnetic-field config from the shared input file.
   class(prism_magnetic_field_at_center_domain_object), intent(inout) :: self
   type(file_ini),                                      intent(in)    :: file_parameters
   integer(I4P)                                                       :: error

   call file_parameters%get(section_name=MFCD_SECTION_NAME, option_name='save_history', &
                            val=self%do_save_history, error=error)
   if (error > 0) self%do_save_history = .false.
   if (.not. self%do_save_history) return

   call file_parameters%get(section_name=MFCD_SECTION_NAME, option_name='history_save', &
                            val=self%history_save, error=error)
   if (error > 0) self%history_save = 10_I4P
   endsubroutine load_from_file

   subroutine save_history(self, it, time, blocks_number, is_to_open, is_to_close)
   !< Save domain-center magnetic-field history.
   class(prism_magnetic_field_at_center_domain_object), intent(inout)        :: self
   integer(I4P),                                      intent(in)           :: it
   real(R8P),                                         intent(in)           :: time
   integer(I4P),                                      intent(in)           :: blocks_number
   logical,                                           intent(in), optional :: is_to_open
   logical,                                           intent(in), optional :: is_to_close
   logical                                                                 :: is_to_open_
   logical                                                                 :: is_to_close_

   if (.not. self%do_save_history) return
   if (mpih%myrank /= 0) return

   is_to_open_ = .false. ; if (present(is_to_open)) is_to_open_ = is_to_open
   is_to_close_ = .false. ; if (present(is_to_close)) is_to_close_ = is_to_close
   if (is_to_open_) then
      open(newunit=self%history_unit, file=self%output_basename//'-magnetic_field_at_center_domain.dat')
      write(self%history_unit,'(A)') &
            '%VARIABLES="it" "blocks_number" "time" "B_x [T]" "B_y [T]" "B_z [T]" ' // &
            '"x_sample" "y_sample" "z_sample" "center_sample_distance"'
   endif
   write(self%history_unit, '(A)') trim(str(it                         ))//' '//&
                                   trim(str(blocks_number              ))//' '//&
                                   trim(str(time                       ))//' '//&
                                   trim(str(self%magnetic_field(1)     ))//' '//&
                                   trim(str(self%magnetic_field(2)     ))//' '//&
                                   trim(str(self%magnetic_field(3)     ))//' '//&
                                   trim(str(self%sample_point(1)       ))//' '//&
                                   trim(str(self%sample_point(2)       ))//' '//&
                                   trim(str(self%sample_point(3)       ))//' '//&
                                   trim(str(self%distance              ))
   if (is_to_close_) close(self%history_unit)
   endsubroutine save_history

endmodule adam_prism_magnetic_field_at_center_domain_object
