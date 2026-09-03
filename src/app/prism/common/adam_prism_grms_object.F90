!< ADAM, PRISM Grms diagnostic object.
module adam_prism_grms_object
!< ADAM, PRISM Grms diagnostic object.

use :: adam_mpih_global, only : mpih
use :: finer
use :: penf

implicit none
private
public :: prism_grms_object

character(len=4), parameter :: GRMS_SECTION_NAME="grms"

type :: prism_grms_object
   !< PRISM Grms diagnostic configuration and history writer.
   logical   :: do_save_history = .false.          !< Enable Grms history output.
   integer(I4P) :: history_save = 10_I4P           !< Grms history output cadence.
   integer(I4P) :: history_unit = -1_I4P           !< Grms history file unit.
   logical   :: use_cylindrical_region = .false.   !< Restrict the diagnostic to an input cylinder.
   real(R8P) :: center(3) = [0._R8P,0._R8P,0._R8P] !< Cylinder center.
   real(R8P) :: axis(3) = [0._R8P,0._R8P,1._R8P]   !< Cylinder axis, normalized at load time.
   real(R8P) :: radius = -1._R8P                   !< Cylinder radius.
   real(R8P) :: length = -1._R8P                   !< Cylinder length.
   character(:), allocatable :: output_basename    !< Basename of output files.
   real(R8P)    :: grms_domain_B = 0._R8P                  !< Grms over the selected domain.
   real(R8P)    :: grms_3db_B = 0._R8P                     !< Grms over the -3 dB subset of the selected domain.
   real(R8P)    :: grms_domain_B_normalized = 0._R8P       !< Grms over the selected domain, normalized by reference_B.
   real(R8P)    :: grms_3db_B_normalized = 0._R8P          !< Grms over the -3 dB subset, normalized by reference_B.
   real(R8P)    :: absdiff_domain_B = 0._R8P               !< Mean abs(B_amp-reference_B) over the selected domain.
   real(R8P)    :: absdiff_3db_B = 0._R8P                  !< Mean abs(B_amp-reference_B) over the -3 dB subset.
   real(R8P)    :: absdiff_domain_B_normalized = 0._R8P    !< absdiff_domain_B normalized by reference_B.
   real(R8P)    :: absdiff_3db_B_normalized = 0._R8P       !< absdiff_3db_B normalized by reference_B.
   real(R8P)    :: reference_B = 0._R8P                    !< Reference rotating magnetic field amplitude.
   real(R8P)    :: threshold_B = 0._R8P                    !< -3 dB threshold amplitude.
   real(R8P)    :: domain_measure = 0._R8P                 !< Measure of the selected domain.
   real(R8P)    :: measure_3db = 0._R8P                    !< Measure of the -3 dB subset.
   integer(I8P) :: domain_cells_number = 0_I8P             !< Number of cells in the selected domain.
   integer(I8P) :: cells_number_3db = 0_I8P                !< Number of cells in the -3 dB subset.
   contains
      procedure, pass(self) :: initialize
      procedure, pass(self) :: load_from_file
      procedure, pass(self) :: save_history
endtype prism_grms_object

contains

   subroutine initialize(self, file_parameters, output_basename, verbose)
   !< Initialize the Grms diagnostic.
   class(prism_grms_object), intent(inout)        :: self
   type(file_ini),           intent(in)           :: file_parameters
   character(*),             intent(in)           :: output_basename
   logical,                  intent(in), optional :: verbose
   logical                                        :: verbose_

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   self%output_basename = trim(output_basename)
   call self%load_from_file(file_parameters=file_parameters)
   if (verbose_ .and. self%do_save_history) then
      print '(A)', mpih%myrankstr//'prism_grms_object initialized'
   endif
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load Grms config from the shared input file.
   class(prism_grms_object), intent(inout)        :: self
   type(file_ini),           intent(in)           :: file_parameters
   logical,                  intent(in), optional :: go_on_fail
   logical                                        :: go_on_fail_
   character(99)                                  :: buff_c
   integer(I4P)                                   :: error
   real(R8P)                                      :: axis_norm

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=GRMS_SECTION_NAME, option_name='save_history', &
                            val=self%do_save_history, error=error)
   if (error > 0) self%do_save_history = .false.
   call file_parameters%get(section_name=GRMS_SECTION_NAME, option_name='history_save', &
                            val=self%history_save, error=error)

   if (.not. self%do_save_history) return

   call file_parameters%get(section_name=GRMS_SECTION_NAME, option_name='region_type', val=buff_c, error=error)
   if (error > 0) then
      self%use_cylindrical_region = .false.
      return
   endif

   select case(trim(adjustl(buff_c)))
   case('all')
      self%use_cylindrical_region = .false.
   case('cylinder')
      self%use_cylindrical_region = .true.
      call file_parameters%get(section_name=GRMS_SECTION_NAME, option_name='center_x', val=self%center(1), error=error)
      if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load [grms].(center_x)')
      call file_parameters%get(section_name=GRMS_SECTION_NAME, option_name='center_y', val=self%center(2), error=error)
      if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load [grms].(center_y)')
      call file_parameters%get(section_name=GRMS_SECTION_NAME, option_name='center_z', val=self%center(3), error=error)
      if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load [grms].(center_z)')
      call file_parameters%get(section_name=GRMS_SECTION_NAME, option_name='axis_x', val=self%axis(1), error=error)
      if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load [grms].(axis_x)')
      call file_parameters%get(section_name=GRMS_SECTION_NAME, option_name='axis_y', val=self%axis(2), error=error)
      if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load [grms].(axis_y)')
      call file_parameters%get(section_name=GRMS_SECTION_NAME, option_name='axis_z', val=self%axis(3), error=error)
      if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load [grms].(axis_z)')
      call file_parameters%get(section_name=GRMS_SECTION_NAME, option_name='radius', val=self%radius, error=error)
      if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load [grms].(radius)')
      call file_parameters%get(section_name=GRMS_SECTION_NAME, option_name='length', val=self%length, error=error)
      if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load [grms].(length)')
      axis_norm = sqrt(sum(self%axis**2))
      if (axis_norm <= 0._R8P) call mpih%error_stop(msg=': [grms].axis must be non-zero')
      if (self%radius <= 0._R8P) call mpih%error_stop(msg=': [grms].radius must be positive')
      if (self%length <= 0._R8P) call mpih%error_stop(msg=': [grms].length must be positive')
      self%axis = self%axis / axis_norm
   case default
      call mpih%error_stop(msg=': [grms].region_type must be "all" or "cylinder"')
   endselect
   endsubroutine load_from_file

   subroutine save_history(self, it, time, blocks_number, is_to_open, is_to_close)
   !< Save Grms history.
   class(prism_grms_object), intent(inout)        :: self
   integer(I4P),             intent(in)           :: it
   real(R8P),                intent(in)           :: time
   integer(I4P),             intent(in)           :: blocks_number
   logical,                  intent(in), optional :: is_to_open
   logical,                  intent(in), optional :: is_to_close
   logical                                        :: is_to_open_
   logical                                        :: is_to_close_

   if (.not. self%do_save_history) return
   if (mpih%myrank /= 0) return

   is_to_open_ = .false. ; if (present(is_to_open)) is_to_open_ = is_to_open
   is_to_close_ = .false. ; if (present(is_to_close)) is_to_close_ = is_to_close
   if (is_to_open_) then
      open(newunit=self%history_unit, file=self%output_basename//'-grms_history.dat')
      write(self%history_unit,'(A)') &
            '%VARIABLES="it" "blocks_number" "time" "G_rms_domain [T/m]" "G_rms_-3dB [T/m]" ' // &
            '"G_rms_domain/B_ref [1/m]" "G_rms_-3dB/B_ref [1/m]" ' // &
            '"mean_abs(B-B_ref)_domain [T]" "mean_abs(B-B_ref)_-3dB [T]" ' // &
            '"mean_abs(B-B_ref)_domain/B_ref" "mean_abs(B-B_ref)_-3dB/B_ref" "B_ref [T]" ' // &
            '"B_-3dB [T]" "measure_domain" "cells_domain" "measure_-3dB" "cells_-3dB"'
   endif
   write(self%history_unit, '(A)') trim(str(it                     ))//' '//&
                                   trim(str(blocks_number          ))//' '//&
                                   trim(str(time                   ))//' '//&
                                   trim(str(self%grms_domain_B     ))//' '//&
                                   trim(str(self%grms_3db_B        ))//' '//&
                                   trim(str(self%grms_domain_B_normalized   ))//' '//&
                                   trim(str(self%grms_3db_B_normalized      ))//' '//&
                                   trim(str(self%absdiff_domain_B           ))//' '//&
                                   trim(str(self%absdiff_3db_B              ))//' '//&
                                   trim(str(self%absdiff_domain_B_normalized))//' '//&
                                   trim(str(self%absdiff_3db_B_normalized   ))//' '//&
                                   trim(str(self%reference_B                ))//' '//&
                                   trim(str(self%threshold_B                ))//' '//&
                                   trim(str(self%domain_measure             ))//' '//&
                                   trim(str(self%domain_cells_number        ))//' '//&
                                   trim(str(self%measure_3db                ))//' '//&
                                   trim(str(self%cells_number_3db           ))
   if (is_to_close_) close(self%history_unit)
   endsubroutine save_history

endmodule adam_prism_grms_object
