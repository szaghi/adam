!< ADAM, AMR marker class definition.
module adam_amr_marker_object
!< ADAM, AMR marker class definition.

!< AMR marker class is a simple object to handle informations concerning AMR marking.

use PENF

implicit none
private
public :: amr_marker_object

type :: amr_marker_object
   !< AMR marker class definition.
   integer(I4P) :: mode         !< Marker mode.
   integer(I4P) :: solid        !< Solid number.
   real(R8P)    :: delta_fine   !< Fine cell space step.
   real(R8P)    :: delta_coarse !< Coarse cell space step.
   integer(I4P) :: ivar         !< ivar.
   real(R8P)    :: tol          !< Tolerance.
endtype amr_marker_object

endmodule adam_amr_marker_object
