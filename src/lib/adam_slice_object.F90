!< ADAM, slice class definition.
module adam_slice_object
!< ADAM, slice class definition.

!< Slice class is a simple object to handle the slicing of ADAM solution fields.

use PENF

implicit none
private
public :: slice_object

type :: slice_object
   !< Slice object.
   character(99)          :: slice_itype           !< Slice interpolation type.
   integer(I4P)           :: slice_save            !< Iteration interval between subsequent data-slice saves.
   integer(I4P)           :: slice_nijk(3)         !< Slice number of points.
   real(R8P)              :: slice_emin(3)         !< Slice minimum extents.
   real(R8P)              :: slice_emax(3)         !< Slice maximum extents.
   real(R8P), allocatable :: slice_points(:,:,:,:) !< Slice points coordinates [3,ni,nj,nk].
endtype slice_object

endmodule adam_slice_object
