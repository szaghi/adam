!< ADAM, general parameters.
module adam_parameters
!< ADAM, general parameters.

use PENF

implicit none
private
public :: TO_BE_REFINED,   &
          TO_BE_DEREFINED, &
          TO_NOT_TOUCH

integer(I4P), parameter :: TO_BE_REFINED=1_I4P    !< Flag for node/block to be refined.
integer(I4P), parameter :: TO_BE_DEREFINED=-1_I4P !< Flag for node/block to be derefined.
integer(I4P), parameter :: TO_NOT_TOUCH=0_I4P     !< Flag for node/block to be untouched.
endmodule adam_parameters
