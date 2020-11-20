!< ADAM, general parameters.
module adam_parameters
!< ADAM, general parameters.

use PENF

implicit none
private
public :: TO_BE_REFINED,   &
          TO_BE_DEREFINED, &
          TO_NOT_TOUCH
public :: delta_neighbor

integer(I4P), parameter :: TO_BE_REFINED=1_I4P    !< Flag for node/block to be refined.
integer(I4P), parameter :: TO_BE_DEREFINED=-1_I4P !< Flag for node/block to be derefined.
integer(I4P), parameter :: TO_NOT_TOUCH=0_I4P     !< Flag for node/block to be untouched.

integer(I4P), parameter :: delta_neighbor(3, 26) = reshape([-1,  0,  0, &! face 1
                                                             1,  0,  0, &! face 2
                                                             0, -1,  0, &! face 3
                                                             0,  1,  0, &! face 4
                                                             0,  0, -1, &! face 5
                                                             0,  0,  1, &! face 6
                                                            -1, -1,  0, &! edge 7
                                                             1, -1,  0, &! edge 8
                                                            -1,  1,  0, &! edge 9
                                                             1,  1,  0, &! edge 10
                                                            -1,  0, -1, &! edge 11
                                                             1,  0, -1, &! edge 12
                                                            -1,  0,  1, &! edge 13
                                                             1,  0,  1, &! edge 14
                                                             0, -1, -1, &! edge 15
                                                             0,  1, -1, &! edge 16
                                                             0, -1,  1, &! edge 17
                                                             0,  1,  1, &! edge 18
                                                            -1, -1, -1, &! corner 19
                                                             1, -1, -1, &! corner 20
                                                            -1,  1, -1, &! corner 21
                                                             1,  1, -1, &! corner 22
                                                            -1, -1,  1, &! corner 23
                                                             1, -1,  1, &! corner 24
                                                            -1,  1,  1, &! corner 25
                                                             1,  1,  1  &! corner 26
                                                           ], [3,26]) !< Neighor map.
endmodule adam_parameters
