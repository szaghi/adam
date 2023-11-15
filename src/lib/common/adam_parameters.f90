!< ADAM, general parameters.
module adam_parameters
!< ADAM, general parameters.

use PENF

implicit none
save
private
public :: BC_PERIODIC
public :: TO_BE_REFINED,   &
          TO_BE_DEREFINED, &
          TO_NOT_TOUCH
public :: FEC_1_6_ARRAY
public :: FEC_TO_DELTA
public :: DELTA_TO_FEC

integer(I4P), parameter :: BC_PERIODIC = -1_I4P !< Flag (reserved) for periodical boundary conditions.

integer(I4P), parameter :: TO_BE_REFINED=1_I4P    !< Flag for node/block to be refined.
integer(I4P), parameter :: TO_BE_DEREFINED=-1_I4P !< Flag for node/block to be derefined.
integer(I4P), parameter :: TO_NOT_TOUCH=0_I4P     !< Flag for node/block to be untouched.

integer(I4P), parameter :: FEC_1_6_ARRAY(26) = [1, & ! 1
                                                2, & ! 2
                                                3, & ! 3
                                                4, & ! 4
                                                5, & ! 5
                                                6, & ! 6
                                                1, & ! 7
                                                2, & ! 8
                                                1, & ! 9
                                                2, & ! 10
                                                1, & ! 11
                                                2, & ! 12
                                                1, & ! 13
                                                2, & ! 14
                                                3, & ! 15
                                                4, & ! 16
                                                3, & ! 17
                                                4, & ! 18
                                                1, & ! 19
                                                2, & ! 20
                                                1, & ! 21
                                                2, & ! 22
                                                1, & ! 23
                                                2, & ! 24
                                                1, & ! 25
                                                2  & ! 26
                                                ]  !< Mapping fec1-26 to fec1-6 for boundaries.
! FEC_1_6_ARRAY([1,7,9,11,13,19,21,23,25])  = 1
! FEC_1_6_ARRAY([2,8,10,12,14,20,22,24,26]) = 2
! FEC_1_6_ARRAY([3,15,17])                  = 3
! FEC_1_6_ARRAY([4,16,18])                  = 4
! FEC_1_6_ARRAY([5])                        = 5
! FEC_1_6_ARRAY([6])                        = 6

integer(I4P), parameter :: FEC_TO_DELTA(3, 26) = reshape([-1,  0,  0, &! face 1
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

integer(I4P), parameter :: DELTA_TO_FEC(-1:1,-1:1,-1:1) = reshape([19, & ! -1, -1, -1   1
                                                                   15, & !  0, -1, -1   2
                                                                   20, & !  1, -1, -1   3
                                                                   11, & ! -1,  0, -1   4
                                                                   5 , & !  0,  0, -1   5
                                                                   12, & !  1,  0, -1   6
                                                                   21, & ! -1,  1, -1   7
                                                                   16, & !  0,  1, -1   8
                                                                   22, & !  1,  1, -1   9
                                                                    7, & ! -1, -1,  0   10
                                                                    3, & !  0, -1,  0   11
                                                                    8, & !  1, -1,  0   12
                                                                    1, & ! -1,  0,  0   13
                                                                   -1, & !  0,  0,  0
                                                                    2, & !  1,  0,  0   14
                                                                    9, & ! -1,  1,  0   15
                                                                    4, & !  0,  1,  0   16
                                                                   10, & !  1,  1,  0   17
                                                                   23, & ! -1, -1,  1   18
                                                                   17, & !  0, -1,  1   19
                                                                   24, & !  1, -1,  1   20
                                                                   13, & ! -1,  0,  1   21
                                                                    6, & !  0,  0,  1   22
                                                                   14, & !  1,  0,  1   23
                                                                   25, & ! -1,  1,  1   24
                                                                   18, & !  0,  1,  1   25
                                                                   26] & !  1,  1,  1   26
                                                                   , [3,3,3]) !< Delta to fec map.
endmodule adam_parameters
