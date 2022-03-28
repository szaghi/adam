!< ADAM, general parameters.
module adam_parameters
!< ADAM, general parameters.

use PENF

implicit none
private
public :: BC_PERIODIC
public :: TO_BE_REFINED,   &
          TO_BE_DEREFINED, &
          TO_NOT_TOUCH
public :: fec_to_delta
public :: delta_to_fec
public :: assign_allocatable

integer(I4P), parameter :: BC_PERIODIC = -1_I4P !< Flag (reserved) for periodical boundary conditions.

integer(I4P), parameter :: TO_BE_REFINED=1_I4P    !< Flag for node/block to be refined.
integer(I4P), parameter :: TO_BE_DEREFINED=-1_I4P !< Flag for node/block to be derefined.
integer(I4P), parameter :: TO_NOT_TOUCH=0_I4P     !< Flag for node/block to be untouched.

integer(I4P), parameter :: fec_to_delta(3, 26) = reshape([-1,  0,  0, &! face 1
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

integer(I4P), parameter :: delta_to_fec(-1:1,-1:1,-1:1) = reshape([19, & ! -1, -1, -1   1
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

interface assign_allocatable
   !< Safe assign allocatable arrays, generic interface.
   module procedure assign_allocatable_char   !< Safe assign allocatable character.
   module procedure assign_allocatable_I8P_1D !< Safe assign allocatable arrays, I8P 1D type.
   module procedure assign_allocatable_I8P_2D !< Safe assign allocatable arrays, I8P 2D type.
   module procedure assign_allocatable_I4P_1D !< Safe assign allocatable arrays, I4P 1D type.
   module procedure assign_allocatable_I4P_2D !< Safe assign allocatable arrays, I4P 2D type.
   module procedure assign_allocatable_R8P_1D !< Safe assign allocatable arrays, R8P 1D type.
   module procedure assign_allocatable_R8P_2D !< Safe assign allocatable arrays, R8P 2D type.
   module procedure assign_allocatable_R8P_3D !< Safe assign allocatable arrays, R8P 3D type.
   module procedure assign_allocatable_R8P_4D !< Safe assign allocatable arrays, R8P 4D type.
   module procedure assign_allocatable_R8P_5D !< Safe assign allocatable arrays, R8P 5D type.
   module procedure assign_allocatable_R8P_6D !< Safe assign allocatable arrays, R8P 6D type.
endinterface assign_allocatable

contains
   ! private procedures
   pure subroutine assign_allocatable_char(lhs, rhs)
   !< Safe assign allocatable arrays, I8P 1D type.
   character(:), allocatable, intent(inout) :: lhs !< Left hand side.
   character(:), allocatable, intent(in)    :: rhs !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_char

   pure subroutine assign_allocatable_I8P_1D(lhs, rhs)
   !< Safe assign allocatable arrays, I8P 1D type.
   integer(I8P), allocatable, intent(inout) :: lhs(:) !< Left hand side.
   integer(I8P), allocatable, intent(in)    :: rhs(:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_I8P_1D

   pure subroutine assign_allocatable_I8P_2D(lhs, rhs)
   !< Safe assign allocatable arrays, I8P 2D type.
   integer(I8P), allocatable, intent(inout) :: lhs(:,:) !< Left hand side.
   integer(I8P), allocatable, intent(in)    :: rhs(:,:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_I8P_2D

   pure subroutine assign_allocatable_I4P_1D(lhs, rhs)
   !< Safe assign allocatable arrays, I4P 1D type.
   integer(I4P), allocatable, intent(inout) :: lhs(:) !< Left hand side.
   integer(I4P), allocatable, intent(in)    :: rhs(:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
         ! print*, 'cazzo lhs ', lhs, size(lhs)
         ! print*, 'cazzo rhs ', rhs, size(rhs)
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_I4P_1D

   pure subroutine assign_allocatable_I4P_2D(lhs, rhs)
   !< Safe assign allocatable arrays, I4P 2D type.
   integer(I4P), allocatable, intent(inout) :: lhs(:,:) !< Left hand side.
   integer(I4P), allocatable, intent(in)    :: rhs(:,:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_I4P_2D

   pure subroutine assign_allocatable_R8P_1D(lhs, rhs)
   !< Safe assign allocatable arrays, R8P 1D type.
   real(R8P), allocatable, intent(inout) :: lhs(:) !< Left hand side.
   real(R8P), allocatable, intent(in)    :: rhs(:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_R8P_1D

   pure subroutine assign_allocatable_R8P_2D(lhs, rhs)
   !< Safe assign allocatable arrays, R8P 2D type.
   real(R8P), allocatable, intent(inout) :: lhs(:,:) !< Left hand side.
   real(R8P), allocatable, intent(in)    :: rhs(:,:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_R8P_2D

   pure subroutine assign_allocatable_R8P_3D(lhs, rhs)
   !< Safe assign allocatable arrays, R8P 3D type.
   real(R8P), allocatable, intent(inout) :: lhs(:,:,:) !< Left hand side.
   real(R8P), allocatable, intent(in)    :: rhs(:,:,:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_R8P_3D

   pure subroutine assign_allocatable_R8P_4D(lhs, rhs)
   !< Safe assign allocatable arrays, R8P 4D type.
   real(R8P), allocatable, intent(inout) :: lhs(:,:,:,:) !< Left hand side.
   real(R8P), allocatable, intent(in)    :: rhs(:,:,:,:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_R8P_4D

   pure subroutine assign_allocatable_R8P_5D(lhs, rhs)
   !< Safe assign allocatable arrays, R8P 5D type.
   real(R8P), allocatable, intent(inout) :: lhs(:,:,:,:,:) !< Left hand side.
   real(R8P), allocatable, intent(in)    :: rhs(:,:,:,:,:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_R8P_5D

   pure subroutine assign_allocatable_R8P_6D(lhs, rhs)
   !< Safe assign allocatable arrays, R8P 6D type.
   real(R8P), allocatable, intent(inout) :: lhs(:,:,:,:,:,:) !< Left hand side.
   real(R8P), allocatable, intent(in)    :: rhs(:,:,:,:,:,:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_R8P_6D
endmodule adam_parameters
