!< ADAM, CPU-memory handling library.
module adam_memory_cpu_lib
!< ADAM, CPU-memory handling library.

use adam_mpi_lib
use PENF
use memorysaver, only : getmemory
use, intrinsic :: iso_c_binding

implicit none
save
private
public :: alloc_var_cpu
public :: assign_allocatable_cpu

interface alloc_var_cpu
!< Allocate CPU variable with memory checking.
module procedure alloc_var_cpu_R8P_1D, &
                 alloc_var_cpu_R8P_2D, &
                 alloc_var_cpu_R8P_5D, &
                 alloc_var_cpu_I4P_1D, &
                 alloc_var_cpu_I8P_1D, &
                 alloc_var_cpu_I8P_2D, &
                 alloc_var_cpu_I8P_3D, &
                 alloc_var_cpu_I4P_2D
endinterface alloc_var_cpu

interface assign_allocatable_cpu
!< Assign CPU variable with memory checking.
module procedure  assign_allocatable_cpu_R8P_2D, &
                  assign_allocatable_cpu_I8P_2D, &
                  assign_allocatable_cpu_I4P_1D
endinterface assign_allocatable_cpu

contains
   subroutine alloc_var_cpu_R8P_1D(var, ulb, msg, verbose)
   !< Allocate CPU variable with memory checking (kind R8P, rank 1).
   real(R8P), allocatable, intent(inout)         :: var(:)              !< Varibale to be allocate on CPU.
   integer(I4P),           intent(in)            :: ulb(2)              !< Upper/lower bounds of variable.
   character(*),           intent(in), optional  :: msg                 !< Message to be printed in verbose mode.
   logical,                intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                     :: msg_                !< Message to be printed in verbose mode, local var.
   logical                                       :: verbose_            !< Flag to activate verbose mode, local var.
   integer(C_LONG)                               :: mem_free, mem_total !< CPU memory.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (allocated(var)) deallocate(var)
   if (verbose_) then
      call getmemory(mem_total, mem_free)
      print '(A)', myrankstr//msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   allocate(var(ulb(1):ulb(2)))
   if (verbose_) then
      call getmemory(mem_total, mem_free)
      print '(A)', myrankstr//msg_//'free/total memory AFTER  allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   endsubroutine alloc_var_cpu_R8P_1D

   subroutine alloc_var_cpu_R8P_2D(var, ulb, msg, verbose)
   !< Allocate CPU variable with memory checking (kind R8P, rank 2).
   real(R8P), allocatable, intent(inout)         :: var(:,:)            !< Varibale to be allocate on CPU.
   integer(I4P),           intent(in)            :: ulb(2,2)            !< Upper/lower bounds of variable.
   character(*),           intent(in), optional  :: msg                 !< Message to be printed in verbose mode.
   logical,                intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                     :: msg_                !< Message to be printed in verbose mode, local var.
   logical                                       :: verbose_            !< Flag to activate verbose mode, local var.
   integer(C_LONG)                               :: mem_free, mem_total !< CPU memory.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (allocated(var)) deallocate(var)
   if (verbose_) then
      call getmemory(mem_total, mem_free)
      print '(A)', myrankstr//msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   allocate(var(ulb(1,1):ulb(2,1), ulb(1,2):ulb(2,2)))
   if (verbose_) then
      call getmemory(mem_total, mem_free)
      print '(A)', myrankstr//msg_//'free/total memory AFTER  allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   endsubroutine alloc_var_cpu_R8P_2D

   subroutine alloc_var_cpu_R8P_5D(var, ulb, msg, verbose)
   !< Allocate CPU variable with memory checking (kind R8P, rank 5).
   real(R8P), allocatable, intent(inout)         :: var(:,:,:,:,:)      !< Varibale to be allocate on CPU.
   integer(I4P),           intent(in)            :: ulb(2,5)            !< Upper/lower bounds of variable.
   character(*),           intent(in), optional  :: msg                 !< Message to be printed in verbose mode.
   logical,                intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                     :: msg_                !< Message to be printed in verbose mode, local var.
   logical                                       :: verbose_            !< Flag to activate verbose mode, local var.
   integer(C_LONG)                               :: mem_free, mem_total !< CPU memory.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (allocated(var)) deallocate(var)
   if (verbose_) then
      call getmemory(mem_total, mem_free)
      print '(A)', myrankstr//msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   allocate(var(ulb(1,1):ulb(2,1), ulb(1,2):ulb(2,2), ulb(1,3):ulb(2,3), ulb(1,4):ulb(2,4), ulb(1,5):ulb(2,5)))
   if (verbose_) then
      call getmemory(mem_total, mem_free)
      print '(A)', myrankstr//msg_//'free/total memory AFTER  allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   endsubroutine alloc_var_cpu_R8P_5D

   subroutine alloc_var_cpu_I4P_1D(var, ulb, msg, verbose)
   !< Allocate CPU variable with memory checking (kind I4P, rank 1).
   integer(I4P), allocatable, intent(inout)         :: var(:)              !< Varibale to be allocate on CPU.
   integer(I4P),              intent(in)            :: ulb(2)              !< Upper/lower bounds of variable.
   character(*),              intent(in), optional  :: msg                 !< Message to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                        :: msg_                !< Message to be printed in verbose mode, local var.
   logical                                          :: verbose_            !< Flag to activate verbose mode, local var.
   integer(C_LONG)                                  :: mem_free, mem_total !< CPU memory.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (allocated(var)) deallocate(var)
   if (verbose_) then
      call getmemory(mem_total, mem_free)
      print '(A)', myrankstr//msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   allocate(var(ulb(1):ulb(2)))
   if (verbose_) then
      call getmemory(mem_total, mem_free)
      print '(A)', myrankstr//msg_//'free/total memory AFTER  allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   endsubroutine alloc_var_cpu_I4P_1D

   subroutine alloc_var_cpu_I8P_1D(var, ulb, msg, verbose)
   !< Allocate CPU variable with memory checking (kind I8P, rank 1).
   integer(I8P), allocatable, intent(inout)         :: var(:)              !< Varibale to be allocate on CPU.
   integer(I4P),              intent(in)            :: ulb(2)              !< Upper/lower bounds of variable.
   character(*),              intent(in), optional  :: msg                 !< Message to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                        :: msg_                !< Message to be printed in verbose mode, local var.
   logical                                          :: verbose_            !< Flag to activate verbose mode, local var.
   integer(C_LONG)                                  :: mem_free, mem_total !< CPU memory.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (allocated(var)) deallocate(var)
   if (verbose_) then
      call getmemory(mem_total, mem_free)
      print '(A)', myrankstr//msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   allocate(var(ulb(1):ulb(2)))
   if (verbose_) then
      call getmemory(mem_total, mem_free)
      print '(A)', myrankstr//msg_//'free/total memory AFTER  allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   endsubroutine alloc_var_cpu_I8P_1D

   subroutine alloc_var_cpu_I8P_2D(var, ulb, msg, verbose)
   !< Allocate CPU variable with memory checking (kind I8P, rank 2).
   integer(I8P), allocatable, intent(inout)         :: var(:,:)            !< Varibale to be allocate on CPU.
   integer(I4P),              intent(in)            :: ulb(2,2)            !< Upper/lower bounds of variable.
   character(*),              intent(in), optional  :: msg                 !< Message to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                        :: msg_                !< Message to be printed in verbose mode, local var.
   logical                                          :: verbose_            !< Flag to activate verbose mode, local var.
   integer(C_LONG)                                  :: mem_free, mem_total !< CPU memory.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (allocated(var)) deallocate(var)
   if (verbose_) then
      call getmemory(mem_total, mem_free)
      print '(A)', myrankstr//msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   allocate(var(ulb(1,1):ulb(2,1), ulb(1,2):ulb(2,2)))
   if (verbose_) then
      call getmemory(mem_total, mem_free)
      print '(A)', myrankstr//msg_//'free/total memory AFTER  allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   endsubroutine alloc_var_cpu_I8P_2D

   subroutine alloc_var_cpu_I8P_3D(var, ulb, msg, verbose)
   !< Allocate CPU variable with memory checking (kind I8P, rank 3).
   integer(I8P), allocatable, intent(inout)         :: var(:,:,:)          !< Varibale to be allocate on CPU.
   integer(I4P),              intent(in)            :: ulb(2,3)            !< Upper/lower bounds of variable.
   character(*),              intent(in), optional  :: msg                 !< Message to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                        :: msg_                !< Message to be printed in verbose mode, local var.
   logical                                          :: verbose_            !< Flag to activate verbose mode, local var.
   integer(C_LONG)                                  :: mem_free, mem_total !< CPU memory.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (allocated(var)) deallocate(var)
   if (verbose_) then
      call getmemory(mem_total, mem_free)
      print '(A)', myrankstr//msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   allocate(var(ulb(1,1):ulb(2,1), ulb(1,2):ulb(2,2), ulb(1,3):ulb(2,3)))
   if (verbose_) then
      call getmemory(mem_total, mem_free)
      print '(A)', myrankstr//msg_//'free/total memory AFTER  allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   endsubroutine alloc_var_cpu_I8P_3D

   subroutine alloc_var_cpu_I4P_2D(var, ulb, msg, verbose)
   !< Allocate CPU variable with memory checking (kind I4P, rank 2).
   integer(I4P), allocatable, intent(inout)         :: var(:,:)            !< Varibale to be allocate on CPU.
   integer(I4P),              intent(in)            :: ulb(2,2)            !< Upper/lower bounds of variable.
   character(*),              intent(in), optional  :: msg                 !< Message to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                        :: msg_                !< Message to be printed in verbose mode, local var.
   logical                                          :: verbose_            !< Flag to activate verbose mode, local var.
   integer(C_LONG)                                  :: mem_free, mem_total !< CPU memory.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (allocated(var)) deallocate(var)
   if (verbose_) then
      call getmemory(mem_total, mem_free)
      print '(A)', myrankstr//msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   allocate(var(ulb(1,1):ulb(2,1), ulb(1,2):ulb(2,2)))
   if (verbose_) then
      call getmemory(mem_total, mem_free)
      print '(A)', myrankstr//msg_//'free/total memory AFTER  allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   endsubroutine alloc_var_cpu_I4P_2D

   subroutine assign_allocatable_cpu_R8P_2D(lhs, rhs, msg, verbose)
   !< Assign CPU variable with memory checking (kind R8P, rank 2).
   !< Variable is returned not allocated if right hand side is not allocated.
   real(R8P), allocatable, intent(inout)         :: lhs(:,:) !< Left hand side of assignement.
   real(R8P), allocatable, intent(in)            :: rhs(:,:) !< Right hand side of assignement.
   character(*),           intent(in), optional  :: msg      !< Message to be printed in verbose mode.
   logical,                intent(in), optional  :: verbose  !< Flag to activate verbose mode.

   if (allocated(lhs)) deallocate(lhs)
   if (allocated(rhs)) then
      if (size(rhs, dim=1)*size(rhs, dim=2)>0) then
         call alloc_var_cpu(var=lhs,                                                  &
                                 ulb=reshape([lbound(rhs,dim=1),ubound(rhs,dim=1),         &
                                              lbound(rhs,dim=2),ubound(rhs,dim=2)],[2,2]), &
                                 msg=msg, verbose=verbose)
         lhs = rhs
      endif
   endif
   endsubroutine assign_allocatable_cpu_R8P_2D

   subroutine assign_allocatable_cpu_I8P_2D(lhs, rhs, msg, verbose)
   !< Assign CPU variable with memory checking (kind I8P, rank 2).
   !< Variable is returned not allocated if right hand side is not allocated.
   integer(I8P), allocatable, intent(inout)         :: lhs(:,:) !< Left hand side of assignement.
   integer(I8P), allocatable, intent(in)            :: rhs(:,:) !< Right hand side of assignement.
   character(*),              intent(in), optional  :: msg      !< Message to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose  !< Flag to activate verbose mode.

   if (allocated(lhs)) deallocate(lhs)
   if (allocated(rhs)) then
      if (size(rhs, dim=1)*size(rhs, dim=2)>0) then
         call alloc_var_cpu(var=lhs,                                                  &
                            ulb=reshape([lbound(rhs,dim=1),ubound(rhs,dim=1),         &
                                         lbound(rhs,dim=2),ubound(rhs,dim=2)],[2,2]), &
                            msg=msg, verbose=verbose)
         lhs = rhs
      endif
   endif
   endsubroutine assign_allocatable_cpu_I8P_2D

   subroutine assign_allocatable_cpu_I4P_1D(lhs, rhs, msg, verbose)
   !< Assign CPU variable with memory checking (kind I4P, rank 1).
   !< Variable is returned not allocated if right hand side is not allocated.
   integer(I4P), allocatable, intent(inout)         :: lhs(:)  !< Left hand side of assignement.
   integer(I4P), allocatable, intent(in)            :: rhs(:)  !< Right hand side of assignement.
   character(*),              intent(in), optional  :: msg     !< Message to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose !< Flag to activate verbose mode.

   if (allocated(lhs)) deallocate(lhs)
   if (allocated(rhs)) then
      if (size(rhs, dim=1)>0) then
         call alloc_var_cpu(var=lhs,                                   &
                            ulb=[lbound(rhs,dim=1),ubound(rhs,dim=1)], &
                            msg=msg, verbose=verbose)
         lhs = rhs
      endif
   endif
   endsubroutine assign_allocatable_cpu_I4P_1D
endmodule adam_memory_cpu_lib
