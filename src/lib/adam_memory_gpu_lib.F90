!< ADAM, GPU-memory handling library.
module adam_memory_gpu_lib
!< ADAM, GPU-memory handling library.

use PENF
use CUDAFOR

implicit none
save
private
! public :: alloc_var_gpu
public :: assign_allocatable_gpu
public :: save_memory_gpu_status

interface alloc_var_gpu
!< Allocate GPU variable with memory checking.
module procedure alloc_var_gpu_R8P_1D, &
                 alloc_var_gpu_R8P_2D, &
                 alloc_var_gpu_R8P_5D, &
                 alloc_var_gpu_I4P_1D, &
                 alloc_var_gpu_I8P_1D, &
                 alloc_var_gpu_I8P_2D, &
                 alloc_var_gpu_I8P_3D
endinterface alloc_var_gpu

interface assign_allocatable_gpu
!< Assign GPU variable with memory checking.
module procedure  assign_allocatable_gpu_R8P_1D, &
                  assign_allocatable_gpu_R8P_2D, &
                  assign_allocatable_gpu_I4P_1D, &
                  assign_allocatable_gpu_I8P_2D, &
                  assign_allocatable_gpu_I8P_3D
endinterface assign_allocatable_gpu

interface transpose_a
   module procedure transpose_a_R8P_2D !< Transpose array (kind R8P, rank 2).
endinterface transpose_a

contains
   subroutine alloc_var_gpu_R8P_1D(var, ulb, msg, verbose)
   !< Allocate GPU variable with memory checking (kind R8P, rank 1).
   real(R8P), allocatable, intent(inout), device :: var(:)              !< Varibale to be allocate on GPU.
   integer(I4P),           intent(in)            :: ulb(2)              !< Upper/lower bounds of variable.
   character(*),           intent(in), optional  :: msg                 !< Message to be printed in verbose mode.
   logical,                intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                     :: msg_                !< Message to be printed in verbose mode, local var.
   logical                                       :: verbose_            !< Flag to activate verbose mode, local var.
   integer(cuda_count_kind)                      :: mem_free, mem_total !< Device memory.
   integer(I4P)                                  :: error               !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (allocated(var)) deallocate(var)
   if (verbose_) then
      error = cudaMemGetInfo(mem_free, mem_total)
      print '(A)', msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   allocate(var(ulb(1):ulb(2)))
   if (verbose_) then
      error = cudaMemGetInfo(mem_free, mem_total)
      print '(A)', msg_//'free/total memory AFTER  allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   endsubroutine alloc_var_gpu_R8P_1D

   subroutine alloc_var_gpu_R8P_2D(var, ulb, msg, verbose)
   !< Allocate GPU variable with memory checking (kind R8P, rank 2).
   real(R8P), allocatable, intent(inout), device :: var(:,:)            !< Varibale to be allocate on GPU.
   integer(I4P),           intent(in)            :: ulb(2,2)            !< Upper/lower bounds of variable.
   character(*),           intent(in), optional  :: msg                 !< Message to be printed in verbose mode.
   logical,                intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                     :: msg_                !< Message to be printed in verbose mode, local var.
   logical                                       :: verbose_            !< Flag to activate verbose mode, local var.
   integer(cuda_count_kind)                      :: mem_free, mem_total !< Device memory.
   integer(I4P)                                  :: error               !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (allocated(var)) deallocate(var)
   if (verbose_) then
      error = cudaMemGetInfo(mem_free, mem_total)
      print '(A)', msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   allocate(var(ulb(1,1):ulb(2,1), ulb(1,2):ulb(2,2)))
   if (verbose_) then
      error = cudaMemGetInfo(mem_free, mem_total)
      print '(A)', msg_//'free/total memory AFTER  allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   endsubroutine alloc_var_gpu_R8P_2D

   subroutine alloc_var_gpu_R8P_5D(var, ulb, msg, verbose)
   !< Allocate GPU variable with memory checking (kind R8P, rank 5).
   real(R8P), allocatable, intent(inout), device :: var(:,:,:,:,:)      !< Varibale to be allocate on GPU.
   integer(I4P),           intent(in)            :: ulb(2,5)            !< Upper/lower bounds of variable.
   character(*),           intent(in), optional  :: msg                 !< Message to be printed in verbose mode.
   logical,                intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                     :: msg_                !< Message to be printed in verbose mode, local var.
   logical                                       :: verbose_            !< Flag to activate verbose mode, local var.
   integer(cuda_count_kind)                      :: mem_free, mem_total !< Device memory.
   integer(I4P)                                  :: error               !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (allocated(var)) deallocate(var)
   if (verbose_) then
      error = cudaMemGetInfo(mem_free, mem_total)
      print '(A)', msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   allocate(var(ulb(1,1):ulb(2,1), ulb(1,2):ulb(2,2), ulb(1,3):ulb(2,3), ulb(1,4):ulb(2,4), ulb(1,5):ulb(2,5)))
   if (verbose_) then
      error = cudaMemGetInfo(mem_free, mem_total)
      print '(A)', msg_//'free/total memory AFTER  allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   endsubroutine alloc_var_gpu_R8P_5D

   subroutine alloc_var_gpu_I4P_1D(var, ulb, msg, verbose)
   !< Allocate GPU variable with memory checking (kind I4P, rank 1).
   integer(I4P), allocatable, intent(inout), device :: var(:)              !< Varibale to be allocate on GPU.
   integer(I4P),              intent(in)            :: ulb(2)              !< Upper/lower bounds of variable.
   character(*),              intent(in), optional  :: msg                 !< Message to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                        :: msg_                !< Message to be printed in verbose mode, local var.
   logical                                          :: verbose_            !< Flag to activate verbose mode, local var.
   integer(cuda_count_kind)                         :: mem_free, mem_total !< Device memory.
   integer(I4P)                                     :: error               !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (allocated(var)) deallocate(var)
   if (verbose_) then
      error = cudaMemGetInfo(mem_free, mem_total)
      print '(A)', msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   allocate(var(ulb(1):ulb(2)))
   if (verbose_) then
      error = cudaMemGetInfo(mem_free, mem_total)
      print '(A)', msg_//'free/total memory AFTER  allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   endsubroutine alloc_var_gpu_I4P_1D

   subroutine alloc_var_gpu_I8P_1D(var, ulb, msg, verbose)
   !< Allocate GPU variable with memory checking (kind I8P, rank 1).
   integer(I8P), allocatable, intent(inout), device :: var(:)              !< Varibale to be allocate on GPU.
   integer(I4P),              intent(in)            :: ulb(2)              !< Upper/lower bounds of variable.
   character(*),              intent(in), optional  :: msg                 !< Message to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                        :: msg_                !< Message to be printed in verbose mode, local var.
   logical                                          :: verbose_            !< Flag to activate verbose mode, local var.
   integer(cuda_count_kind)                         :: mem_free, mem_total !< Device memory.
   integer(I4P)                                     :: error               !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (allocated(var)) deallocate(var)
   if (verbose_) then
      error = cudaMemGetInfo(mem_free, mem_total)
      print '(A)', msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   allocate(var(ulb(1):ulb(2)))
   if (verbose_) then
      error = cudaMemGetInfo(mem_free, mem_total)
      print '(A)', msg_//'free/total memory AFTER  allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   endsubroutine alloc_var_gpu_I8P_1D

   subroutine alloc_var_gpu_I8P_2D(var, ulb, msg, verbose)
   !< Allocate GPU variable with memory checking (kind I8P, rank 2).
   integer(I8P), allocatable, intent(inout), device :: var(:,:)            !< Varibale to be allocate on GPU.
   integer(I4P),              intent(in)            :: ulb(2,2)            !< Upper/lower bounds of variable.
   character(*),              intent(in), optional  :: msg                 !< Message to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                        :: msg_                !< Message to be printed in verbose mode, local var.
   logical                                          :: verbose_            !< Flag to activate verbose mode, local var.
   integer(cuda_count_kind)                         :: mem_free, mem_total !< Device memory.
   integer(I4P)                                     :: error               !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (allocated(var)) deallocate(var)
   if (verbose_) then
      error = cudaMemGetInfo(mem_free, mem_total)
      print '(A)', msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   allocate(var(ulb(1,1):ulb(2,1), ulb(1,2):ulb(2,2)))
   if (verbose_) then
      error = cudaMemGetInfo(mem_free, mem_total)
      print '(A)', msg_//'free/total memory AFTER  allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   endsubroutine alloc_var_gpu_I8P_2D

   subroutine alloc_var_gpu_I8P_3D(var, ulb, msg, verbose)
   !< Allocate GPU variable with memory checking (kind I8P, rank 3).
   integer(I8P), allocatable, intent(inout), device :: var(:,:,:)          !< Varibale to be allocate on GPU.
   integer(I4P),              intent(in)            :: ulb(2,3)            !< Upper/lower bounds of variable.
   character(*),              intent(in), optional  :: msg                 !< Message to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                        :: msg_                !< Message to be printed in verbose mode, local var.
   logical                                          :: verbose_            !< Flag to activate verbose mode, local var.
   integer(cuda_count_kind)                         :: mem_free, mem_total !< Device memory.
   integer(I4P)                                     :: error               !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (allocated(var)) deallocate(var)
   if (verbose_) then
      error = cudaMemGetInfo(mem_free, mem_total)
      print '(A)', msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   allocate(var(ulb(1,1):ulb(2,1), ulb(1,2):ulb(2,2), ulb(1,3):ulb(2,3)))
   if (verbose_) then
      error = cudaMemGetInfo(mem_free, mem_total)
      print '(A)', msg_//'free/total memory AFTER  allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   endif
   endsubroutine alloc_var_gpu_I8P_3D

   subroutine assign_allocatable_gpu_R8P_1D(lhs, rhs, msg, verbose)
   !< Assign GPU variable with memory checking (kind R8P, rank 1).
   !< Variable is returned not allocated if right hand side is not allocated.
   real(R8P), allocatable, intent(inout), device :: lhs(:)  !< Left hand side of assignement.
   real(R8P), allocatable, intent(in)            :: rhs(:)  !< Right hand side of assignement.
   character(*),           intent(in), optional  :: msg     !< Message to be printed in verbose mode.
   logical,                intent(in), optional  :: verbose !< Flag to activate verbose mode.

   if (allocated(lhs)) deallocate(lhs)
   if (allocated(rhs)) then
      if (size(rhs, dim=1)>0) then
         call alloc_var_gpu(var=lhs, ulb=[lbound(rhs,dim=1),ubound(rhs,dim=1)], msg=msg, verbose=verbose)
         lhs = rhs
      endif
   endif
   endsubroutine assign_allocatable_gpu_R8P_1D

   subroutine assign_allocatable_gpu_R8P_2D(lhs, rhs, transposed, msg, verbose)
   !< Assign GPU variable with memory checking (kind R8P, rank 2).
   !< Variable is returned not allocated if right hand side is not allocated.
   real(R8P), allocatable, intent(inout), device :: lhs(:,:)    !< Left hand side of assignement.
   real(R8P), allocatable, intent(in)            :: rhs(:,:)    !< Right hand side of assignement.
   logical,                intent(in), optional  :: transposed  !< Assign trasposed rhs.
   character(*),           intent(in), optional  :: msg         !< Message to be printed in verbose mode.
   logical,                intent(in), optional  :: verbose     !< Flag to activate verbose mode.
   logical                                       :: transposed_ !< Assign trasposed rhs, local var.
   integer(I4P)                                  :: ulb(2,2)    !< Upper/lower bounds of variable.
   real(R8P), allocatable                        :: rhst(:,:)   !< Right hand side transposed.

   if (allocated(lhs)) deallocate(lhs)
   if (allocated(rhs)) then
      if (size(rhs, dim=1)*size(rhs, dim=2)>0) then
         transposed_ = .false. ; if (present(transposed)) transposed_ = transposed
         ulb=reshape([lbound(rhs, dim=1),ubound(rhs, dim=1),lbound(rhs, dim=2),ubound(rhs, dim=2)], [2,2])
         if (transposed_) then
            allocate(rhst(ulb(1,2):ulb(2,2),ulb(1,1):ulb(2,1)))
            call alloc_var_gpu(var=lhs, ulb=reshape([ulb(1,2),ulb(2,2),ulb(1,1),ulb(2,1)],[2,2]), msg=msg, verbose=verbose)
            call transpose_a(ii=ulb(:,1), jj=ulb(:,2), a=rhs, t=rhst)
            lhs = rhst
            deallocate(rhst)
         else
            call alloc_var_gpu(var=lhs, ulb=ulb, msg=msg, verbose=verbose)
            lhs = rhs
         endif
      endif
   endif
   endsubroutine assign_allocatable_gpu_R8P_2D

   subroutine assign_allocatable_gpu_I4P_1D(lhs, rhs, msg, verbose)
   !< Assign GPU variable with memory checking (kind I4P, rank 1).
   !< Variable is returned not allocated if right hand side is not allocated.
   integer(I4P), allocatable, intent(inout), device :: lhs(:)  !< Varibale to be allocate on GPU.
   integer(I4P), allocatable, intent(in)            :: rhs(:)  !< Right hand side of assignement.
   character(*),              intent(in), optional  :: msg     !< Message to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose !< Flag to activate verbose mode.

   if (allocated(lhs)) deallocate(lhs)
   if (allocated(rhs)) then
      if (size(rhs, dim=1)>0) then
         call alloc_var_gpu(var=lhs, ulb=[lbound(rhs,dim=1),ubound(rhs,dim=1)], msg=msg, verbose=verbose)
         lhs = rhs
      endif
   endif
   endsubroutine assign_allocatable_gpu_I4P_1D

   subroutine assign_allocatable_gpu_I8P_2D(lhs, rhs, msg, verbose)
   !< Assign GPU variable with memory checking (kind I8P, rank 2).
   !< Variable is returned not allocated if right hand side is not allocated.
   integer(I8P), allocatable, intent(inout), device :: lhs(:,:) !< Left hand side of assignement.
   integer(I8P), allocatable, intent(in)            :: rhs(:,:) !< Right hand side of assignement.
   character(*),              intent(in), optional  :: msg      !< Message to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose  !< Flag to activate verbose mode.

   if (allocated(lhs)) deallocate(lhs)
   if (allocated(rhs)) then
      if (size(rhs, dim=1)*size(rhs, dim=2)>0) then
         call alloc_var_gpu(var=lhs,                                                                                     &
                            ulb=reshape([lbound(rhs,dim=1),ubound(rhs,dim=1),lbound(rhs,dim=2),ubound(rhs,dim=2)],[2,2]),&
                            msg=msg, verbose=verbose)
         lhs = rhs
      endif
   endif
   endsubroutine assign_allocatable_gpu_I8P_2D

   subroutine assign_allocatable_gpu_I8P_3D(lhs, rhs, msg, verbose)
   !< Assign GPU variable with memory checking (kind I8P, rank 3).
   !< Variable is returned not allocated if right hand side is not allocated.
   integer(I8P), allocatable, intent(inout), device :: lhs(:,:,:) !< Left hand side of assignement.
   integer(I8P), allocatable, intent(in)            :: rhs(:,:,:) !< Right hand side of assignement.
   character(*),              intent(in), optional  :: msg        !< Message to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose    !< Flag to activate verbose mode.

   if (allocated(lhs)) deallocate(lhs)
   if (allocated(rhs)) then
      if (size(rhs, dim=1)*size(rhs, dim=2)*size(rhs, dim=3)>0) then
         call alloc_var_gpu(var=lhs,                                                  &
                            ulb=reshape([lbound(rhs,dim=1),ubound(rhs,dim=1),         &
                                         lbound(rhs,dim=2),ubound(rhs,dim=2),         &
                                         lbound(rhs,dim=3),ubound(rhs,dim=3)],[2,3]), &
                            msg=msg, verbose=verbose)
         lhs = rhs
      endif
   endif
   endsubroutine assign_allocatable_gpu_I8P_3D

   subroutine save_memory_gpu_status(file_name, tag)
   !< Save the current CPU-memory status into a file.
   !< File is accessed in append position.
   character(*), intent(in)           :: file_name           !< File name.
   character(*), intent(in), optional :: tag                 !< Tag of current status.
   character(:), allocatable          :: tag_                !< Tag of current status, local var.
   integer(cuda_count_kind)           :: mem_free, mem_total !< Device memory.
   integer(I4P)                       :: file_unit           !< File unit.
   integer(I4P)                       :: error               !< Error traping flag.

   tag_ = '' ; if (present(tag)) tag_ = trim(tag)
   error = cudaMemGetInfo(mem_free, mem_total)
   open(newunit=file_unit, file=trim(file_name), position="append")
   write(file_unit,*) tag_, mem_free, mem_total
   close(file_unit)
   endsubroutine save_memory_gpu_status

   subroutine transpose_a_R8P_2D(ii, jj, a, t)
   !< Transpose array (kind R8P, rank 2).
   integer(I4P), intent(in)  :: ii(2), jj(2)               !< Array bounds.
   real(R8P),    intent(in)  :: a(ii(1):ii(2),jj(1):jj(2)) !< Input array.
   real(R8P),    intent(out) :: t(jj(1):jj(2),ii(1):ii(2)) !< Transposed array.
   integer(I4P)              :: i, j                       !< Counter.

   do j=jj(1), jj(2)
      do i=ii(1), ii(2)
         t(j,i) = a(i,j)
      enddo
   enddo
   endsubroutine transpose_a_R8P_2D
endmodule adam_memory_gpu_lib
