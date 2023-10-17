!< ADAM, GMP memory handling library.
module adam_memory_gmp_lib
!< ADAM, GMP memory handling library.

use penf
use adam_gmp_utils
use omp_lib, only : omp_get_initial_device

implicit none
save
private
! public :: alloc_var_gpu
public :: assign_allocatable_gpu
!public :: save_memory_gpu_status

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
   subroutine alloc_var_gpu_R8P_1D(var, ulb, omp_dev, varname, verbose)
   !< Allocate GPU variable with memory checking (kind R8P, rank 1).
   real(R8P), pointer, intent(inout)         :: var(:)              !< Varibale to be allocate on GPU.
   integer(I4P),       intent(in)            :: ulb(2)              !< Upper/lower bounds of variable.
   integer(I4P),       intent(in)            :: omp_dev             !< OpenMP device ID.
   character(*),       intent(in), optional  :: varname             !< Variable name to be printed in verbose mode.
   logical,            intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                 :: varname_            !< Variable name to be printed in verbose mode, local var.
   logical                                   :: verbose_            !< Flag to activate verbose mode, local var.
   integer(I4P)                              :: error               !< Error traping flag.

   varname_ = ''      ; if (present(varname)) varname_ = varname
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   !if (verbose_) then
   !   error = cudaMemGetInfo(mem_free, mem_total)
   !   print '(A)', msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   !endif
   call omp_target_alloc_f(fptr_dev=var, ubounds=[ulb(2)], lbounds=[ulb(1)], omp_dev=omp_dev, ierr=error)
   if (error/=0) then
      print '(A)', 'Error allocating variable '//varname_
      error stop
   else
      if (verbose_) print '(A)', varname_//' : allocation OK!'
   endif
   endsubroutine alloc_var_gpu_R8P_1D

   subroutine alloc_var_gpu_R8P_2D(var, ulb, omp_dev, varname, verbose)
   !< Allocate GPU variable with memory checking (kind R8P, rank 2).
   real(R8P), pointer, intent(inout)         :: var(:,:)            !< Varibale to be allocate on GPU.
   integer(I4P),       intent(in)            :: ulb(2,2)            !< Upper/lower bounds of variable.
   integer(I4P),       intent(in)            :: omp_dev             !< OpenMP device ID.
   character(*),       intent(in), optional  :: varname             !< Variable name to be printed in verbose mode.
   logical,            intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                 :: varname_            !< Variable name to be printed in verbose mode, local var.
   logical                                   :: verbose_            !< Flag to activate verbose mode, local var.
   integer(I4P)                              :: error               !< Error traping flag.

   varname_ = ''      ; if (present(varname)) varname_ = varname
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   !if (verbose_) then
   !   error = cudaMemGetInfo(mem_free, mem_total)
   !   print '(A)', msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   !endif
   call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, ierr=error)
   if (error/=0) then
      print '(A)', 'Error allocating variable '//varname_
      error stop
   else
      if (verbose_) print '(A)', varname_//' : allocation OK!'
   endif
   endsubroutine alloc_var_gpu_R8P_2D

   subroutine alloc_var_gpu_R8P_5D(var, ulb, omp_dev, varname, verbose)
   !< Allocate GPU variable with memory checking (kind R8P, rank 5).
   real(R8P), pointer, intent(inout)         :: var(:,:,:,:,:)      !< Varibale to be allocate on GPU.
   integer(I4P),       intent(in)            :: ulb(2,5)            !< Upper/lower bounds of variable.
   integer(I4P),       intent(in)            :: omp_dev             !< OpenMP device ID.
   character(*),       intent(in), optional  :: varname             !< Variable name to be printed in verbose mode.
   logical,            intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                 :: varname_            !< Variable name to be printed in verbose mode, local var.
   logical                                   :: verbose_            !< Flag to activate verbose mode, local var.
   integer(I4P)                              :: error               !< Error traping flag.

   varname_ = ''      ; if (present(varname)) varname_ = varname
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   !if (verbose_) then
   !   error = cudaMemGetInfo(mem_free, mem_total)
   !   print '(A)', msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   !endif
   call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, ierr=error)
   if (error/=0) then
      print '(A)', 'Error allocating variable '//varname_
      error stop
   else
      if (verbose_) print '(A)', varname_//' : allocation OK!'
   endif
   endsubroutine alloc_var_gpu_R8P_5D

   subroutine alloc_var_gpu_I4P_1D(var, ulb, omp_dev, varname, verbose)
   !< Allocate GPU variable with memory checking (kind I4P, rank 1).
   integer(I4P), pointer, intent(inout)         :: var(:)              !< Varibale to be allocate on GPU.
   integer(I4P),          intent(in)            :: ulb(2)              !< Upper/lower bounds of variable.
   integer(I4P),          intent(in)            :: omp_dev             !< OpenMP device ID.
   character(*),          intent(in), optional  :: varname             !< Variable name to be printed in verbose mode.
   logical,               intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                    :: varname_            !< Variable name to be printed in verbose mode, local var.
   logical                                      :: verbose_            !< Flag to activate verbose mode, local var.
   integer(I4P)                                 :: error               !< Error traping flag.

   varname_ = ''      ; if (present(varname)) varname_ = varname
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   !if (verbose_) then
   !   error = cudaMemGetInfo(mem_free, mem_total)
   !   print '(A)', msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   !endif
   call omp_target_alloc_f(fptr_dev=var, ubounds=[ulb(2)], lbounds=[ulb(1)], omp_dev=omp_dev, ierr=error)
   if (error/=0) then
      print '(A)', 'Error allocating variable '//varname_
      error stop
   else
      if (verbose_) print '(A)', varname_//' : allocation OK!'
   endif
   endsubroutine alloc_var_gpu_I4P_1D

   subroutine alloc_var_gpu_I8P_1D(var, ulb, omp_dev, varname, verbose)
   !< Allocate GPU variable with memory checking (kind I8P, rank 1).
   integer(I8P), pointer, intent(inout)         :: var(:)              !< Varibale to be allocate on GPU.
   integer(I4P),          intent(in)            :: ulb(2)              !< Upper/lower bounds of variable.
   integer(I4P),          intent(in)            :: omp_dev             !< OpenMP device ID.
   character(*),          intent(in), optional  :: varname             !< Variable name to be printed in verbose mode.
   logical,               intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                    :: varname_            !< Variable name to be printed in verbose mode, local var.
   logical                                      :: verbose_            !< Flag to activate verbose mode, local var.
   integer(I4P)                                 :: error               !< Error traping flag.

   varname_ = ''      ; if (present(varname)) varname_ = varname
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   !if (verbose_) then
   !   error = cudaMemGetInfo(mem_free, mem_total)
   !   print '(A)', msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   !endif
   call omp_target_alloc_f(fptr_dev=var, ubounds=[ulb(2)], lbounds=[ulb(1)], omp_dev=omp_dev, ierr=error)
   if (error/=0) then
      print '(A)', 'Error allocating variable '//varname_
      error stop
   else
      if (verbose_) print '(A)', varname_//' : allocation OK!'
   endif
   endsubroutine alloc_var_gpu_I8P_1D

   subroutine alloc_var_gpu_I8P_2D(var, ulb, omp_dev, varname, verbose)
   !< Allocate GPU variable with memory checking (kind I8P, rank 2).
   integer(I8P), pointer, intent(inout)         :: var(:,:)            !< Varibale to be allocate on GPU.
   integer(I4P),          intent(in)            :: ulb(2,2)            !< Upper/lower bounds of variable.
   integer(I4P),          intent(in)            :: omp_dev             !< OpenMP device ID.
   character(*),          intent(in), optional  :: varname             !< Variable name to be printed in verbose mode.
   logical,               intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                    :: varname_            !< Variable name to be printed in verbose mode, local var.
   logical                                      :: verbose_            !< Flag to activate verbose mode, local var.
   integer(I4P)                                 :: error               !< Error traping flag.

   varname_ = ''      ; if (present(varname)) varname_ = varname
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   !if (verbose_) then
   !   error = cudaMemGetInfo(mem_free, mem_total)
   !   print '(A)', msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   !endif
   call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, ierr=error)
   if (error/=0) then
      print '(A)', 'Error allocating variable '//varname_
      error stop
   else
      if (verbose_) print '(A)', varname_//' : allocation OK!'
   endif
   endsubroutine alloc_var_gpu_I8P_2D

   subroutine alloc_var_gpu_I8P_3D(var, ulb, omp_dev, varname, verbose)
   !< Allocate GPU variable with memory checking (kind I8P, rank 2).
   integer(I8P), pointer, intent(inout)         :: var(:,:,:)          !< Varibale to be allocate on GPU.
   integer(I4P),          intent(in)            :: ulb(2,3)            !< Upper/lower bounds of variable.
   integer(I4P),          intent(in)            :: omp_dev             !< OpenMP device ID.
   character(*),          intent(in), optional  :: varname             !< Variable name to be printed in verbose mode.
   logical,               intent(in), optional  :: verbose             !< Flag to activate verbose mode.
   character(:), allocatable                    :: varname_            !< Variable name to be printed in verbose mode, local var.
   logical                                      :: verbose_            !< Flag to activate verbose mode, local var.
   integer(I4P)                                 :: error               !< Error traping flag.

   varname_ = ''      ; if (present(varname)) varname_ = varname
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   !if (verbose_) then
   !   error = cudaMemGetInfo(mem_free, mem_total)
   !   print '(A)', msg_//'free/total memory BEFORE allocate:'//trim(str([mem_free,mem_total]))//'[bytes]'
   !endif
   call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=ulb(1,:), omp_dev=omp_dev, ierr=error)
   if (error/=0) then
      print '(A)', 'Error allocating variable '//varname_
      error stop
   else
      if (verbose_) print '(A)', varname_//' : allocation OK!'
   endif
   endsubroutine alloc_var_gpu_I8P_3D

   subroutine assign_allocatable_gpu_R8P_1D(lhs, rhs, omp_dev, varname, verbose)
   !< Assign GPU variable with memory checking (kind R8P, rank 1).
   !< Variable is returned not allocated if right hand side is not allocated.
   real(R8P), pointer,     intent(inout)         :: lhs(:)    !< Left hand side of assignement.
   real(R8P), allocatable, intent(in), target    :: rhs(:)    !< Right hand side of assignement.
   integer(I4P),           intent(in)            :: omp_dev   !< Flag to activate verbose mode.
   character(*),           intent(in), optional  :: varname   !< Variable name to be printed in verbose mode.
   logical,                intent(in), optional  :: verbose   !< Flag to activate verbose mode.
   character(:), allocatable                     :: varname_  !< Variable name to be printed in verbose mode, local var.
   logical                                       :: verbose_  !< Flag to activate verbose mode, local var.
   integer(I4P)                                  :: error     !< Output of memcpy.

   varname_ = ''      ; if (present(varname)) varname_ = varname
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(lhs)) call omp_target_free_f(lhs,omp_dev)
   if (allocated(rhs)) then
      if (size(rhs, dim=1)>0) then
         call alloc_var_gpu(var=lhs, ulb=[lbound(rhs,dim=1),ubound(rhs,dim=1)], omp_dev=omp_dev, varname=varname_, verbose=verbose_)
         error = omp_target_memcpy_f(fptr_dst=lhs, fptr_src=rhs, dst_off=0_I4P, src_off=0_I4P, &
                                     omp_dst_dev=omp_dev, omp_src_dev=int(omp_get_initial_device(),I4P))
         if (error/=0) then
            print *, 'Error copying variable '//varname_
            error stop
         endif
      endif
   endif
   endsubroutine assign_allocatable_gpu_R8P_1D

   subroutine assign_allocatable_gpu_R8P_2D(lhs, rhs, omp_dev, transposed, varname, verbose)
   !< Assign GPU variable with memory checking (kind R8P, rank 2).
   !< Variable is returned not allocated if right hand side is not allocated.
   real(R8P), pointer,     intent(inout)         :: lhs(:,:)    !< Left hand side of assignement.
   real(R8P), allocatable, intent(in), target    :: rhs(:,:)    !< Right hand side of assignement.
   integer(I4P),           intent(in)            :: omp_dev !< Flag to activate verbose mode.
   logical,                intent(in), optional  :: transposed  !< Assign trasposed rhs.
   character(*),           intent(in), optional  :: varname     !< Message to be printed in verbose mode.
   logical,                intent(in), optional  :: verbose     !< Flag to activate verbose mode.
   logical                                       :: transposed_ !< Assign trasposed rhs, local var.
   integer(I4P)                                  :: ulb(2,2)    !< Upper/lower bounds of variable.
   real(R8P), allocatable                        :: rhst(:,:)   !< Right hand side transposed.
   character(:), allocatable                     :: varname_    !< Variable name to be printed in verbose mode, local var.
   logical                                       :: verbose_    !< Flag to activate verbose mode, local var.
   integer(I4P)                                  :: error       !< Output of memcpy.

   varname_ = ''      ; if (present(varname)) varname_ = varname
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(lhs)) call omp_target_free_f(lhs,omp_dev)
   if (allocated(rhs)) then
      if (size(rhs, dim=1)*size(rhs, dim=2)>0) then
         transposed_ = .false. ; if (present(transposed)) transposed_ = transposed
         ulb=reshape([lbound(rhs, dim=1),ubound(rhs, dim=1),lbound(rhs, dim=2),ubound(rhs, dim=2)], [2,2])
         if (transposed_) then
            allocate(rhst(ulb(1,2):ulb(2,2),ulb(1,1):ulb(2,1)))
            call alloc_var_gpu(var=lhs, ulb=reshape([ulb(1,2),ulb(2,2),ulb(1,1),ulb(2,1)],[2,2]), omp_dev=omp_dev, varname=varname_, verbose=verbose_)
            call transpose_a(ii=ulb(:,1), jj=ulb(:,2), a=rhs, t=rhst)
            error = omp_target_memcpy_f(fptr_dst=lhs, fptr_src=rhst, dst_off=0_I4P, src_off=0_I4P, &
                                        omp_dst_dev=omp_dev, omp_src_dev=int(omp_get_initial_device(),I4P))
            if (error/=0) then
               print *, 'Error copying variable '//varname_
               error stop
            endif
            deallocate(rhst)
         else
            call alloc_var_gpu(var=lhs, ulb=ulb, varname=varname_, omp_dev=omp_dev, verbose=verbose_)
            error = omp_target_memcpy_f(fptr_dst=lhs, fptr_src=rhs, dst_off=0_I4P, src_off=0_I4P, &
                                        omp_dst_dev=omp_dev, omp_src_dev=int(omp_get_initial_device(),I4P))
            if (error/=0) then
               print *, 'Error copying variable '//varname_
               error stop
            endif
         endif
      endif
   endif
   endsubroutine assign_allocatable_gpu_R8P_2D

   subroutine assign_allocatable_gpu_I4P_1D(lhs, rhs, omp_dev, varname, verbose)
   !< Assign GPU variable with memory checking (kind I4P, rank 1).
   !< Variable is returned not allocated if right hand side is not allocated.
   Integer(I4P), pointer,     intent(inout)         :: lhs(:)    !< Left hand side of assignement.
   Integer(I4P), allocatable, intent(in), target    :: rhs(:)    !< Right hand side of assignement.
   integer(I4P),              intent(in)            :: omp_dev   !< Flag to activate verbose mode.
   character(*),              intent(in), optional  :: varname   !< Variable name to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose   !< Flag to activate verbose mode.
   character(:), allocatable                        :: varname_  !< Variable name to be printed in verbose mode, local var.
   logical                                          :: verbose_  !< Flag to activate verbose mode, local var.
   integer(I4P)                                     :: error     !< Output of memcpy.

   varname_ = ''      ; if (present(varname)) varname_ = varname
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(lhs)) call omp_target_free_f(lhs,omp_dev)
   if (allocated(rhs)) then
      if (size(rhs, dim=1)>0) then
         call alloc_var_gpu(var=lhs, ulb=[lbound(rhs,dim=1),ubound(rhs,dim=1)], omp_dev=omp_dev, varname=varname_, verbose=verbose_)
         error = omp_target_memcpy_f(fptr_dst=lhs, fptr_src=rhs, dst_off=0_I4P, src_off=0_I4P, &
                                     omp_dst_dev=omp_dev, omp_src_dev=int(omp_get_initial_device(),I4P))
         if (error/=0) then
            print *, 'Error copying variable '//varname_
            error stop
         endif
      endif
   endif
   endsubroutine assign_allocatable_gpu_I4P_1D

   subroutine assign_allocatable_gpu_I8P_2D(lhs, rhs, omp_dev, varname, verbose)
   !< Assign GPU variable with memory checking (kind I8P, rank 2).
   !< Variable is returned not allocated if right hand side is not allocated.
   Integer(I8P), pointer,     intent(inout)         :: lhs(:,:) !< Left hand side of assignement.
   Integer(I8P), allocatable, intent(in), target    :: rhs(:,:) !< Right hand side of assignement.
   integer(I4P),              intent(in)            :: omp_dev  !< Flag to activate verbose mode.
   character(*),              intent(in), optional  :: varname  !< Variable name to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose  !< Flag to activate verbose mode.
   character(:), allocatable                        :: varname_ !< Variable name to be printed in verbose mode, local var.
   logical                                          :: verbose_ !< Flag to activate verbose mode, local var.
   integer(I4P)                                     :: error    !< Output of memcpy.

   varname_ = ''      ; if (present(varname)) varname_ = varname
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(lhs)) call omp_target_free_f(lhs,omp_dev)
   if (allocated(rhs)) then
      if (size(rhs, dim=1)*size(rhs, dim=2)>0) then
         call alloc_var_gpu(var=lhs,                                                                                     &
                            ulb=reshape([lbound(rhs,dim=1),ubound(rhs,dim=1),lbound(rhs,dim=2),ubound(rhs,dim=2)],[2,2]),&
                            omp_dev=omp_dev, varname=varname_, verbose=verbose_)
         error = omp_target_memcpy_f(fptr_dst=lhs, fptr_src=rhs, dst_off=0_I4P, src_off=0_I4P, &
                                     omp_dst_dev=omp_dev, omp_src_dev=int(omp_get_initial_device(),I4P))
         if (error/=0) then
            print *, 'Error copying variable '//varname_
            error stop
         endif
      endif
   endif
   endsubroutine assign_allocatable_gpu_I8P_2D

   subroutine assign_allocatable_gpu_I8P_3D(lhs, rhs, omp_dev, varname, verbose)
   !< Assign GPU variable with memory checking (kind I8P, rank 3).
   !< Variable is returned not allocated if right hand side is not allocated.
   Integer(I8P), pointer,     intent(inout)         :: lhs(:,:,:) !< Left hand side of assignement.
   Integer(I8P), allocatable, intent(in), target    :: rhs(:,:,:) !< Right hand side of assignement.
   integer(I4P),              intent(in)            :: omp_dev    !< Flag to activate verbose mode.
   character(*),              intent(in), optional  :: varname    !< Variable name to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose    !< Flag to activate verbose mode.
   character(:), allocatable                        :: varname_   !< Variable name to be printed in verbose mode, local var.
   logical                                          :: verbose_   !< Flag to activate verbose mode, local var.
   integer(I4P)                                     :: error      !< Output of memcpy.

   varname_ = ''      ; if (present(varname)) varname_ = varname
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(lhs)) call omp_target_free_f(lhs,omp_dev)
   if (allocated(rhs)) then
      if (size(rhs, dim=1)*size(rhs, dim=2)*size(rhs, dim=3)>0) then
         call alloc_var_gpu(var=lhs,                                                  &
                            ulb=reshape([lbound(rhs,dim=1),ubound(rhs,dim=1),         &
                                         lbound(rhs,dim=2),ubound(rhs,dim=2),         &
                                         lbound(rhs,dim=3),ubound(rhs,dim=3)],[2,3]), &
                            omp_dev=omp_dev, varname=varname_, verbose=verbose_)
         error = omp_target_memcpy_f(fptr_dst=lhs, fptr_src=rhs, dst_off=0_I4P, src_off=0_I4P, &
                                     omp_dst_dev=omp_dev, omp_src_dev=int(omp_get_initial_device(),I4P))
         if (error/=0) then
            print *, 'Error copying variable '//varname_
            error stop
         endif
      endif
   endif
   endsubroutine assign_allocatable_gpu_I8P_3D

   !subroutine save_memory_gpu_status(file_name, tag)
   !!< Save the current GPU-memory status into a file.
   !!< File is accessed in append position.
   !character(*), intent(in)           :: file_name           !< File name.
   !character(*), intent(in), optional :: tag                 !< Tag of current status.
   !character(:), allocatable          :: tag_                !< Tag of current status, local var.
   !integer(cuda_count_kind)           :: mem_free, mem_total !< Device memory.
   !integer(I4P)                       :: file_unit           !< File unit.
   !integer(I4P)                       :: error               !< Error traping flag.

   !tag_ = '' ; if (present(tag)) tag_ = trim(tag)
   !error = cudaMemGetInfo(mem_free, mem_total)
   !open(newunit=file_unit, file=trim(file_name), position="append")
   !write(file_unit,*) tag_, mem_free, mem_total
   !close(file_unit)
   !endsubroutine save_memory_gpu_status

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
endmodule adam_memory_gmp_lib
