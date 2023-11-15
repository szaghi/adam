!< ADAM, GMP memory handling library.
module adam_memory_gmp_library
!< ADAM, GMP memory handling library.

use penf
use adam_gmp_utils
use omp_lib, only : omp_get_initial_device

implicit none
save
private
public :: alloc_var_gpu
public :: assign_allocatable_gpu
public :: save_memory_gpu_status

interface alloc_var_gpu
!< Allocate GPU variable with memory checking.
module procedure alloc_var_gpu_R8P_1D, &
                 alloc_var_gpu_R8P_2D, &
                 alloc_var_gpu_R8P_3D, &
                 alloc_var_gpu_R8P_4D, &
                 alloc_var_gpu_R8P_5D, &
                 alloc_var_gpu_R8P_6D, &
                 alloc_var_gpu_I4P_1D, &
                 alloc_var_gpu_I4P_2D, &
                 alloc_var_gpu_I4P_5D, &
                 alloc_var_gpu_I8P_1D, &
                 alloc_var_gpu_I8P_2D, &
                 alloc_var_gpu_I8P_3D
endinterface alloc_var_gpu

interface assign_allocatable_gpu
!< Assign GPU variable with memory checking.
module procedure  assign_allocatable_gpu_R8P_1D,               &
                  assign_allocatable_gpu_R8P_2D,               &
                  assign_allocatable_gpu_R8P_3D,               &
                  assign_allocatable_gpu_R8P_4D,               &
                  assign_allocatable_gpu_I4P_1D,               &
                  ! assign_allocatable_gpu_I4P_1D_rhs_allocated, &
                  assign_allocatable_gpu_I4P_5D,               &
                  assign_allocatable_gpu_I8P_2D,               &
                  assign_allocatable_gpu_I8P_3D
endinterface assign_allocatable_gpu

interface transpose_a
   module procedure transpose_a_R8P_2D !< Transpose array (kind R8P, rank 2).
endinterface transpose_a

contains
   subroutine alloc_var_gpu_R8P_1D(var, ulb, omp_dev, init_val, msg, verbose)
   !< Allocate GPU variable with memory checking (kind R8P, rank 1).
   real(R8P), pointer, intent(inout)         :: var(:)   !< Varibale to be allocate on GPU.
   integer(I4P),       intent(in)            :: ulb(2)   !< Upper/lower bounds of variable.
   integer(I4P),       intent(in)            :: omp_dev  !< OpenMP device ID.
   real(R8P),          intent(in), optional  :: init_val !< Initialization value.
   character(*),       intent(in), optional  :: msg      !< Message to be printed in verbose mode.
   logical,            intent(in), optional  :: verbose  !< Flag to activate verbose mode.
   character(:), allocatable                 :: msg_     !< Message to be printed in verbose mode, local var.
   logical                                   :: verbose_ !< Flag to activate verbose mode, local var.
   integer(I4P)                              :: error    !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   if (present(init_val)) then
      call omp_target_alloc_f(fptr_dev=var, ubounds=[ulb(2)], lbounds=[ulb(1)], omp_dev=omp_dev, init_value=init_val, ierr=error)
   else
      call omp_target_alloc_f(fptr_dev=var, ubounds=[ulb(2)], lbounds=[ulb(1)], omp_dev=omp_dev, ierr=error)
   endif
   if (error/=0) then
      print '(A)', 'Error in alloc_var_gpu_R8P_1D: '//msg_
      error stop
   else
      if (verbose_) print '(A)', msg_
   endif
   endsubroutine alloc_var_gpu_R8P_1D

   subroutine alloc_var_gpu_R8P_2D(var, ulb, omp_dev, init_val, msg, verbose)
   !< Allocate GPU variable with memory checking (kind R8P, rank 2).
   real(R8P), pointer, intent(inout)         :: var(:,:) !< Varibale to be allocate on GPU.
   integer(I4P),       intent(in)            :: ulb(2,2) !< Upper/lower bounds of variable.
   integer(I4P),       intent(in)            :: omp_dev  !< OpenMP device ID.
   real(R8P),          intent(in), optional  :: init_val !< Initialization value.
   character(*),       intent(in), optional  :: msg      !< Message to be printed in verbose mode.
   logical,            intent(in), optional  :: verbose  !< Flag to activate verbose mode.
   character(:), allocatable                 :: msg_     !< Message to be printed in verbose mode, local var.
   logical                                   :: verbose_ !< Flag to activate verbose mode, local var.
   integer(I4P)                              :: error    !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   if (present(init_val)) then
      call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, init_value=init_val, ierr=error)
   else
      call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, ierr=error)
   endif
   if (error/=0) then
      print '(A)', 'Error in alloc_var_gpu_R8P_2D: '//msg_
      error stop
   else
      if (verbose_) print '(A)', msg_
   endif
   endsubroutine alloc_var_gpu_R8P_2D

   subroutine alloc_var_gpu_R8P_3D(var, ulb, omp_dev, init_val, msg, verbose)
   !< Allocate GPU variable with memory checking (kind R8P, rank 3).
   real(R8P), pointer, intent(inout)         :: var(:,:,:) !< Varibale to be allocate on GPU.
   integer(I4P),       intent(in)            :: ulb(2,3)   !< Upper/lower bounds of variable.
   integer(I4P),       intent(in)            :: omp_dev    !< OpenMP device ID.
   real(R8P),          intent(in), optional  :: init_val   !< Initialization value.
   character(*),       intent(in), optional  :: msg        !< Message to be printed in verbose mode.
   logical,            intent(in), optional  :: verbose    !< Flag to activate verbose mode.
   character(:), allocatable                 :: msg_       !< Message to be printed in verbose mode, local var.
   logical                                   :: verbose_   !< Flag to activate verbose mode, local var.
   integer(I4P)                              :: error      !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   if (present(init_val)) then
      call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, init_value=init_val, ierr=error)
   else
      call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, ierr=error)
   endif
   if (error/=0) then
      print '(A)', 'Error in alloc_var_gpu_R8P_3D: '//msg_
      error stop
   else
      if (verbose_) print '(A)', msg_
   endif
   endsubroutine alloc_var_gpu_R8P_3D

   subroutine alloc_var_gpu_R8P_4D(var, ulb, omp_dev, init_val, msg, verbose)
   !< Allocate GPU variable with memory checking (kind R8P, rank 4).
   real(R8P), pointer, intent(inout)         :: var(:,:,:,:) !< Varibale to be allocate on GPU.
   integer(I4P),       intent(in)            :: ulb(2,4)     !< Upper/lower bounds of variable.
   integer(I4P),       intent(in)            :: omp_dev      !< OpenMP device ID.
   real(R8P),          intent(in), optional  :: init_val     !< Initialization value.
   character(*),       intent(in), optional  :: msg          !< Message to be printed in verbose mode.
   logical,            intent(in), optional  :: verbose      !< Flag to activate verbose mode.
   character(:), allocatable                 :: msg_         !< Message to be printed in verbose mode, local var.
   logical                                   :: verbose_     !< Flag to activate verbose mode, local var.
   integer(I4P)                              :: error        !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   if (present(init_val)) then
      call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, init_value=init_val, ierr=error)
   else
      call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, ierr=error)
   endif
   if (error/=0) then
      print '(A)', 'Error in alloc_var_gpu_R8P_4D: '//msg_
      error stop
   else
      if (verbose_) print '(A)', msg_
   endif
   endsubroutine alloc_var_gpu_R8P_4D

   subroutine alloc_var_gpu_R8P_5D(var, ulb, omp_dev, init_val, msg, verbose)
   !< Allocate GPU variable with memory checking (kind R8P, rank 5).
   real(R8P), pointer, intent(inout)         :: var(:,:,:,:,:) !< Varibale to be allocate on GPU.
   integer(I4P),       intent(in)            :: ulb(2,5)       !< Upper/lower bounds of variable.
   integer(I4P),       intent(in)            :: omp_dev        !< OpenMP device ID.
   real(R8P),          intent(in), optional  :: init_val       !< Initialization value.
   character(*),       intent(in), optional  :: msg            !< Message to be printed in verbose mode.
   logical,            intent(in), optional  :: verbose        !< Flag to activate verbose mode.
   character(:), allocatable                 :: msg_           !< Message to be printed in verbose mode, local var.
   logical                                   :: verbose_       !< Flag to activate verbose mode, local var.
   integer(I4P)                              :: error          !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   if (present(init_val)) then
      call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, init_value=init_val, ierr=error)
   else
      call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, ierr=error)
   endif
   if (error/=0) then
      print '(A)', 'Error in alloc_var_gpu_R8P_5D: '//msg_
      error stop
   else
      if (verbose_) print '(A)', msg_
   endif
   endsubroutine alloc_var_gpu_R8P_5D

   subroutine alloc_var_gpu_R8P_6D(var, ulb, omp_dev, init_val, msg, verbose)
   !< Allocate GPU variable with memory checking (kind R8P, rank 6).
   real(R8P), pointer, intent(inout)         :: var(:,:,:,:,:,:) !< Varibale to be allocate on GPU.
   integer(I4P),       intent(in)            :: ulb(2,6)         !< Upper/lower bounds of variable.
   integer(I4P),       intent(in)            :: omp_dev          !< OpenMP device ID.
   real(R8P),          intent(in), optional  :: init_val         !< Initialization value.
   character(*),       intent(in), optional  :: msg              !< Message to be printed in verbose mode.
   logical,            intent(in), optional  :: verbose          !< Flag to activate verbose mode.
   character(:), allocatable                 :: msg_             !< Message to be printed in verbose mode, local var.
   logical                                   :: verbose_         !< Flag to activate verbose mode, local var.
   integer(I4P)                              :: error            !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   if (present(init_val)) then
      call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, init_value=init_val, ierr=error)
   else
      call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, ierr=error)
   endif
   if (error/=0) then
      print '(A)', 'Error in alloc_var_gpu_R8P_6D: '//msg_
      error stop
   else
      if (verbose_) print '(A)', msg_
   endif
   endsubroutine alloc_var_gpu_R8P_6D

   subroutine alloc_var_gpu_I4P_1D(var, ulb, omp_dev, init_val, msg, verbose)
   !< Allocate GPU variable with memory checking (kind I4P, rank 1).
   integer(I4P), pointer, intent(inout)         :: var(:)   !< Varibale to be allocate on GPU.
   integer(I4P),          intent(in)            :: ulb(2)   !< Upper/lower bounds of variable.
   integer(I4P),          intent(in)            :: omp_dev  !< OpenMP device ID.
   integer(I4P),          intent(in), optional  :: init_val !< Initialization value.
   character(*),          intent(in), optional  :: msg      !< Message to be printed in verbose mode.
   logical,               intent(in), optional  :: verbose  !< Flag to activate verbose mode.
   character(:), allocatable                    :: msg_     !< Message to be printed in verbose mode, local var.
   logical                                      :: verbose_ !< Flag to activate verbose mode, local var.
   integer(I4P)                                 :: error    !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   if (present(init_val)) then
      call omp_target_alloc_f(fptr_dev=var, ubounds=[ulb(2)], lbounds=[ulb(1)], omp_dev=omp_dev, init_value=init_val, ierr=error)
   else
      call omp_target_alloc_f(fptr_dev=var, ubounds=[ulb(2)], lbounds=[ulb(1)], omp_dev=omp_dev, ierr=error)
   endif
   if (error/=0) then
      print '(A)', 'Error in alloc_var_gpu_I4P_1D: '//msg_
      error stop
   else
      if (verbose_) print '(A)', msg_
   endif
   endsubroutine alloc_var_gpu_I4P_1D

   subroutine alloc_var_gpu_I4P_2D(var, ulb, omp_dev, init_val, msg, verbose)
   !< Allocate GPU variable with memory checking (kind I4P, rank 2).
   integer(I4P), pointer, intent(inout)         :: var(:,:) !< Varibale to be allocate on GPU.
   integer(I4P),          intent(in)            :: ulb(2,2) !< Upper/lower bounds of variable.
   integer(I4P),          intent(in)            :: omp_dev  !< OpenMP device ID.
   integer(I4P),          intent(in), optional  :: init_val !< Initialization value.
   character(*),          intent(in), optional  :: msg      !< Message to be printed in verbose mode.
   logical,               intent(in), optional  :: verbose  !< Flag to activate verbose mode.
   character(:), allocatable                    :: msg_     !< Message to be printed in verbose mode, local var.
   logical                                      :: verbose_ !< Flag to activate verbose mode, local var.
   integer(I4P)                                 :: error    !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   if (present(init_val)) then
      call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, init_value=init_val, ierr=error)
   else
      call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, ierr=error)
   endif
   if (error/=0) then
      print '(A)', 'Error in alloc_var_gpu_I4P_2D: '//msg_
      error stop
   else
      if (verbose_) print '(A)', msg_
   endif
   endsubroutine alloc_var_gpu_I4P_2D

   subroutine alloc_var_gpu_I4P_5D(var, ulb, omp_dev, init_val, msg, verbose)
   !< Allocate GPU variable with memory checking (kind I4P, rank 5).
   integer(I4P), pointer, intent(inout)         :: var(:,:,:,:,:) !< Varibale to be allocate on GPU.
   integer(I4P),          intent(in)            :: ulb(2,5)       !< Upper/lower bounds of variable.
   integer(I4P),          intent(in)            :: omp_dev        !< OpenMP device ID.
   integer(I4P),          intent(in), optional  :: init_val       !< Initialization value.
   character(*),          intent(in), optional  :: msg            !< Message to be printed in verbose mode.
   logical,               intent(in), optional  :: verbose        !< Flag to activate verbose mode.
   character(:), allocatable                    :: msg_           !< Message to be printed in verbose mode, local var.
   logical                                      :: verbose_       !< Flag to activate verbose mode, local var.
   integer(I4P)                                 :: error          !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   if (present(init_val)) then
      call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, init_value=init_val, ierr=error)
   else
      call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, ierr=error)
   endif
   if (error/=0) then
      print '(A)', 'Error in alloc_var_gpu_I4P_5D: '//msg_
      error stop
   else
      if (verbose_) print '(A)', msg_
   endif
   endsubroutine alloc_var_gpu_I4P_5D

   subroutine alloc_var_gpu_I8P_1D(var, ulb, omp_dev, init_val, msg, verbose)
   !< Allocate GPU variable with memory checking (kind I8P, rank 1).
   integer(I8P), pointer, intent(inout)         :: var(:)   !< Varibale to be allocate on GPU.
   integer(I4P),          intent(in)            :: ulb(2)   !< Upper/lower bounds of variable.
   integer(I4P),          intent(in)            :: omp_dev  !< OpenMP device ID.
   integer(I8P),          intent(in), optional  :: init_val !< Initialization value.
   character(*),          intent(in), optional  :: msg      !< Message to be printed in verbose mode.
   logical,               intent(in), optional  :: verbose  !< Flag to activate verbose mode.
   character(:), allocatable                    :: msg_     !< Message to be printed in verbose mode, local var.
   logical                                      :: verbose_ !< Flag to activate verbose mode, local var.
   integer(I4P)                                 :: error    !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   if (present(init_val)) then
      call omp_target_alloc_f(fptr_dev=var, ubounds=[ulb(2)], lbounds=[ulb(1)], omp_dev=omp_dev, init_value=init_val, ierr=error)
   else
      call omp_target_alloc_f(fptr_dev=var, ubounds=[ulb(2)], lbounds=[ulb(1)], omp_dev=omp_dev, ierr=error)
   endif
   if (error/=0) then
      print '(A)', 'Error in alloc_var_gpu_I8P_1D: '//msg_
      error stop
   else
      if (verbose_) print '(A)', msg_
   endif
   endsubroutine alloc_var_gpu_I8P_1D

   subroutine alloc_var_gpu_I8P_2D(var, ulb, omp_dev, init_val, msg, verbose)
   !< Allocate GPU variable with memory checking (kind I8P, rank 2).
   integer(I8P), pointer, intent(inout)         :: var(:,:) !< Varibale to be allocate on GPU.
   integer(I4P),          intent(in)            :: ulb(2,2) !< Upper/lower bounds of variable.
   integer(I4P),          intent(in)            :: omp_dev  !< OpenMP device ID.
   integer(I8P),          intent(in), optional  :: init_val !< Initialization value.
   character(*),          intent(in), optional  :: msg      !< Message to be printed in verbose mode.
   logical,               intent(in), optional  :: verbose  !< Flag to activate verbose mode.
   character(:), allocatable                    :: msg_     !< Message to be printed in verbose mode, local var.
   logical                                      :: verbose_ !< Flag to activate verbose mode, local var.
   integer(I4P)                                 :: error    !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   if (present(init_val)) then
      call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, init_value=init_val, ierr=error)
   else
      call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=[ulb(1,:)], omp_dev=omp_dev, ierr=error)
   endif
   if (error/=0) then
      print '(A)', 'Error in alloc_var_gpu_I8P_2D: '//msg_
      error stop
   else
      if (verbose_) print '(A)', msg_
   endif
   endsubroutine alloc_var_gpu_I8P_2D

   subroutine alloc_var_gpu_I8P_3D(var, ulb, omp_dev, init_val, msg, verbose)
   !< Allocate GPU variable with memory checking (kind I8P, rank 2).
   integer(I8P), pointer, intent(inout)         :: var(:,:,:) !< Varibale to be allocate on GPU.
   integer(I4P),          intent(in)            :: ulb(2,3)   !< Upper/lower bounds of variable.
   integer(I4P),          intent(in)            :: omp_dev    !< OpenMP device ID.
   integer(I8P),          intent(in), optional  :: init_val   !< Initialization value.
   character(*),          intent(in), optional  :: msg        !< Message to be printed in verbose mode.
   logical,               intent(in), optional  :: verbose    !< Flag to activate verbose mode.
   character(:), allocatable                    :: msg_       !< Message to be printed in verbose mode, local var.
   logical                                      :: verbose_   !< Flag to activate verbose mode, local var.
   integer(I4P)                                 :: error      !< Error traping flag.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(var)) call omp_target_free_f(var,omp_dev)
   if (present(init_val)) then
      call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=ulb(1,:), omp_dev=omp_dev, init_value=init_val, ierr=error)
   else
      call omp_target_alloc_f(fptr_dev=var, ubounds=ulb(2,:), lbounds=ulb(1,:), omp_dev=omp_dev, ierr=error)
   endif
   if (error/=0) then
      print '(A)', 'Error in alloc_var_gpu_I8P_3D: '//msg_
      error stop
   else
      if (verbose_) print '(A)', msg_
   endif
   endsubroutine alloc_var_gpu_I8P_3D

   subroutine assign_allocatable_gpu_R8P_1D(lhs, rhs, omp_dev, msg, verbose)
   !< Assign GPU variable with memory checking (kind R8P, rank 1).
   !< Variable is returned not allocated if right hand side is not allocated.
   real(R8P), pointer,     intent(inout)         :: lhs(:)    !< Left hand side of assignement.
   real(R8P), allocatable, intent(in), target    :: rhs(:)    !< Right hand side of assignement.
   integer(I4P),           intent(in)            :: omp_dev   !< Flag to activate verbose mode.
   character(*),           intent(in), optional  :: msg       !< Message to be printed in verbose mode.
   logical,                intent(in), optional  :: verbose   !< Flag to activate verbose mode.
   character(:), allocatable                     :: msg_      !< Message to be printed in verbose mode, local var.
   logical                                       :: verbose_  !< Flag to activate verbose mode, local var.
   integer(I4P)                                  :: error     !< Output of memcpy.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(lhs)) call omp_target_free_f(lhs,omp_dev)
   if (allocated(rhs)) then
      if (size(rhs, dim=1)>0) then
         call alloc_var_gpu(var=lhs, ulb=[lbound(rhs,dim=1),ubound(rhs,dim=1)], omp_dev=omp_dev, msg=msg_, verbose=verbose_)
         error = omp_target_memcpy_f(fptr_dst=lhs, fptr_src=rhs, dst_off=0_I4P, src_off=0_I4P, &
                                     omp_dst_dev=omp_dev, omp_src_dev=int(omp_get_initial_device(),I4P))
         if (error/=0) then
            print '(A)', 'Error in assign_allocatable_gpu_R8P_1D (omp_target_memcpy_f): '//msg_
            error stop
         endif
      endif
   endif
   endsubroutine assign_allocatable_gpu_R8P_1D

   subroutine assign_allocatable_gpu_R8P_2D(lhs, rhs, omp_dev, transposed, msg, verbose)
   !< Assign GPU variable with memory checking (kind R8P, rank 2).
   !< Variable is returned not allocated if right hand side is not allocated.
   real(R8P), pointer,     intent(inout)         :: lhs(:,:)    !< Left hand side of assignement.
   real(R8P), allocatable, intent(in), target    :: rhs(:,:)    !< Right hand side of assignement.
   integer(I4P),           intent(in)            :: omp_dev     !< Flag to activate verbose mode.
   logical,                intent(in), optional  :: transposed  !< Assign trasposed rhs.
   character(*),           intent(in), optional  :: msg         !< Message to be printed in verbose mode.
   logical,                intent(in), optional  :: verbose     !< Flag to activate verbose mode.
   logical                                       :: transposed_ !< Assign trasposed rhs, local var.
   integer(I4P)                                  :: ulb(2,2)    !< Upper/lower bounds of variable.
   real(R8P), allocatable                        :: rhst(:,:)   !< Right hand side transposed.
   character(:), allocatable                     :: msg_        !< Message to be printed in verbose mode, local var.
   logical                                       :: verbose_    !< Flag to activate verbose mode, local var.
   integer(I4P)                                  :: error       !< Output of memcpy.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(lhs)) call omp_target_free_f(lhs,omp_dev)
   if (allocated(rhs)) then
      if (size(rhs, dim=1)*size(rhs, dim=2)>0) then
         transposed_ = .false. ; if (present(transposed)) transposed_ = transposed
         ulb=reshape([lbound(rhs, dim=1),ubound(rhs, dim=1),lbound(rhs, dim=2),ubound(rhs, dim=2)], [2,2])
         if (transposed_) then
            allocate(rhst(ulb(1,2):ulb(2,2),ulb(1,1):ulb(2,1)))
            call alloc_var_gpu(var=lhs, ulb=reshape([ulb(1,2),ulb(2,2),ulb(1,1),ulb(2,1)],[2,2]), omp_dev=omp_dev, &
                               msg=msg_, verbose=verbose_)
            call transpose_a(ii=ulb(:,1), jj=ulb(:,2), a=rhs, t=rhst)
            error = omp_target_memcpy_f(fptr_dst=lhs, fptr_src=rhst, dst_off=0_I4P, src_off=0_I4P, &
                                        omp_dst_dev=omp_dev, omp_src_dev=int(omp_get_initial_device(),I4P))
            if (error/=0) then
               print '(A)', 'Error in assign_allocatable_gpu_R8P_2D (transposed, omp_target_memcpy_f): '//msg_
               error stop
            endif
            deallocate(rhst)
         else
            call alloc_var_gpu(var=lhs, ulb=ulb, msg=msg_, omp_dev=omp_dev, verbose=verbose_)
            error = omp_target_memcpy_f(fptr_dst=lhs, fptr_src=rhs, dst_off=0_I4P, src_off=0_I4P, &
                                        omp_dst_dev=omp_dev, omp_src_dev=int(omp_get_initial_device(),I4P))
            if (error/=0) then
               print '(A)', 'Error in assign_allocatable_gpu_R8P_2D (omp_target_memcpy_f): '//msg_
               error stop
            endif
         endif
      endif
   endif
   endsubroutine assign_allocatable_gpu_R8P_2D

   subroutine assign_allocatable_gpu_R8P_3D(lhs, rhs, omp_dev, msg, verbose)
   !< Assign GPU variable with memory checking (kind R8P, rank 3).
   !< Variable is returned not allocated if right hand side is not allocated.
   real(R8P), pointer,     intent(inout)         :: lhs(:,:,:) !< Left hand side of assignement.
   real(R8P), allocatable, intent(in), target    :: rhs(:,:,:) !< Right hand side of assignement.
   integer(I4P),           intent(in)            :: omp_dev    !< Flag to activate verbose mode.
   character(*),           intent(in), optional  :: msg        !< Message to be printed in verbose mode.
   logical,                intent(in), optional  :: verbose    !< Flag to activate verbose mode.
   character(:), allocatable                     :: msg_       !< Message to be printed in verbose mode, local var.
   logical                                       :: verbose_   !< Flag to activate verbose mode, local var.
   integer(I4P)                                  :: error      !< Output of memcpy.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(lhs)) call omp_target_free_f(lhs,omp_dev)
   if (allocated(rhs)) then
      if (size(rhs,dim=1)*size(rhs,dim=2)*size(rhs,dim=3)>0) then
         call alloc_var_gpu(var=lhs, &
                            ulb=reshape([lbound(rhs,dim=1),ubound(rhs,dim=1),         &
                                         lbound(rhs,dim=2),ubound(rhs,dim=2),         &
                                         lbound(rhs,dim=3),ubound(rhs,dim=3)],[2,3]), &
                            omp_dev=omp_dev, msg=msg_, verbose=verbose_)
         error = omp_target_memcpy_f(fptr_dst=lhs, fptr_src=rhs, dst_off=0_I4P, src_off=0_I4P, &
                                     omp_dst_dev=omp_dev, omp_src_dev=int(omp_get_initial_device(),I4P))
         if (error/=0) then
            print '(A)', 'Error in assign_allocatable_gpu_R8P_3D (omp_target_memcpy_f): '//msg_
            error stop
         endif
      endif
   endif
   endsubroutine assign_allocatable_gpu_R8P_3D

   subroutine assign_allocatable_gpu_R8P_4D(lhs, rhs, omp_dev, msg, verbose)
   !< Assign GPU variable with memory checking (kind R8P, rank 4).
   !< Variable is returned not allocated if right hand side is not allocated.
   real(R8P), pointer,     intent(inout)         :: lhs(:,:,:,:) !< Left hand side of assignement.
   real(R8P), allocatable, intent(in), target    :: rhs(:,:,:,:) !< Right hand side of assignement.
   integer(I4P),           intent(in)            :: omp_dev      !< Flag to activate verbose mode.
   character(*),           intent(in), optional  :: msg          !< Message to be printed in verbose mode.
   logical,                intent(in), optional  :: verbose      !< Flag to activate verbose mode.
   character(:), allocatable                     :: msg_         !< Message to be printed in verbose mode, local var.
   logical                                       :: verbose_     !< Flag to activate verbose mode, local var.
   integer(I4P)                                  :: error        !< Output of memcpy.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(lhs)) call omp_target_free_f(lhs,omp_dev)
   if (allocated(rhs)) then
      if (size(rhs,dim=1)*size(rhs,dim=2)*size(rhs,dim=3)*size(rhs,dim=4)>0) then
         call alloc_var_gpu(var=lhs, &
                            ulb=reshape([lbound(rhs,dim=1),ubound(rhs,dim=1),         &
                                         lbound(rhs,dim=2),ubound(rhs,dim=2),         &
                                         lbound(rhs,dim=3),ubound(rhs,dim=3),         &
                                         lbound(rhs,dim=4),ubound(rhs,dim=4)],[2,4]), &
                            omp_dev=omp_dev, msg=msg_, verbose=verbose_)
         error = omp_target_memcpy_f(fptr_dst=lhs, fptr_src=rhs, dst_off=0_I4P, src_off=0_I4P, &
                                     omp_dst_dev=omp_dev, omp_src_dev=int(omp_get_initial_device(),I4P))
         if (error/=0) then
            print '(A)', 'Error in assign_allocatable_gpu_R8P_4D (omp_target_memcpy_f): '//msg_
            error stop
         endif
      endif
   endif
   endsubroutine assign_allocatable_gpu_R8P_4D

   subroutine assign_allocatable_gpu_I4P_1D(lhs, rhs, omp_dev, msg, verbose)
   !< Assign GPU variable with memory checking (kind I4P, rank 1).
   !< Variable is returned not allocated if right hand side is not allocated.
   Integer(I4P), pointer,     intent(inout)         :: lhs(:)    !< Left hand side of assignement.
   Integer(I4P), allocatable, intent(in), target    :: rhs(:)    !< Right hand side of assignement.
   integer(I4P),              intent(in)            :: omp_dev   !< Flag to activate verbose mode.
   character(*),              intent(in), optional  :: msg       !< Message to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose   !< Flag to activate verbose mode.
   character(:), allocatable                        :: msg_      !< Message to be printed in verbose mode, local var.
   logical                                          :: verbose_  !< Flag to activate verbose mode, local var.
   integer(I4P)                                     :: error     !< Output of memcpy.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(lhs)) call omp_target_free_f(lhs,omp_dev)
   if (allocated(rhs)) then
      if (size(rhs, dim=1)>0) then
         call alloc_var_gpu(var=lhs, ulb=[lbound(rhs,dim=1),ubound(rhs,dim=1)], omp_dev=omp_dev, msg=msg_, verbose=verbose_)
         error = omp_target_memcpy_f(fptr_dst=lhs, fptr_src=rhs, dst_off=0_I4P, src_off=0_I4P, &
                                     omp_dst_dev=omp_dev, omp_src_dev=int(omp_get_initial_device(),I4P))
         if (error/=0) then
            print '(A)', 'Error in assign_allocatable_gpu_I4P_1D (omp_target_memcpy_f): '//msg_
            error stop
         endif
      endif
   endif
   endsubroutine assign_allocatable_gpu_I4P_1D

   subroutine assign_allocatable_gpu_I4P_1D_rhs_allocated(lhs, rhsa, omp_dev, msg, verbose)
   !< Assign GPU variable with memory checking (kind I4P, rank 1).
   !< Variable is returned not allocated if right hand side is not allocated.
   Integer(I4P), pointer, intent(inout)         :: lhs(:)   !< Left hand side of assignement.
   Integer(I4P),          intent(in), target    :: rhsa(:)  !< Right hand side of assignement.
   integer(I4P),          intent(in)            :: omp_dev  !< Flag to activate verbose mode.
   character(*),          intent(in), optional  :: msg      !< Message to be printed in verbose mode.
   logical,               intent(in), optional  :: verbose  !< Flag to activate verbose mode.
   character(:), allocatable                    :: msg_     !< Message to be printed in verbose mode, local var.
   logical                                      :: verbose_ !< Flag to activate verbose mode, local var.
   integer(I4P)                                 :: error    !< Output of memcpy.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(lhs)) call omp_target_free_f(lhs,omp_dev)
   if (size(rhsa, dim=1)>0) then
      call alloc_var_gpu(var=lhs, ulb=[lbound(rhsa,dim=1),ubound(rhsa,dim=1)], omp_dev=omp_dev, msg=msg_, verbose=verbose_)
      error = omp_target_memcpy_f(fptr_dst=lhs, fptr_src=rhsa, dst_off=0_I4P, src_off=0_I4P, &
                                  omp_dst_dev=omp_dev, omp_src_dev=int(omp_get_initial_device(),I4P))
      if (error/=0) then
         print '(A)', 'Error in assign_allocatable_gpu_I4P_1D_rhs_allocated (omp_target_memcpy_f): '//msg_
         error stop
      endif
   endif
   endsubroutine assign_allocatable_gpu_I4P_1D_rhs_allocated

   subroutine assign_allocatable_gpu_I4P_5D(lhs, rhs, omp_dev, msg, verbose)
   !< Assign GPU variable with memory checking (kind I4P, rank 5).
   !< Variable is returned not allocated if right hand side is not allocated.
   integer(I4P), pointer,     intent(inout)         :: lhs(:,:,:,:,:) !< Left hand side of assignement.
   integer(I4P), allocatable, intent(in), target    :: rhs(:,:,:,:,:) !< Right hand side of assignement.
   integer(I4P),              intent(in)            :: omp_dev        !< Flag to activate verbose mode.
   character(*),              intent(in), optional  :: msg            !< Message to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose        !< Flag to activate verbose mode.
   character(:), allocatable                        :: msg_           !< Message to be printed in verbose mode, local var.
   logical                                          :: verbose_       !< Flag to activate verbose mode, local var.
   integer(I4P)                                     :: error          !< Output of memcpy.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(lhs)) call omp_target_free_f(lhs,omp_dev)
   if (allocated(rhs)) then
      if (size(rhs,dim=1)*size(rhs,dim=2)*size(rhs,dim=3)*size(rhs,dim=4)*size(rhs,dim=5)>0) then
         call alloc_var_gpu(var=lhs, &
                            ulb=reshape([lbound(rhs,dim=1),ubound(rhs,dim=1),         &
                                         lbound(rhs,dim=2),ubound(rhs,dim=2),         &
                                         lbound(rhs,dim=3),ubound(rhs,dim=3),         &
                                         lbound(rhs,dim=4),ubound(rhs,dim=4),         &
                                         lbound(rhs,dim=5),ubound(rhs,dim=5)],[2,5]), &
                            omp_dev=omp_dev, msg=msg_, verbose=verbose_)
         error = omp_target_memcpy_f(fptr_dst=lhs, fptr_src=rhs, dst_off=0_I4P, src_off=0_I4P, &
                                     omp_dst_dev=omp_dev, omp_src_dev=int(omp_get_initial_device(),I4P))
         if (error/=0) then
            print '(A)', 'Error in assign_allocatable_gpu_I4P_5D (omp_target_memcpy_f): '//msg_
            error stop
         endif
      endif
   endif
   endsubroutine assign_allocatable_gpu_I4P_5D

   subroutine assign_allocatable_gpu_I8P_2D(lhs, rhs, omp_dev, msg, verbose)
   !< Assign GPU variable with memory checking (kind I8P, rank 2).
   !< Variable is returned not allocated if right hand side is not allocated.
   Integer(I8P), pointer,     intent(inout)         :: lhs(:,:) !< Left hand side of assignement.
   Integer(I8P), allocatable, intent(in), target    :: rhs(:,:) !< Right hand side of assignement.
   integer(I4P),              intent(in)            :: omp_dev  !< Flag to activate verbose mode.
   character(*),              intent(in), optional  :: msg      !< Message to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose  !< Flag to activate verbose mode.
   character(:), allocatable                        :: msg_     !< Message to be printed in verbose mode, local var.
   logical                                          :: verbose_ !< Flag to activate verbose mode, local var.
   integer(I4P)                                     :: error    !< Output of memcpy.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(lhs)) call omp_target_free_f(lhs,omp_dev)
   if (allocated(rhs)) then
      if (size(rhs, dim=1)*size(rhs, dim=2)>0) then
         call alloc_var_gpu(var=lhs,                                                                                     &
                            ulb=reshape([lbound(rhs,dim=1),ubound(rhs,dim=1),lbound(rhs,dim=2),ubound(rhs,dim=2)],[2,2]),&
                            omp_dev=omp_dev, msg=msg_, verbose=verbose_)
         error = omp_target_memcpy_f(fptr_dst=lhs, fptr_src=rhs, dst_off=0_I4P, src_off=0_I4P, &
                                     omp_dst_dev=omp_dev, omp_src_dev=int(omp_get_initial_device(),I4P))
         if (error/=0) then
            print '(A)', 'Error in assign_allocatable_gpu_I8P_2D (omp_target_memcpy_f): '//msg_
            error stop
         endif
      endif
   endif
   endsubroutine assign_allocatable_gpu_I8P_2D

   subroutine assign_allocatable_gpu_I8P_3D(lhs, rhs, omp_dev, msg, verbose)
   !< Assign GPU variable with memory checking (kind I8P, rank 3).
   !< Variable is returned not allocated if right hand side is not allocated.
   Integer(I8P), pointer,     intent(inout)         :: lhs(:,:,:) !< Left hand side of assignement.
   Integer(I8P), allocatable, intent(in), target    :: rhs(:,:,:) !< Right hand side of assignement.
   integer(I4P),              intent(in)            :: omp_dev    !< Flag to activate verbose mode.
   character(*),              intent(in), optional  :: msg        !< Message to be printed in verbose mode.
   logical,                   intent(in), optional  :: verbose    !< Flag to activate verbose mode.
   character(:), allocatable                        :: msg_       !< Message to be printed in verbose mode, local var.
   logical                                          :: verbose_   !< Flag to activate verbose mode, local var.
   integer(I4P)                                     :: error      !< Output of memcpy.

   msg_     = ''      ; if (present(msg    )) msg_     = msg
   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (associated(lhs)) call omp_target_free_f(lhs,omp_dev)
   if (allocated(rhs)) then
      if (size(rhs,dim=1)*size(rhs,dim=2)*size(rhs,dim=3)>0) then
         call alloc_var_gpu(var=lhs,                                                  &
                            ulb=reshape([lbound(rhs,dim=1),ubound(rhs,dim=1),         &
                                         lbound(rhs,dim=2),ubound(rhs,dim=2),         &
                                         lbound(rhs,dim=3),ubound(rhs,dim=3)],[2,3]), &
                            omp_dev=omp_dev, msg=msg_, verbose=verbose_)
         error = omp_target_memcpy_f(fptr_dst=lhs, fptr_src=rhs, dst_off=0_I4P, src_off=0_I4P, &
                                     omp_dst_dev=omp_dev, omp_src_dev=int(omp_get_initial_device(),I4P))
         if (error/=0) then
            print '(A)', 'Error in assign_allocatable_gpu_I8P_3D (omp_target_memcpy_f): '//msg_
            error stop
         endif
      endif
   endif
   endsubroutine assign_allocatable_gpu_I8P_3D

   subroutine save_memory_gpu_status(file_name, tag)
   !< Save the current GPU-memory status into a file.
   !< File is accessed in append position.
   character(*), intent(in)           :: file_name           !< File name.
   character(*), intent(in), optional :: tag                 !< Tag of current status.
   character(:), allocatable          :: tag_                !< Tag of current status, local var.
   integer(I4P)                       :: mem_free, mem_total !< Process memory.
   integer(I4P)                       :: file_unit           !< File unit.
   integer(I4P)                       :: error               !< Error traping flag.

   tag_ = '' ; if (present(tag)) tag_ = trim(tag)
   ! error = cudaMemGetInfo(mem_free, mem_total)
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
endmodule adam_memory_gmp_library
