!< ADAM, PRISM coil source definition, FNL backend.

#include "fundal.H"

module adam_prism_fnl_coil_object
!< ADAM, PRISM coil source definition, FNL backend.

! ADAM classes, libraries, parameters
use :: adam_common_library
! ADAM FNL classes, libraries, parameters
use :: adam_fnl_library
! PRISM common classes, libraries, parameters
use :: adam_prism_common_library
! third party modules
use :: fundal
use :: penf

implicit none
private
public :: prism_fnl_coil_object

type :: prism_fnl_coil_object
   real(R8P),    pointer :: A_gpu(:)=>null()               !< Current amplitude (A)
   real(R8P),    pointer :: f_gpu(:)=>null()               !< Current frequency, if AC (Hz)
   real(R8P),    pointer :: phase_gpu(:)=>null()           !< Current initial phase, if AC
   real(R8P),    pointer :: J_vec_gpu(:,:,:,:,:,:)=>null() !< Matrice contenente versori corrente spire (se assente = 0)
   ! Plasma current on the coil support, device twin of prism_coil_object%(support, plasma_j): the same cell list (pushed
   ! with J_vec), the plasma current kept per slot and summed by the coil stamp (register_plasma_dev, restore_plasma_dev).
   integer(I4P), pointer :: support_gpu(:,:)=>null()       !< Cells (b,i,j,k) with J_vec /= 0, interior first, (4,n).
   real(R8P),    pointer :: plasma_j_gpu(:,:,:)=>null()    !< Plasma current on the support, (3,n,slot).
   integer(I4P)          :: support_number=0_I4P           !< Number of support cells.
   integer(I4P)          :: support_interior=0_I4P         !< Number of interior support cells (the first ones).
   integer(I4P), allocatable :: support_host(:,:)          !< Host copy of the device support (to detect a change).
   logical               :: plasma_j_active(2)=.false.     !< Slot holds a plasma current.
   contains
      ! public methods
      procedure, pass(self) :: copy_cpu_gpu      !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu      !< Copy data from GPU to CPU.
      procedure, pass(self) :: destroy           !< Free device data owned by this object.
      procedure, pass(self) :: initialize        !< Initialize class from global singletons.
      procedure, pass(self) :: capture_plasma_dev  !< Keep the plasma current of the coil support.
      procedure, pass(self) :: detect_plasma_dev   !< Plasma on the interior coil support.
      procedure, pass(self) :: restore_plasma_dev  !< Write the kept plasma current back (inside the coil stamp).
      ! private methods
      procedure, pass(self), private :: free_support !< Free the device support and its plasma slots.
      procedure, pass(self), private :: push_support !< Push the support of the host J_vec and its committed slot.
endtype prism_fnl_coil_object

contains
   ! public methods
   subroutine destroy(self)
   !< Free device data owned by the coil FNL object.
   class(prism_fnl_coil_object), intent(inout) :: self !< Coils.

   if (associated(self%A_gpu)) then
      call dev_free(self%A_gpu, mydev)
      nullify(self%A_gpu)
   endif
   if (associated(self%f_gpu)) then
      call dev_free(self%f_gpu, mydev)
      nullify(self%f_gpu)
   endif
   if (associated(self%phase_gpu)) then
      call dev_free(self%phase_gpu, mydev)
      nullify(self%phase_gpu)
   endif
   if (associated(self%j_vec_gpu)) then
      call dev_free(self%j_vec_gpu, mydev)
      nullify(self%j_vec_gpu)
   endif
   call self%free_support
   endsubroutine destroy

   subroutine copy_cpu_gpu(self, coil, grid, blocks_number, buf6D, db6, hb6, verbose)
   !< Copy data from CPU to GPU.
   class(prism_fnl_coil_object), intent(inout)           :: self         !< The field.
   class(prism_coil_object),     intent(in)              :: coil         !< Coils on host.
   type(grid_object),            intent(in)              :: grid         !< Grid (sibling realm component, threaded in).
   integer(I4P),                 intent(in)              :: blocks_number !< Number of blocks in use.
   real(R8P),                    intent(inout), optional :: buf6D(1:,                                 &
                                                                  1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,&
                                                                  1:,1:) !< Buffer (host memory, device shape), rank 6D.
   integer(I4P),                 intent(in),    optional :: db6(2,6)     !< Device data bounds, rank 6.
   integer(I4P),                 intent(in),    optional :: hb6(2,6)     !< Host data bounds, rank 6.
   logical,                      intent(in),    optional :: verbose      !< Flag to activate verbose mode.
   logical                                               :: verbose_     !< Flag to activate verbose mode, local var.
   integer(I4P)                                          :: db6_(2,6)    !< Device data bounds, rank 6.
   integer(I4P)                                          :: hb6_(2,6)    !< Host data bounds, rank 6.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   ! Zero-coil case (issue #22 F1): the host coil arrays are unallocated and the device
   ! twins are zero-size — any memcpy would fault (SIGSEGV in cuMemcpy on the first
   ! source-free case ever run on FNL, rmf-amr-fd-pulse). Nothing to move: return.
   if (coil%total_coils_number == 0) return
   if (verbose_) call mpih_fnl%print_message('prism_fnl_coil_object%copy_cpu_gpu start')
   call dev_memcpy_to_device(src=coil%coil_amplitude    ,dst=self%A_gpu    )
   call dev_memcpy_to_device(src=coil%f                 ,dst=self%f_gpu    )
   call dev_memcpy_to_device(src=coil%phase             ,dst=self%phase_gpu)
   if (present(buf6D)) then
      if (present(db6) .and. present(hb6)) then
         call dev_memcpy_to_device(bb=db6,ij=[1,5],tb=hb6,dst=self%j_vec_gpu,src=coil%j_vec,buf=buf6D)
      else
         db6_(1,:) = lbound(self%j_vec_gpu ) ; db6_(2,:) = ubound(self%j_vec_gpu )
         hb6_(1,:) = lbound(coil%j_vec) ; hb6_(2,:) = ubound(coil%j_vec)
         call dev_memcpy_to_device(bb=db6_,ij=[1,5],tb=hb6_,dst=self%j_vec_gpu,src=coil%j_vec,buf=buf6D)
      endif
   else
      call dev_assign_to_device(src=coil%j_vec,dst=self%j_vec_gpu,ij=[1,5])
   endif
   call self%push_support(coil=coil, blocks_number=blocks_number)
   if (verbose_) call mpih_fnl%print_message('prism_fnl_coil_object%copy_cpu_gpu finish')
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self, coil, grid, buf6D, db6, hb6, verbose)
   !< Copy data from GPU to CPU.
   class(prism_fnl_coil_object), intent(inout)           :: self         !< The field.
   class(prism_coil_object),     intent(inout)           :: coil         !< Coils on host.
   type(grid_object),            intent(in)              :: grid         !< Grid (sibling realm component, threaded in).
   real(R8P),                    intent(inout), optional :: buf6D(1:,                                 &
                                                                  1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,&
                                                                  1:,1:) !< Buffer (host memory, device shape), rank 6D.
   integer(I4P),                 intent(in),    optional :: db6(2,6)     !< Device data bounds, rank 6.
   integer(I4P),                 intent(in),    optional :: hb6(2,6)     !< Host data bounds, rank 6.
   logical,                      intent(in),    optional :: verbose      !< Flag to activate verbose mode.
   logical                                               :: verbose_     !< Flag to activate verbose mode, local var.
   integer(I4P)                                          :: db6_(2,6)    !< Device data bounds, rank 6.
   integer(I4P)                                          :: hb6_(2,6)    !< Host data bounds, rank 6.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   ! Zero-coil case (issue #22 F1): see copy_cpu_gpu — nothing to move, and the
   ! host destinations are unallocated. Return before any device memcpy.
   if (coil%total_coils_number == 0) return
   if (verbose_) call mpih_fnl%print_message('prism_fnl_coil_object%copy_gpu_cpu start')
   call dev_memcpy_from_device(src=self%A_gpu    ,dst=coil%coil_amplitude)
   call dev_memcpy_from_device(src=self%f_gpu    ,dst=coil%f             )
   call dev_memcpy_from_device(src=self%phase_gpu,dst=coil%phase         )
   if (present(buf6D)) then
      if (present(db6) .and. present(hb6)) then
         call dev_memcpy_from_device(bb=db6,ij=[1,5],tb=hb6,dst=coil%j_vec,src=self%j_vec_gpu,buf=buf6D)
      else
         db6_(1,:) = lbound(self%j_vec_gpu) ; db6_(2,:) = ubound(self%j_vec_gpu )
         hb6_(1,:) = lbound(coil%j_vec    ) ; hb6_(2,:) = ubound(coil%j_vec     )
         call dev_memcpy_from_device(bb=db6_,ij=[1,5],tb=hb6_,dst=coil%j_vec,src=self%j_vec_gpu,buf=buf6D)
      endif
   else
      call dev_assign_from_device(src=self%j_vec_gpu,dst=coil%j_vec,ij=[1,5])
   endif
   if (verbose_) call mpih_fnl%print_message('prism_fnl_coil_object%copy_gpu_cpu finish')
   endsubroutine copy_gpu_cpu

   subroutine initialize(self, coil, field, grid, buf6D, db6, hb6)
   !< Initialize class from program-scope `field` (adam_field_global) and `grid` (adam_grid_global) singletons.
   !< Requires `mpih_fnl` (adam_fnl_mpih_global), `field` and `grid` singletons to be ready.
   class(prism_fnl_coil_object), intent(inout)      :: self !< Coils.
   class(prism_coil_object),     intent(in), target :: coil !< Coils on host.
   type(field_object),           intent(in)         :: field !< Field (sibling realm component, threaded in).
   type(grid_object),            intent(in)         :: grid !< Grid (sibling realm component, threaded in).
   real(R8P),                    intent(inout), optional :: buf6D(1:,                                 &
                                                                  1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,&
                                                                  1:,1:) !< Buffer (host memory, device shape), rank 6D.
   integer(I4P),                 intent(in),    optional :: db6(2,6) !< Device data bounds, rank 6.
   integer(I4P),                 intent(in),    optional :: hb6(2,6) !< Host data bounds, rank 6.
   integer(I4P)                                     :: nv             !< Counter.
   integer(I4P)                                     :: ierr           !< Error status.
   integer(I8P)                                     :: j_vec_bytes    !< Requested j_vec_gpu bytes.
   integer(I8P)                                     :: j_vec_elements !< Requested j_vec_gpu elements.
   real(R8P)                                        :: j_vec_gib      !< Requested j_vec_gpu GiB.
   real(R8P)                                        :: dev_free_gib   !< Device memory reported at device init, GiB.
   character(len=:), allocatable                    :: alloc_msg      !< Allocation diagnostic.

   call self%destroy()
   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, ngc=>grid%ngc, nb=>field%nb, blocks_number=>field%blocks_number, &
             nc=>coil%total_coils_number)
   print '(A)', mpih_fnl%myrankstr//'prism_fnl_coil_object%initialize start'
   ! Zero-coil case (issue #22 F1): host j_vec is unallocated (size() on it is illegal)
   ! and every device twin would be zero-size. Leave the device pointers null; all the
   ! copy paths return early on total_coils_number == 0.
   if (nc == 0) then
      print '(A)', mpih_fnl%myrankstr//'prism_fnl_coil_object%initialize finish (no coils: nothing on device)'
      return
   endif
   nv = size(coil%j_vec,dim=1)
   j_vec_elements = int(nb, I8P) * int(ni + 2_I4P * ngc, I8P) * int(nj + 2_I4P * ngc, I8P) * &
                    int(nk + 2_I4P * ngc, I8P) * int(nv, I8P) * int(nc, I8P)
   j_vec_bytes = j_vec_elements * int(storage_size(0._R8P) / 8, I8P)
   j_vec_gib = real(j_vec_bytes, R8P) / 1024._R8P**3
   dev_free_gib = real(mpih_fnl%dev_memory_avail, R8P) / 1024._R8P**3
   call dev_alloc(fptr_dev=self%A_gpu        ,ubounds=[nc                           ],lbounds=[0                      ],ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate A_gpu in prism_fnl_coil_object%initialize')
   call dev_alloc(fptr_dev=self%f_gpu        ,ubounds=[nc                           ],lbounds=[0                      ],ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate f_gpu in prism_fnl_coil_object%initialize')
   call dev_alloc(fptr_dev=self%phase_gpu    ,ubounds=[nc                           ],lbounds=[0                      ],ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate phase_gpu in prism_fnl_coil_object%initialize')
   call dev_alloc(fptr_dev=self%j_vec_gpu    ,ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv,nc],lbounds=[1,1-ngc,1-ngc,1-ngc,1,1],ierr=ierr)
   if (ierr /= 0_I4P) then
      alloc_msg = ': failed to allocate j_vec_gpu in prism_fnl_coil_object%initialize'// &
                  ' ierr='//trim(str(ierr))// &
                  ' requested_GiB='//trim(str(j_vec_gib))// &
                  ' device_init_free_GiB='//trim(str(dev_free_gib))// &
                  ' device='//trim(str(mpih_fnl%mydev))// &
                  ' nb='//trim(str(nb))// &
                  ' blocks_number='//trim(str(blocks_number))// &
                  ' ni='//trim(str(ni))//' nj='//trim(str(nj))//' nk='//trim(str(nk))// &
                  ' ngc='//trim(str(ngc))//' nv_coil='//trim(str(nv))//' coils='//trim(str(nc))
      call mpih_fnl%error_stop(msg=alloc_msg)
   endif
   if (present(buf6D) .and. present(db6) .and. present(hb6)) then
      call self%copy_cpu_gpu(coil=coil, grid=grid, blocks_number=blocks_number, buf6D=buf6D, db6=db6, hb6=hb6)
   elseif (present(buf6D)) then
      call self%copy_cpu_gpu(coil=coil, grid=grid, blocks_number=blocks_number, buf6D=buf6D)
   else
      call self%copy_cpu_gpu(coil=coil, grid=grid, blocks_number=blocks_number)
   endif
   print '(A)', mpih_fnl%myrankstr//'prism_fnl_coil_object%initialize finish'
   endassociate
   endsubroutine initialize

   subroutine detect_plasma_dev(self, ngc, q_gpu, var_jx, var_jy, var_jz, check_current, check_charge, overlap_current, &
                                overlap_charge)
   !< Device twin of prism_coil_object%detect_plasma: plasma current (and charge) on the interior support cells right after
   !< a PIC deposit into q_gpu. Local to this rank.
   class(prism_fnl_coil_object), intent(in)  :: self            !< Coils, device data.
   integer(I4P),                 intent(in)  :: ngc             !< Number of ghost cells.
   real(R8P),                    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Deposited field (device).
   integer(I4P),                 intent(in)  :: var_jx          !< Index of Jx.
   integer(I4P),                 intent(in)  :: var_jy          !< Index of Jy.
   integer(I4P),                 intent(in)  :: var_jz          !< Index of Jz.
   logical,                      intent(in)  :: check_current   !< The J rows hold a fresh deposit.
   logical,                      intent(in)  :: check_charge    !< The charge row holds a fresh deposit.
   logical,                      intent(out) :: overlap_current !< Plasma current on the support.
   logical,                      intent(out) :: overlap_charge  !< Plasma charge on the support.
   integer(I4P)                              :: oc, oq          !< Reduction flags.
   integer(I4P)                              :: nv_q            !< Width of q (charge is the last row).

   overlap_current = .false.
   overlap_charge  = .false.
   if (self%support_interior <= 0_I4P) return
   nv_q = int(size(q_gpu, dim=5), I4P)
   call detect_kernel(n=self%support_interior, ngc=ngc, nv_q=nv_q, var_jx=var_jx, var_jy=var_jy, var_jz=var_jz,     &
                      do_current=merge(1_I4P, 0_I4P, check_current), do_charge=merge(1_I4P, 0_I4P, check_charge), &
                      support_gpu=self%support_gpu, q_gpu=q_gpu, overlap_current=oc, overlap_charge=oq)
   overlap_current = oc /= 0_I4P
   overlap_charge  = oq /= 0_I4P
   contains
      subroutine detect_kernel(n, ngc, nv_q, var_jx, var_jy, var_jz, do_current, do_charge, support_gpu, q_gpu, &
                               overlap_current, overlap_charge)
      !< Plasma current or charge on the first n support cells, device reduction.
      integer(I4P), intent(in)  :: n, ngc, nv_q, var_jx, var_jy, var_jz, do_current, do_charge !< Sizes and switches.
      integer(I4P), intent(in)  :: support_gpu(1:,1:)                                           !< Support cells.
      real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)                            !< Deposited field.
      integer(I4P), intent(out) :: overlap_current, overlap_charge                              !< Reduction flags.
      integer(I4P)              :: c, b, i, j, k                                                !< Counters.

      overlap_current = 0_I4P
      overlap_charge  = 0_I4P
      !$acc parallel loop independent gang vector DEVICEVAR(support_gpu, q_gpu) &
      !$acc& firstprivate(n, nv_q, var_jx, var_jy, var_jz, do_current, do_charge) &
      !$acc& reduction(max: overlap_current, overlap_charge)
      !$omp OMPLOOP DEVICEPTR(support_gpu, q_gpu) &
      !$omp& firstprivate(n, nv_q, var_jx, var_jy, var_jz, do_current, do_charge) &
      !$omp& reduction(max: overlap_current, overlap_charge)
      do c=1, n
         b = support_gpu(1,c) ; i = support_gpu(2,c) ; j = support_gpu(3,c) ; k = support_gpu(4,c)
         if (do_current /= 0_I4P) then
            if ((q_gpu(b,i,j,k,var_jx) /= 0._R8P) .or. (q_gpu(b,i,j,k,var_jy) /= 0._R8P) .or. &
                (q_gpu(b,i,j,k,var_jz) /= 0._R8P)) overlap_current = 1_I4P
         endif
         if (do_charge /= 0_I4P) then
            if (q_gpu(b,i,j,k,nv_q) /= 0._R8P) overlap_charge = 1_I4P
         endif
      enddo
      endsubroutine detect_kernel
   endsubroutine detect_plasma_dev

   subroutine capture_plasma_dev(self, ngc, q_gpu, slot, var_jx, var_jy, var_jz, active)
   !< Device twin of prism_coil_object%capture_plasma: keep the J rows of the support cells (interior and ghost, the ghosts
   !< refreshed by the caller) in `slot`, or empty the slot (no plasma current on the support: the stamp is unchanged).
   class(prism_fnl_coil_object), intent(inout) :: self   !< Coils, device data.
   integer(I4P),                 intent(in)    :: ngc    !< Number of ghost cells.
   real(R8P),                    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Deposited field (device).
   integer(I4P),                 intent(in)    :: slot   !< COIL_SLOT_COMMITTED or COIL_SLOT_STAGE.
   integer(I4P),                 intent(in)    :: var_jx !< Index of Jx.
   integer(I4P),                 intent(in)    :: var_jy !< Index of Jy.
   integer(I4P),                 intent(in)    :: var_jz !< Index of Jz.
   logical,                      intent(in)    :: active !< Plasma current on the support (of some rank).

   if (self%support_number <= 0_I4P) return
   self%plasma_j_active(slot) = active
   if (.not.active) return
   call capture_kernel(n=self%support_number, ngc=ngc, slot=slot, var_jx=var_jx, var_jy=var_jy, var_jz=var_jz, &
                       support_gpu=self%support_gpu, q_gpu=q_gpu, plasma_j_gpu=self%plasma_j_gpu)
   contains
      subroutine capture_kernel(n, ngc, slot, var_jx, var_jy, var_jz, support_gpu, q_gpu, plasma_j_gpu)
      !< Keep the J rows of the support cells in the slot, device kernel.
      integer(I4P), intent(in)    :: n, ngc, slot, var_jx, var_jy, var_jz !< Sizes and indexes.
      integer(I4P), intent(in)    :: support_gpu(1:,1:)                   !< Support cells.
      real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Deposited field.
      real(R8P),    intent(inout) :: plasma_j_gpu(1:,1:,1:)               !< Plasma current slots.
      integer(I4P)                :: c, b, i, j, k                        !< Counters.

      !$acc parallel loop independent gang vector DEVICEVAR(support_gpu, q_gpu, plasma_j_gpu) &
      !$acc& firstprivate(n, slot, var_jx, var_jy, var_jz)
      !$omp OMPLOOP DEVICEPTR(support_gpu, q_gpu, plasma_j_gpu) firstprivate(n, slot, var_jx, var_jy, var_jz)
      do c=1, n
         b = support_gpu(1,c) ; i = support_gpu(2,c) ; j = support_gpu(3,c) ; k = support_gpu(4,c)
         plasma_j_gpu(1,c,slot) = q_gpu(b,i,j,k,var_jx)
         plasma_j_gpu(2,c,slot) = q_gpu(b,i,j,k,var_jy)
         plasma_j_gpu(3,c,slot) = q_gpu(b,i,j,k,var_jz)
      enddo
      endsubroutine capture_kernel
   endsubroutine capture_plasma_dev

   subroutine restore_plasma_dev(self, ngc, q_gpu, slot, var_jx, var_jy, var_jz)
   !< Device twin of prism_coil_object%restore_plasma: inside the coil stamp, after the support was cleared, write back the
   !< plasma current kept in `slot` (nothing when the slot is empty); the coil current is then added on top.
   class(prism_fnl_coil_object), intent(in)    :: self   !< Coils, device data.
   integer(I4P),                 intent(in)    :: ngc    !< Number of ghost cells.
   real(R8P),                    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field (device).
   integer(I4P),                 intent(in)    :: slot   !< COIL_SLOT_COMMITTED or COIL_SLOT_STAGE.
   integer(I4P),                 intent(in)    :: var_jx !< Index of Jx.
   integer(I4P),                 intent(in)    :: var_jy !< Index of Jy.
   integer(I4P),                 intent(in)    :: var_jz !< Index of Jz.

   if (self%support_number <= 0_I4P) return
   if (.not.self%plasma_j_active(slot)) return
   call restore_kernel(n=self%support_number, ngc=ngc, slot=slot, var_jx=var_jx, var_jy=var_jy, var_jz=var_jz, &
                       support_gpu=self%support_gpu, plasma_j_gpu=self%plasma_j_gpu, q_gpu=q_gpu)
   contains
      subroutine restore_kernel(n, ngc, slot, var_jx, var_jy, var_jz, support_gpu, plasma_j_gpu, q_gpu)
      !< Write the slot back into the J rows of the support cells, device kernel.
      integer(I4P), intent(in)    :: n, ngc, slot, var_jx, var_jy, var_jz !< Sizes and indexes.
      integer(I4P), intent(in)    :: support_gpu(1:,1:)                   !< Support cells.
      real(R8P),    intent(in)    :: plasma_j_gpu(1:,1:,1:)               !< Plasma current slots.
      real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Field.
      integer(I4P)                :: c, b, i, j, k                        !< Counters.

      !$acc parallel loop independent gang vector DEVICEVAR(support_gpu, plasma_j_gpu, q_gpu) &
      !$acc& firstprivate(n, slot, var_jx, var_jy, var_jz)
      !$omp OMPLOOP DEVICEPTR(support_gpu, plasma_j_gpu, q_gpu) firstprivate(n, slot, var_jx, var_jy, var_jz)
      do c=1, n
         b = support_gpu(1,c) ; i = support_gpu(2,c) ; j = support_gpu(3,c) ; k = support_gpu(4,c)
         q_gpu(b,i,j,k,var_jx) = plasma_j_gpu(1,c,slot)
         q_gpu(b,i,j,k,var_jy) = plasma_j_gpu(2,c,slot)
         q_gpu(b,i,j,k,var_jz) = plasma_j_gpu(3,c,slot)
      enddo
      endsubroutine restore_kernel
   endsubroutine restore_plasma_dev

   ! private methods
   subroutine free_support(self)
   !< Free the device support and its plasma slots.
   class(prism_fnl_coil_object), intent(inout) :: self !< Coils, device data.

   if (associated(self%support_gpu)) then
      call dev_free(self%support_gpu, mydev)
      nullify(self%support_gpu)
   endif
   if (associated(self%plasma_j_gpu)) then
      call dev_free(self%plasma_j_gpu, mydev)
      nullify(self%plasma_j_gpu)
   endif
   if (allocated(self%support_host)) deallocate(self%support_host)
   self%support_number = 0_I4P
   self%support_interior = 0_I4P
   self%plasma_j_active = .false.
   endsubroutine free_support

   subroutine push_support(self, coil, blocks_number)
   !< Push the support of the host J_vec (the list prism_coil_object%support_cells builds, in the host order). An unchanged
   !< support keeps the device slots (a re-push of the same J_vec); a changed one (AMR) starts with empty slots, filled
   !< again by the next deposit. The committed slot of the host is uploaded when it was kept on this same support (the t0
   !< deposit, registered on the host by compute_coils_current_time_zero): the caller then invalidates the host support,
   !< so a later push cannot upload it again.
   class(prism_fnl_coil_object), intent(inout) :: self            !< Coils, device data.
   class(prism_coil_object),     intent(in)    :: coil            !< Coils on host.
   integer(I4P),                 intent(in)    :: blocks_number   !< Number of blocks in use.
   integer(I4P), allocatable                   :: cells(:,:)      !< Support cells (4,n).
   real(R8P),    allocatable                   :: plasma_j(:,:,:) !< Plasma slots, host staging.
   logical                                     :: same_device     !< The device support is this one already.
   logical                                     :: same_host       !< The host slots refer to this support.
   integer(I4P)                                :: n               !< Number of support cells.
   integer(I4P)                                :: interior        !< Number of interior support cells.

   call coil%support_cells(blocks_number=blocks_number, cells=cells, interior=interior)
   n = int(size(cells, dim=2), I4P)
   same_device = .false.
   if (allocated(self%support_host)) then
      if (size(self%support_host, dim=2) == n) same_device = all(self%support_host == cells)
   endif
   same_host = .false.
   if (coil%support_valid .and. allocated(coil%support)) then
      if (size(coil%support, dim=2) == n) same_host = all(coil%support == cells)
   endif
   same_host = same_host .and. coil%plasma_j_active(COIL_SLOT_COMMITTED)
   if (same_device .and. (.not.same_host)) return
   call self%free_support
   self%support_host = cells
   self%support_number = n
   self%support_interior = interior
   if (n <= 0_I4P) return
   allocate(plasma_j(1:3, 1:n, 1:2))
   plasma_j = 0._R8P
   if (same_host) then
      plasma_j(:,:,COIL_SLOT_COMMITTED) = coil%plasma_j(:,:,COIL_SLOT_COMMITTED)
      self%plasma_j_active(COIL_SLOT_COMMITTED) = .true.
   endif
   call dev_assign_to_device(src=cells, dst=self%support_gpu)
   call dev_assign_to_device(src=plasma_j, dst=self%plasma_j_gpu)
   endsubroutine push_support
endmodule adam_prism_fnl_coil_object
