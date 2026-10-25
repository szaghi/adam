!< ADAM, PRISM field probes, FNL backend.

#include "fundal.H"

module adam_prism_fnl_probes_object
!< ADAM, PRISM field probes, FNL backend.
!<
!< Device twin of the probe stencils of `prism_probes_object`: the probes are sampled on the device and only the
!< `nv*probes_number` local sums are copied to the host, which then reduces and writes them as on the CPU.

! PRISM common classes, libraries, parameters
use :: adam_prism_common_library
! third party modules
use :: fundal
use :: penf

implicit none
private
public :: prism_fnl_probes_object

type :: prism_fnl_probes_object
   integer(I4P)          :: entries = 0_I4P          !< Local stencil cells.
   integer(I4P), pointer :: first_gpu(:)=>null()     !< First local cell of each probe [probes_number+1].
   integer(I4P), pointer :: b_gpu(:)=>null()         !< Block of each local cell [entries].
   integer(I4P), pointer :: i_gpu(:)=>null()         !< i of each local cell [entries].
   integer(I4P), pointer :: j_gpu(:)=>null()         !< j of each local cell [entries].
   integer(I4P), pointer :: k_gpu(:)=>null()         !< k of each local cell [entries].
   real(R8P),    pointer :: weight_gpu(:)=>null()    !< Weight of each local cell [entries].
   real(R8P),    pointer :: value_gpu(:)=>null()     !< Local sums, v fastest [nv*probes_number].
   real(R8P), allocatable :: value(:)                !< Host copy of the local sums [nv*probes_number].
   contains
      ! public methods
      procedure, pass(self) :: destroy    !< Free device data owned by this object.
      procedure, pass(self) :: initialize !< Copy the host stencils to the device.
      procedure, pass(self) :: sample     !< Sample the device fields into the host probes (local sums).
endtype prism_fnl_probes_object

contains
   ! public methods
   subroutine destroy(self)
   !< Free device data owned by this object.
   class(prism_fnl_probes_object), intent(inout) :: self !< Probes, device.

   if (associated(self%first_gpu)) then
      call dev_free(self%first_gpu, mydev)
      nullify(self%first_gpu)
   endif
   if (associated(self%b_gpu)) then
      call dev_free(self%b_gpu, mydev)
      nullify(self%b_gpu)
   endif
   if (associated(self%i_gpu)) then
      call dev_free(self%i_gpu, mydev)
      nullify(self%i_gpu)
   endif
   if (associated(self%j_gpu)) then
      call dev_free(self%j_gpu, mydev)
      nullify(self%j_gpu)
   endif
   if (associated(self%k_gpu)) then
      call dev_free(self%k_gpu, mydev)
      nullify(self%k_gpu)
   endif
   if (associated(self%weight_gpu)) then
      call dev_free(self%weight_gpu, mydev)
      nullify(self%weight_gpu)
   endif
   if (associated(self%value_gpu)) then
      call dev_free(self%value_gpu, mydev)
      nullify(self%value_gpu)
   endif
   if (allocated(self%value)) deallocate(self%value)
   self%entries = 0_I4P
   endsubroutine destroy

   subroutine initialize(self, probes, nv)
   !< Copy the host stencils to the device; a rank without stencil cells keeps no device data.
   class(prism_fnl_probes_object), intent(inout) :: self   !< Probes, device.
   type(prism_probes_object),      intent(in)    :: probes !< Probes, host, stencils built.
   integer(I4P),                   intent(in)    :: nv     !< Number of sampled variables.
   integer(I4P), allocatable                     :: cell_row(:) !< One row of the local cells, contiguous.

   call self%destroy
   if (probes%probes_number == 0_I4P) return
   allocate(self%value(nv*probes%probes_number))
   self%value = 0._R8P
   self%entries = size(probes%weight, kind=I4P)
   ! zero-size device arrays are not allocated (a zero-size memcpy faults, see prism_fnl_coil_object%copy_cpu_gpu)
   if (self%entries == 0_I4P) return
   call dev_assign_to_device(src=probes%first,  dst=self%first_gpu)
   ! the rows of `cell` are strided: copy them through contiguous buffers
   cell_row = probes%cell(1,:) ; call dev_assign_to_device(src=cell_row, dst=self%b_gpu)
   cell_row = probes%cell(2,:) ; call dev_assign_to_device(src=cell_row, dst=self%i_gpu)
   cell_row = probes%cell(3,:) ; call dev_assign_to_device(src=cell_row, dst=self%j_gpu)
   cell_row = probes%cell(4,:) ; call dev_assign_to_device(src=cell_row, dst=self%k_gpu)
   call dev_assign_to_device(src=probes%weight, dst=self%weight_gpu)
   call dev_assign_to_device(src=self%value,    dst=self%value_gpu)
   endsubroutine initialize

   subroutine sample(self, probes, ngc, q_gpu)
   !< Sample the device fields into the host probes (local sums, reduced by `prism_probes_object%save_history`).
   class(prism_fnl_probes_object), intent(inout) :: self   !< Probes, device.
   type(prism_probes_object),      intent(inout) :: probes !< Probes, host.
   integer(I4P),                   intent(in)    :: ngc    !< Number of ghost cells.
   real(R8P),                      intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Fields [nb,ni,nj,nk,nv].

   if (probes%probes_number == 0_I4P) return
   if (self%entries > 0_I4P) then
      call sample_dev_kernel(nv=probes%nv, probes_number=probes%probes_number, ngc=ngc, first_gpu=self%first_gpu, &
                             b_gpu=self%b_gpu, i_gpu=self%i_gpu, j_gpu=self%j_gpu, k_gpu=self%k_gpu,              &
                             weight_gpu=self%weight_gpu, q_gpu=q_gpu, value_gpu=self%value_gpu)
      call dev_memcpy_from_device(src=self%value_gpu, dst=self%value)
   endif
   probes%value = reshape(self%value, [probes%nv, probes%probes_number])
   endsubroutine sample

   ! private non TBP
   subroutine sample_dev_kernel(nv, probes_number, ngc, first_gpu, b_gpu, i_gpu, j_gpu, k_gpu, weight_gpu, q_gpu, &
                                value_gpu)
   !< Local weighted sums of the stencil cells, one thread per (variable, probe).
   integer(I4P), intent(in)    :: nv, probes_number, ngc  !< Sizes.
   integer(I4P), intent(in)    :: first_gpu(1:)           !< First local cell of each probe.
   integer(I4P), intent(in)    :: b_gpu(1:), i_gpu(1:), j_gpu(1:), k_gpu(1:) !< Local cells.
   real(R8P),    intent(in)    :: weight_gpu(1:)          !< Local cell weights.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Fields.
   real(R8P),    intent(inout) :: value_gpu(1:)           !< Local sums, v fastest.
   real(R8P)                   :: total                   !< Local sum.
   integer(I4P)                :: p, v, s                 !< Counters.

   !$acc parallel loop collapse(2) DEVICEVAR(first_gpu,b_gpu,i_gpu,j_gpu,k_gpu,weight_gpu,q_gpu,value_gpu) &
   !$acc& firstprivate(nv,probes_number) private(total)
   !$omp OMPLOOP collapse(2) DEVICEPTR(first_gpu,b_gpu,i_gpu,j_gpu,k_gpu,weight_gpu,q_gpu,value_gpu) &
   !$omp& firstprivate(nv,probes_number) private(total)
   do p=1, probes_number
      do v=1, nv
         total = 0._R8P
         !$acc loop seq
         do s=first_gpu(p), first_gpu(p+1)-1
            total = total + weight_gpu(s) * q_gpu(b_gpu(s),i_gpu(s),j_gpu(s),k_gpu(s),v)
         enddo
         value_gpu(v+(p-1)*nv) = total
      enddo
   enddo
   endsubroutine sample_dev_kernel
endmodule adam_prism_fnl_probes_object
