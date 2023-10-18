!< ADAM, MPI handler class definition, NVF backend.
module adam_mpih_nvf_object
!< ADAM, MPI handler class definition, NVF backend.
!<
!< Extend common mpih class adding NVF features.

use adam_mpih_object
use penf
use cudafor
use mpi

implicit none
private
public ::  mpih_nvf_object

type, extends(mpih_object) :: mpih_nvf_object
   !< MPI handler class, NVF backend.
   integer(I4P) :: mydev=0_I4P      !< My GPU rank.
   integer(I4P) :: local_comm=0_I4P !< Local communicator.
   contains
      ! override methods
      procedure, pass(self) :: initialize !< Initialize MPI handler data.
      ! public methods
      procedure, pass(self) :: check_cuda_error        !< Check if CUDA error occurs and abort in case.
      procedure, nopass     :: compute_cuda_dimensions !< Compute CUDA grid dimensions for GPU parallel computations.
      procedure, pass(self) :: print_device_properties !< Pretty print device properties.
endtype mpih_nvf_object

contains
    ! override methods
   subroutine initialize(self, do_mpi_init, do_device_init)
   !< Initialize MPI handler data.
   class(mpih_nvf_object) , intent(inout)        :: self              !< MPI handler.
   logical,                 intent(in), optional :: do_mpi_init       !< Flag to activate MPI init call.
   logical,                 intent(in), optional :: do_device_init    !< Flag to activate device init call.
   type(cudadeviceprop)                          :: device_properties !< Device properties.

   call self%mpih_object%initialize(do_mpi_init=do_mpi_init)
   call self%print_message('mpih_nvf_object%initialize start')
   if (present(do_device_init)) then
      if (do_device_init) then
         call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, self%local_comm, self%error)
         call MPI_COMM_RANK(self%local_comm, self%mydev, self%error)
         self%error = CudaSetDevice(self%mydev)
         self%error = cudaGetDeviceProperties(device_properties, self%mydev)
         call self%print_device_properties(device_properties)
         self%memory_avail = real(device_properties%totalGlobalMem, R8P)/1e9
      endif
   endif
   call self%print_message('mpih_nvf_object%initialize finish')
   endsubroutine initialize

   ! public methods
   subroutine check_cuda_error(self, error_code, msg)
   !< Check if CUDA error occurs and abort in case.
   class(mpih_nvf_object), intent(inout), optional :: self       !< MPI handler.
   integer(I4P),           intent(in),    optional :: error_code !< Abort error code.
   character(*),           intent(in),    optional :: msg        !< Error message.
   character(:), allocatable                       :: msg_       !< Error message, local variable.

   self%error = cudaGetLastError()
   if (self%error /= cudaSuccess) then
      msg_ = cudaGetErrorString(self%error) ; if (present(msg)) msg_ = '"'//msg//'": '//msg_
      call self%abort(error_code=error_code, msg=msg_)
   endif
   endsubroutine

   subroutine compute_cuda_dimensions(grid_x, grid_y, grid, tBlock)
   !< Compute CUDA grid dimensions for GPU parallel computations.
   integer(I4P), intent(in)  :: grid_x, grid_y !< CUDA grid xy dimensions.
   type(dim3), intent(inout) :: grid, tBlock   !< CUDA grid and block.

   tBlock = dim3(32,8,1)
   grid = dim3(ceiling(real(grid_x)/tBlock%x), ceiling(real(grid_y)/tBlock%y), 1)
   endsubroutine compute_cuda_dimensions

   subroutine print_device_properties(self, device_properties)
   !< Pretty print device properties.
   class(mpih_nvf_object), intent(in) :: self              !< MPI handler.
   type(cudadeviceprop),   intent(in) :: device_properties !< Device properties.

   associate(r=>self%myrankstr)
   print'(A)',r//"total global memory:         "//trim(str(real(device_properties%totalGlobalMem)/1e9           ,.true.))//" Gbytes"
   print'(A)',r//"shared mem per block:        "//trim(str(     device_properties%sharedMemPerBlock             ,.true.))//" bytes"
   print'(A)',r//"regs per block:              "//trim(str(     device_properties%regsPerBlock                  ,.true.))
   print'(A)',r//"warp size:                   "//trim(str(     device_properties%warpSize                      ,.true.))
   print'(A)',r//"max threads per block:       "//trim(str(     device_properties%maxThreadsPerBlock            ,.true.))
   print'(A)',r//"max threads dim:             "//trim(str(     device_properties%maxThreadsDim                 ,.true.))
   print'(A)',r//"clock rate:                  "//trim(str(real(device_properties%clockRate)/1e6                ,.true.))//" GHz"
   print'(A)',r//"total const memory:          "//trim(str(     device_properties%totalConstMem                 ,.true.))//" bytes"
   print'(A)',r//"compute capability revision: "//trim(str(    [device_properties%major,device_properties%minor],.true.))
   print'(A)',r//"multi processor count:       "//trim(str(     device_properties%multiProcessorCount           ,.true.))
   print'(A)',r//"L2 cache size:               "//trim(str(     device_properties%l2CacheSize                   ,.true.))
   print'(A)',r//"max threads per SMP:         "//trim(str(     device_properties%maxThreadsPerMultiProcessor   ,.true.))
   print'(A)',r//"device rank:                 "//trim(str(     self%mydev                                      ,.true.))
   endassociate
   endsubroutine print_device_properties
endmodule adam_mpih_nvf_object
