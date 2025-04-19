!< ADAM, MPI handler class definition, NVF backend.
module adam_mpih_nvf_object
!< ADAM, MPI handler class definition, NVF backend.
!<
!< Extend common mpih class adding NVF features.

use adam_mpih_object
use finer
use penf
use cudafor
use mpi

implicit none
private
public :: mpih_nvf_object
public :: CUDA_SILENT_ERR

integer(I4P), parameter :: CUDA_SILENT_ERR = -101_I4P !< CUDA silent error code.

character(len=4), parameter :: INI_SECTION_NAME="cuda" !< INI (config) file section name containing configs.

type, extends(mpih_object) :: mpih_nvf_object
   !< MPI handler class, NVF backend.
   integer(I4P) :: mydev=0_I4P      !< My GPU rank.
   integer(I4P) :: local_comm=0_I4P !< Local communicator.
   integer(I4P) :: iercuda=0_I4P    !< Error trapping flag for CUDAFortran.
   integer(I4P) :: cblk(3)=[32,8,1] !< Default CUDA block dimensions.
   type(dim3)   :: grid, tBlock     !< CUDA grid and block.
   contains
      ! override methods
      procedure, pass(self) :: initialize !< Initialize MPI handler data.
      ! public methods
      procedure, pass(self) :: check_cuda_error        !< Check if CUDA error occurs and abort in case.
      procedure, pass(self) :: description             !< Return pretty-printed object description.
      procedure, pass(self) :: load_from_file          !< Load config from file.
      procedure, pass(self) :: print_device_properties !< Pretty print device properties.
      procedure, pass(self) :: set_cuda_dimensions     !< Compute CUDA grid dimensions for GPU parallel computations.
endtype mpih_nvf_object

contains
    ! override methods
   subroutine initialize(self, do_mpi_init, do_device_init, myrankstr_char_length, verbose)
   !< Initialize MPI handler data.
   class(mpih_nvf_object) , intent(inout)        :: self                  !< MPI handler.
   logical,                 intent(in), optional :: do_mpi_init           !< Flag to activate MPI init call.
   logical,                 intent(in), optional :: do_device_init        !< Flag to activate device init call.
   integer(I4P),            intent(in), optional :: myrankstr_char_length !< MPI ID string length.
   logical,                 intent(in), optional :: verbose               !< Trigger verbose output.
   logical                                       :: verbose_              !< Trigger verbose output, local variable.
   type(cudadeviceprop)                          :: device_properties     !< Device properties.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose

   call self%mpih_object%initialize(do_mpi_init=do_mpi_init, myrankstr_char_length=myrankstr_char_length, verbose=verbose)
   if (verbose_) call self%print_message('mpih_nvf_object%initialize start')
   if (present(do_device_init)) then
      if (do_device_init) then
         call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, self%local_comm, self%error)
         call MPI_COMM_RANK(self%local_comm, self%mydev, self%error)
         self%error = CudaSetDevice(self%mydev)
         self%error = cudaGetDeviceProperties(device_properties, self%mydev)
         if (verbose_) call self%print_device_properties(device_properties)
         self%memory_avail = real(device_properties%totalGlobalMem, R8P)/1e9
      endif
   endif
   if (verbose_) call self%print_message('mpih_nvf_object%initialize finish')
   endsubroutine initialize

   ! public methods
   subroutine check_cuda_error(self, error_code, msg)
   !< Check if CUDA error occurs and abort in case.
   class(mpih_nvf_object), intent(inout)           :: self       !< MPI handler.
   integer(I4P),           intent(in),    optional :: error_code !< Abort error code.
   character(*),           intent(in),    optional :: msg        !< Error message.
   character(:), allocatable                       :: msg_       !< Error message, local variable.

   self%iercuda=cudaDeviceSynchronize()
   self%error = cudaGetLastError()
   if (self%error /= cudaSuccess) then
      msg_ = cudaGetErrorString(self%error) ; if (present(msg)) msg_ = '"'//msg//'": '//msg_
      call self%abort(error_code=error_code, msg=msg_)
   endif
   endsubroutine

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(mpih_nvf_object) , intent(in) :: self             !< MPI handler.
   character(len=:), allocatable       :: desc             !< Description.
   character(len=1), parameter         :: NL=new_line('a') !< New line character.

   desc = self%mpih_object%description()//NL
   desc = desc//self%myrankstr//'MPIH NVF main data'//NL
   desc = desc//self%myrankstr//'  mydev:      '//trim(str(self%mydev     ))//NL
   desc = desc//self%myrankstr//'  local_comm: '//trim(str(self%local_comm))//NL
   desc = desc//self%myrankstr//'  cblk:       '//trim(str(self%cblk      ))
   endfunction description

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(mpih_nvf_object) , intent(inout)        :: self              !< MPI handler.
   type(file_ini),          intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                 intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                       :: go_on_fail_     !< Go on if load fails.
   integer(I4P)                                  :: buffer          !< Buffer.
   integer(I4P)                                  :: err             !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='block_x', val=buffer, error=err)
   if (.not.go_on_fail_.and.err>0) call self%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(block_x)')
   self%cblk(1) = buffer
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='block_y', val=buffer, error=err)
   if (.not.go_on_fail_.and.err>0) call self%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(block_x)')
   self%cblk(2) = buffer
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='block_z', val=buffer, error=err)
   if (.not.go_on_fail_.and.err>0) call self%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(block_x)')
   self%cblk(3) = buffer
   endsubroutine load_from_file

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

   subroutine set_cuda_dimensions(self, cgrd, cblk)
   !< Compute CUDA grid dimensions for GPU parallel computations.
   class(mpih_nvf_object), intent(inout)        :: self    !< MPI handler.
   integer(I4P),           intent(in)           :: cgrd(3) !< CUDA grid dimensions.
   integer(I4P),           intent(in), optional :: cblk(3) !< CUDA block dimensions.

   if (present(cblk)) then
      self%tBlock = dim3(cblk(1),cblk(2),cblk(3))
   else
      self%tBlock = dim3(self%cblk(1),self%cblk(2),self%cblk(3))
   endif
   self%grid = dim3(ceiling(real(cgrd(1))/self%tBlock%x),ceiling(real(cgrd(2))/self%tBlock%y),ceiling(real(cgrd(3))/self%tBlock%z))
   endsubroutine set_cuda_dimensions
endmodule adam_mpih_nvf_object
