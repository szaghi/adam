!< ADAM, PRISM (Plasma Research usIng Simulation Methods) equations system class definition, GPU (FNL) backend.

#include "fundal.H"

module adam_prism_fnl_object
!< ADAM, PRISM (Plasma Research usIng Simulation Methods) equations system class definition, GPU (FNL) backend.

! use adam_prism_fnl_library
! use :: fundal, save_memory_status_gpu=>save_memory_status
! use :: penf, save_memory_status_cpu=>save_memory_status
! use mpi

implicit none
private
public :: prism_fnl_object
public :: NV_PIC, NV_MHD

! integer(I4P), parameter :: NV_PIC = 10_I4P !< Number of PIC variables.
! integer(I4P), parameter :: NV_MHD =  9_I4P !< Number of MHD variables.

type, extends(prism_common_object) :: prism_fnl_object
   !< PRISM equations system class definition, GPU (FNL) backend.
   ! ADAM library objects
   type(mpih_fnl_object)  :: mpih_gpu  !< MPI handler, FNL backend.
   type(field_fnl_object) :: field_gpu !< The field, FNL backend.
   ! device data
   real(R8P), pointer :: q_pic_gpu(:,:,:) !< PIC centered variables.
   !< PIC center variables definition:
   !< q_pic_gpu(b,p,v) where "b" is the block index, "p" is the particle index and
   !< "v" the variable index [blocks_number,particles_number,nv_pic]:
   !< q_pic_gpu(b,p,1) : x coordinate
   !< q_pic_gpu(b,p,2) : y coordinate
   !< q_pic_gpu(b,p,3) : z coordinate
   !< q_pic_gpu(b,p,4) : u velocity component (velocity along x axis)
   !< q_pic_gpu(b,p,5) : v velocity component (velocity along y axis)
   !< q_pic_gpu(b,p,6) : w velocity component (velocity along z axis)
   !< q_pic_gpu(b,p,7) : c particle charge
   !< q_pic_gpu(b,p,8) : m particle mass
   real(R8P), pointer :: q_mhd_gpu(:,:,:,:,:) !< MHD cell centered variables.
   !< PIC center variables definition:
   !< q_mhd_gpu(b,i,j,k,v) where "b" is the block index, "i,j,k" are the cell index and
   !< "v" the variable index [blocks_number,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,nv_mhd]:
   !< q_pic_gpu(b,i,j,k,1) : Ex component (Electric field along x axis)
   !< q_pic_gpu(b,i,j,k,2) : Ey component (Electric field along y axis)
   !< q_pic_gpu(b,i,j,k,3) : Ez component (Electric field along z axis)
   !< q_pic_gpu(b,i,j,k,4) : Bx component (Magnetic field along x axis)
   !< q_pic_gpu(b,i,j,k,5) : By component (Magnetic field along y axis)
   !< q_pic_gpu(b,i,j,k,6) : Bz component (Magnetic field along z axis)
   !< q_pic_gpu(b,i,j,k,7) : Ix component (Current field along x axis)
   !< q_pic_gpu(b,i,j,k,8) : Iy component (Current field along y axis)
   !< q_pic_gpu(b,i,j,k,9) : Iz component (Current field along z axis)

   ! real(R8P), pointer :: dq_gpu(:,:,:,:,:)    !< Eikonal right hand side.
   ! real(R8P), pointer :: flx_gpu(:,:,:,:,:)   !< Fluxes along x.
   ! real(R8P), pointer :: fly_gpu(:,:,:,:,:)   !< Fluxes along y.
   ! real(R8P), pointer :: flz_gpu(:,:,:,:,:)   !< Fluxes along z.
   contains
      ! auxiliary methods
      procedure, pass(self) :: allocate_gpu !< Allocate GPU data.
      procedure, pass(self) :: copy_cpu_gpu !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu !< Copy data from GPU to CPU.
      procedure, pass(self) :: initialize   !< Initialize the equation.
endtype prism_fnl_object

contains
   ! auxiliary methods
   subroutine allocate_gpu(self)
   !< Allocate GPU data.
   class(prism_fnl_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: ierr !< Error status.

   call self%mpih%print_message('prism_fnl_object%allocate_gpu start')
   associate(nv_pic=>self%nv_pic, &
             nv_mhd=>self%nv_mhd, &
             nb=>self%nb,         &
             ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk)
   call dev_alloc(fptr_dev=self%q_pic_gpu, &
                  ubounds=[nb,np,nv_pic], lbounds=[1,1,1], init_value=0._R8P, ierr=ierr)
   call dev_alloc(fptr_dev=self%q_mhd_gpu, &
                  ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv_mhd], lbounds=[1,1-ngc,1-ngc,1-ngc,1], init_value=0._R8P, ierr=ierr)
   endassociate
   call self%mpih%print_message('prism_fnl_object%allocate_gpu finish')
   endsubroutine allocate_gpu

   subroutine copy_cpu_gpu(self)
   !< Copy data from CPU to GPU.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   call self%field_gpu%copy_transpose_cpu_gpu(nv=self%nv_mhd, q_cpu=self%field%q_mhd, q_gpu=self%q_mhd_gpu)
   !call self%field_gpu%copy_transpose_cpu_gpu(nv=self%nv_pic, q_cpu=self%field%q_pic, q_gpu=self%q_mhd_pic)
   call self%field_gpu%copy_cpu_gpu(verbose=.false.)
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self)
   !< Copy data from GPU to CPU.
   class(prism_fnl_object), intent(inout) :: self !< The equation.

   call self%field_gpu%copy_transpose_gpu_cpu(nv=self%nv_mhd, q_gpu=self%q_mhd_gpu, q_cpu=self%field%q_mhd)
   endsubroutine copy_gpu_cpu

   subroutine initialize(self, filename)
   !< Initialize the equation.
   class(prism_fnl_object), intent(inout) :: self     !< The equation.
   character(*),            intent(in)    :: filename !< Input file name.

   call self%mpih_gpu%initialize(do_mpi_init=.true., do_device_init=.true., verbose=.true.)
   call self%mpih_gpu%print_message('prism_fnl_object%initialize start')
   call self%initialize_common(filename=filename, memory_avail=real(self%mpih_gpu%dev_memory_avail,R8P), verbose=.true.)
   ! call self%field_gpu%initialize(field=self%adam%field, nv_aux=self%nv_aux, verbose=.false.)
   call self%allocate_gpu(q_gpu=self%field_gpu%q_gpu)
   call self%mpih%print_message(self%mpih_gpu%description()//new_line('a')//'prism_fnl_object%initialize finish')
   endsubroutine initialize
endmodule adam_prism_fnl_object

