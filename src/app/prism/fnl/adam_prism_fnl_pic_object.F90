!< ADAM, PRISM Particle-in-Cell definition, FNL backend.

#include "fundal.H"

module adam_prism_fnl_pic_object
!< ADAM, PRISM Particle-in-Cell definition, FNL backend.

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
public :: prism_fnl_pic_object

integer(I4P), parameter :: PIC_VARIABLES_NUMBER   = 8_I4P
integer(I4P), parameter :: PIC_FIELDS_NUMBER      = 6_I4P
integer(I4P), parameter :: PIC_NEIGHBOURS_NUMBER  = 4_I4P

type :: prism_fnl_pic_object
   real(R8P),    pointer :: q_pic_gpu(:,:)          => null() !< PIC variables on device [particle, variable].
   real(R8P),    pointer :: pic_fields_gpu(:,:)     => null() !< Fields at particle locations on device [particle, field].
   integer(I4P), pointer :: neighbour_list_gpu(:,:) => null() !< Particle-cell map on device [particle, entry].
   integer(I4P)          :: db2_q_pic(2,2)          = 0_I4P   !< Device bounds for q_pic_gpu.
   integer(I4P)          :: hb2_q_pic(2,2)          = 0_I4P   !< Host bounds for q_pic.
   integer(I4P)          :: db2_pic_fields(2,2)     = 0_I4P   !< Device bounds for pic_fields_gpu.
   integer(I4P)          :: hb2_pic_fields(2,2)     = 0_I4P   !< Host bounds for pic_fields.
   integer(I4P)          :: db2_neighbour_list(2,2) = 0_I4P   !< Device bounds for neighbour_list_gpu.
   integer(I4P)          :: hb2_neighbour_list(2,2) = 0_I4P   !< Host bounds for neighbour_list.
   real(R8P), allocatable :: buf_q_pic_R8P(:,:)           !< Host buffer with device layout for q_pic copies.
   real(R8P), allocatable :: buf_pic_fields_R8P(:,:)      !< Host buffer with device layout for pic_fields copies.
   integer(I4P), allocatable :: buf_neighbour_list_I4P(:,:) !< Host buffer with device layout for neighbour_list copies.
   real(R8P)             :: sigma        = 0.0_R8P           !< Standard deviation for Gaussian weighting.
   real(R8P)             :: cutoff_sigma = 0.0_R8P           !< Gaussian cutoff radius in sigma units.
   integer(I4P)          :: particle_number = 0_I4P          !< Total number of particles.
   integer(I4P)          :: n_ions         = 0_I4P           !< Total ions number.
   integer(I4P)          :: n_electrons    = 0_I4P           !< Total electrons number.
   integer(I4P)          :: n_neutrals     = 0_I4P           !< Total neutrals number.
   procedure(particle_cartesian_grid_index_interface_dev), pass(self), pointer :: particle_cartesian_grid_index_dev => null()
   procedure(particle_weighting_interface_dev),            pass(self), pointer :: particle_weighting_dev            => null()
   procedure(current_weighting_interface_dev),             pass(self), pointer :: current_weighting_dev             => null()
   procedure(field_weighting_interface_dev),               pass(self), pointer :: field_weighting_dev               => null()
contains
   procedure, pass(self) :: copy_cpu_gpu
   procedure, pass(self) :: copy_gpu_cpu
   procedure, pass(self) :: copy_q_pic_gpu_cpu
   procedure, pass(self) :: destroy
   procedure, pass(self) :: initialize
   procedure, pass(self) :: particle_cartesian_grid_index_dev_impl
   procedure, pass(self) :: NGP_charge_weighting_dev
   procedure, pass(self) :: CIC_charge_weighting_dev
   procedure, pass(self) :: TSC_charge_weighting_dev
   procedure, pass(self) :: cubic_charge_weighting_dev
   procedure, pass(self) :: quartic_charge_weighting_dev
   procedure, pass(self) :: quintic_charge_weighting_dev
   procedure, pass(self) :: Gaussian_charge_weighting_dev
   procedure, pass(self) :: NGP_current_weighting_dev
   procedure, pass(self) :: CIC_current_weighting_dev
   procedure, pass(self) :: TSC_current_weighting_dev
   procedure, pass(self) :: cubic_current_weighting_dev
   procedure, pass(self) :: quartic_current_weighting_dev
   procedure, pass(self) :: quintic_current_weighting_dev
   procedure, pass(self) :: Gaussian_current_weighting_dev
   procedure, pass(self) :: zeroD_field_weighting_dev
   procedure, pass(self) :: oneD_field_weighting_dev
   procedure, pass(self) :: twoD_field_weighting_dev
   procedure, pass(self) :: threeD_field_weighting_dev
   procedure, pass(self) :: fourD_field_weighting_dev
   procedure, pass(self) :: fiveD_field_weighting_dev
   procedure, pass(self) :: Gaussian_field_weighting_dev
endtype prism_fnl_pic_object

interface
   subroutine particle_cartesian_grid_index_interface_dev(self, field_fnl, field, grid, q_pic_gpu)
   import :: prism_fnl_pic_object, field_fnl_object, field_object, grid_object, R8P
   class(prism_fnl_pic_object), intent(inout) :: self             !< PIC object on device.
   type(field_fnl_object),      intent(in)    :: field_fnl        !< Device field helper.
   type(field_object),          intent(in)    :: field            !< Host field helper.
   type(grid_object),           intent(in)    :: grid             !< Host grid helper.
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:) !< PIC variables on device.
   endsubroutine particle_cartesian_grid_index_interface_dev

   subroutine particle_weighting_interface_dev(self, field_fnl, field, grid, q_gpu, q_pic_gpu, nv)
   import :: prism_fnl_pic_object, field_fnl_object, field_object, grid_object, I4P, R8P
   class(prism_fnl_pic_object), intent(inout) :: self                     !< PIC object on device.
   type(field_fnl_object),      intent(in)    :: field_fnl                !< Device field helper.
   type(field_object),          intent(in)    :: field                    !< Host field helper.
   type(grid_object),           intent(in)    :: grid                     !< Host grid helper.
   real(R8P),                   intent(inout) :: q_gpu(1:,               &
                                                      1-grid%ngc:,       &
                                                      1-grid%ngc:,       &
                                                      1-grid%ngc:,1:)    !< Field variables on device.
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)         !< PIC variables on device.
   integer(I4P),                intent(in)    :: nv                       !< Number of variables.
   endsubroutine particle_weighting_interface_dev

   subroutine current_weighting_interface_dev(self, field_fnl, field, grid, q_gpu, q_pic_gpu, nv)
   import :: prism_fnl_pic_object, field_fnl_object, field_object, grid_object, I4P, R8P
   class(prism_fnl_pic_object), intent(inout) :: self                     !< PIC object on device.
   type(field_fnl_object),      intent(in)    :: field_fnl                !< Device field helper.
   type(field_object),          intent(in)    :: field                    !< Host field helper.
   type(grid_object),           intent(in)    :: grid                     !< Host grid helper.
   real(R8P),                   intent(inout) :: q_gpu(1:,               &
                                                      1-grid%ngc:,       &
                                                      1-grid%ngc:,       &
                                                      1-grid%ngc:,1:)    !< Field variables on device.
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)         !< PIC variables on device.
   integer(I4P),                intent(in)    :: nv                       !< Number of variables.
   endsubroutine current_weighting_interface_dev

   subroutine field_weighting_interface_dev(self, field_fnl, field, grid, pic_fields_gpu, q_gpu, q_pic_gpu, nv)
   import :: prism_fnl_pic_object, field_fnl_object, field_object, grid_object, I4P, R8P
   class(prism_fnl_pic_object), intent(inout) :: self                     !< PIC object on device.
   type(field_fnl_object),      intent(in)    :: field_fnl                !< Device field helper.
   type(field_object),          intent(in)    :: field                    !< Host field helper.
   type(grid_object),           intent(in)    :: grid                     !< Host grid helper.
   real(R8P),                   intent(inout) :: pic_fields_gpu(1:,1:)    !< Fields at particle locations on device.
   real(R8P),                   intent(in)    :: q_gpu(1:,               &
                                                      1-grid%ngc:,       &
                                                      1-grid%ngc:,       &
                                                      1-grid%ngc:,1:)    !< Field variables on device.
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)         !< PIC variables on device.
   integer(I4P),                intent(in)    :: nv                       !< Number of variables.
   endsubroutine field_weighting_interface_dev
endinterface

contains
   subroutine destroy(self)
   !< Free device and host staging data owned by this object.
   class(prism_fnl_pic_object), intent(inout) :: self !< PIC object on device.

   if (associated(self%q_pic_gpu)) then
      call dev_free(self%q_pic_gpu, mydev)
      nullify(self%q_pic_gpu)
   endif
   if (associated(self%pic_fields_gpu)) then
      call dev_free(self%pic_fields_gpu, mydev)
      nullify(self%pic_fields_gpu)
   endif
   if (associated(self%neighbour_list_gpu)) then
      call dev_free(self%neighbour_list_gpu, mydev)
      nullify(self%neighbour_list_gpu)
   endif
   if (allocated(self%buf_q_pic_R8P)) deallocate(self%buf_q_pic_R8P)
   if (allocated(self%buf_pic_fields_R8P)) deallocate(self%buf_pic_fields_R8P)
   if (allocated(self%buf_neighbour_list_I4P)) deallocate(self%buf_neighbour_list_I4P)
   self%db2_q_pic = 0_I4P
   self%hb2_q_pic = 0_I4P
   self%db2_pic_fields = 0_I4P
   self%hb2_pic_fields = 0_I4P
   self%db2_neighbour_list = 0_I4P
   self%hb2_neighbour_list = 0_I4P
   self%particle_number = 0_I4P
   self%n_ions = 0_I4P
   self%n_electrons = 0_I4P
   self%n_neutrals = 0_I4P
   nullify(self%particle_cartesian_grid_index_dev)
   nullify(self%particle_weighting_dev)
   nullify(self%current_weighting_dev)
   nullify(self%field_weighting_dev)
   endsubroutine destroy

   subroutine copy_cpu_gpu(self, pic, q_pic, pic_fields, verbose)
   !< Copy PIC data from CPU to GPU.
   class(prism_fnl_pic_object), intent(inout)         :: self             !< PIC object on device.
   class(prism_pic_object),     intent(in)            :: pic              !< PIC object on host.
   real(R8P),                   intent(in), target    :: q_pic(1:,1:)     !< PIC variables on host.
   real(R8P),                   intent(in), target    :: pic_fields(1:,1:) !< Fields at particle locations on host.
   logical,                     intent(in), optional  :: verbose          !< Trigger verbose output.
   logical                                            :: verbose_         !< Trigger verbose output, local var.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (self%particle_number == 0) return
   if (verbose_) call mpih_fnl%print_message('prism_fnl_pic_object%copy_cpu_gpu start')
   call dev_memcpy_to_device(bb=self%db2_q_pic,              tb=self%hb2_q_pic,              &
                             dst=self%q_pic_gpu,             src=q_pic,                       &
                             buf=self%buf_q_pic_R8P)
   call dev_memcpy_to_device(bb=self%db2_pic_fields,         tb=self%hb2_pic_fields,         &
                             dst=self%pic_fields_gpu,        src=pic_fields,                  &
                             buf=self%buf_pic_fields_R8P)
   call dev_memcpy_to_device(bb=self%db2_neighbour_list,     tb=self%hb2_neighbour_list,     &
                             dst=self%neighbour_list_gpu,    src=pic%neighbour_list,          &
                             buf=self%buf_neighbour_list_I4P)
   if (verbose_) call mpih_fnl%print_message('prism_fnl_pic_object%copy_cpu_gpu finish')
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self, pic, q_pic, pic_fields, verbose)
   !< Copy PIC data from GPU to CPU.
   class(prism_fnl_pic_object), intent(inout)          :: self              !< PIC object on device.
   class(prism_pic_object),     intent(inout)          :: pic               !< PIC object on host.
   real(R8P),                   intent(inout), target  :: q_pic(1:,1:)      !< PIC variables on host.
   real(R8P),                   intent(inout), target  :: pic_fields(1:,1:) !< Fields at particle locations on host.
   logical,                     intent(in), optional   :: verbose           !< Trigger verbose output.
   logical                                             :: verbose_          !< Trigger verbose output, local var.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (self%particle_number == 0) return
   if (verbose_) call mpih_fnl%print_message('prism_fnl_pic_object%copy_gpu_cpu start')
   call dev_memcpy_from_device(bb=self%db2_q_pic,           tb=self%hb2_q_pic,           &
                               dst=q_pic,                   src=self%q_pic_gpu,           &
                               buf=self%buf_q_pic_R8P)
   call dev_memcpy_from_device(bb=self%db2_pic_fields,      tb=self%hb2_pic_fields,      &
                               dst=pic_fields,              src=self%pic_fields_gpu,      &
                               buf=self%buf_pic_fields_R8P)
   call dev_memcpy_from_device(bb=self%db2_neighbour_list,  tb=self%hb2_neighbour_list,  &
                               dst=pic%neighbour_list,      src=self%neighbour_list_gpu,  &
                               buf=self%buf_neighbour_list_I4P)
   if (verbose_) call mpih_fnl%print_message('prism_fnl_pic_object%copy_gpu_cpu finish')
   endsubroutine copy_gpu_cpu

   subroutine copy_q_pic_gpu_cpu(self, q_pic, verbose)
   !< Copy only PIC particle state from GPU to CPU.
   class(prism_fnl_pic_object), intent(inout)        :: self         !< PIC object on device.
   real(R8P),                   intent(inout), target :: q_pic(1:,1:) !< PIC variables on host.
   logical,                     intent(in), optional :: verbose      !< Trigger verbose output.
   logical                                           :: verbose_     !< Trigger verbose output, local var.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (self%particle_number == 0) return
   if (verbose_) call mpih_fnl%print_message('prism_fnl_pic_object%copy_q_pic_gpu_cpu start')
   call dev_memcpy_from_device(bb=self%db2_q_pic, tb=self%hb2_q_pic, dst=q_pic, src=self%q_pic_gpu, buf=self%buf_q_pic_R8P)
   if (verbose_) call mpih_fnl%print_message('prism_fnl_pic_object%copy_q_pic_gpu_cpu finish')
   endsubroutine copy_q_pic_gpu_cpu

   subroutine initialize(self, pic, q_pic, pic_fields)
   !< Initialize PIC device mirrors and dispatch hooks.
   class(prism_fnl_pic_object), intent(inout) :: self              !< PIC object on device.
   type(prism_pic_object),      intent(in)    :: pic               !< PIC object on host.
   real(R8P),                   intent(in), target :: q_pic(1:,1:)      !< PIC variables on host.
   real(R8P),                   intent(in), target :: pic_fields(1:,1:) !< Fields at particle locations on host.
   integer(I4P)                                :: ierr              !< Error status.

   call mpih_fnl%print_message('prism_fnl_pic_object%initialize start')
   call self%destroy()

   self%particle_number = pic%particle_number
   self%n_ions          = pic%n_ions
   self%n_electrons     = pic%n_electrons
   self%n_neutrals      = pic%n_neutrals
   self%sigma           = pic%sigma
   self%cutoff_sigma    = pic%cutoff_sigma

   self%particle_cartesian_grid_index_dev => particle_cartesian_grid_index_dev_impl

   select case(trim(pic%particle_weighting_model))
   case('NGP')
      self%particle_weighting_dev => NGP_charge_weighting_dev
   case('CIC')
      self%particle_weighting_dev => CIC_charge_weighting_dev
   case('TSC')
      self%particle_weighting_dev => TSC_charge_weighting_dev
   case('cubic')
      self%particle_weighting_dev => cubic_charge_weighting_dev
   case('quartic')
      self%particle_weighting_dev => quartic_charge_weighting_dev
   case('quintic')
      self%particle_weighting_dev => quintic_charge_weighting_dev
   case('Gaussian')
      self%particle_weighting_dev => Gaussian_charge_weighting_dev
   case default
      call mpih_fnl%error_stop(msg=': invalid particle weighting model in prism_fnl_pic_object%initialize')
   endselect

   select case(trim(pic%current_weighting_model))
   case('NGP')
      self%current_weighting_dev => NGP_current_weighting_dev
   case('CIC')
      self%current_weighting_dev => CIC_current_weighting_dev
   case('TSC')
      self%current_weighting_dev => TSC_current_weighting_dev
   case('cubic')
      self%current_weighting_dev => cubic_current_weighting_dev
   case('quartic')
      self%current_weighting_dev => quartic_current_weighting_dev
   case('quintic')
      self%current_weighting_dev => quintic_current_weighting_dev
   case('Gaussian')
      self%current_weighting_dev => Gaussian_current_weighting_dev
   case default
      call mpih_fnl%error_stop(msg=': invalid current weighting model in prism_fnl_pic_object%initialize')
   endselect

   select case(trim(pic%field_weighting_model))
   case('0D')
      self%field_weighting_dev => zeroD_field_weighting_dev
   case('1D')
      self%field_weighting_dev => oneD_field_weighting_dev
   case('2D')
      self%field_weighting_dev => twoD_field_weighting_dev
   case('3D')
      self%field_weighting_dev => threeD_field_weighting_dev
   case('4D')
      self%field_weighting_dev => fourD_field_weighting_dev
   case('5D')
      self%field_weighting_dev => fiveD_field_weighting_dev
   case('Gaussian')
      self%field_weighting_dev => Gaussian_field_weighting_dev
   case default
      call mpih_fnl%error_stop(msg=': invalid field weighting model in prism_fnl_pic_object%initialize')
   endselect

   if (self%particle_number == 0) then
      call mpih_fnl%print_message('prism_fnl_pic_object%initialize finish (no particles: nothing on device)')
      return
   endif

   call dev_alloc(fptr_dev=self%q_pic_gpu,          ubounds=[self%particle_number, PIC_VARIABLES_NUMBER ], &
                  lbounds=[1,1], ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate q_pic_gpu in prism_fnl_pic_object%initialize')
   call dev_alloc(fptr_dev=self%pic_fields_gpu,     ubounds=[self%particle_number, PIC_FIELDS_NUMBER    ], &
                  lbounds=[1,1], ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate pic_fields_gpu in prism_fnl_pic_object%initialize')
   call dev_alloc(fptr_dev=self%neighbour_list_gpu, ubounds=[self%particle_number, PIC_NEIGHBOURS_NUMBER], &
                  lbounds=[1,1], ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate neighbour_list_gpu in prism_fnl_pic_object%initialize')
   allocate(self%buf_q_pic_R8P(1:self%particle_number,1:PIC_VARIABLES_NUMBER), stat=ierr)
   if (ierr /= 0) call mpih_fnl%error_stop(msg=': failed host allocation of buf_q_pic_R8P in prism_fnl_pic_object%initialize')
   allocate(self%buf_pic_fields_R8P(1:self%particle_number,1:PIC_FIELDS_NUMBER), stat=ierr)
   if (ierr /= 0) call mpih_fnl%error_stop(msg=': failed host allocation of buf_pic_fields_R8P in prism_fnl_pic_object%initialize')
   allocate(self%buf_neighbour_list_I4P(1:self%particle_number,1:PIC_NEIGHBOURS_NUMBER), stat=ierr)
   if (ierr /= 0) call mpih_fnl%error_stop(msg=': failed host allocation of buf_neighbour_list_I4P in prism_fnl_pic_object%initialize')
   self%db2_q_pic(1,:)          = [1_I4P, 1_I4P]
   self%db2_q_pic(2,:)          = [self%particle_number, PIC_VARIABLES_NUMBER]
   self%hb2_q_pic(1,:)          = [1_I4P, 1_I4P]
   self%hb2_q_pic(2,:)          = [PIC_VARIABLES_NUMBER, self%particle_number]
   self%db2_pic_fields(1,:)     = [1_I4P, 1_I4P]
   self%db2_pic_fields(2,:)     = [self%particle_number, PIC_FIELDS_NUMBER]
   self%hb2_pic_fields(1,:)     = [1_I4P, 1_I4P]
   self%hb2_pic_fields(2,:)     = [PIC_FIELDS_NUMBER, self%particle_number]
   self%db2_neighbour_list(1,:) = [1_I4P, 1_I4P]
   self%db2_neighbour_list(2,:) = [self%particle_number, PIC_NEIGHBOURS_NUMBER]
   self%hb2_neighbour_list(1,:) = [1_I4P, 1_I4P]
   self%hb2_neighbour_list(2,:) = [PIC_NEIGHBOURS_NUMBER, self%particle_number]
   call self%copy_cpu_gpu(pic=pic, q_pic=q_pic, pic_fields=pic_fields)
   call mpih_fnl%print_message('prism_fnl_pic_object%initialize finish')
   endsubroutine initialize

   subroutine particle_cartesian_grid_index_dev_impl(self, field_fnl, field, grid, q_pic_gpu)
   !< Compute the grid index corresponding to a particle position. Good for cartesian grids only.
   class(prism_fnl_pic_object), intent(inout) :: self             !< PIC object on device.
   type(field_fnl_object),      intent(in)    :: field_fnl        !< Device field helper.
   type(field_object),          intent(in)    :: field            !< Host field helper.
   type(grid_object),           intent(in)    :: grid             !< Host grid helper.
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:) !< PIC variables on device.
   integer(I4P)                               :: n, b             !< Counters.
   integer(I4P)                               :: ni, nj, nk       !< Grid sizes.
   integer(I4P)                               :: blocks_number     !< Number of blocks.
   integer(I4P)                               :: i_p, j_p, k_p     !< Particle cell indices.
   integer(I4P)                               :: block_p           !< Particle block index.
   real(R8P)                                  :: x_p, y_p, z_p     !< Particle coordinates.
   real(R8P)                                  :: dx, dy, dz        !< Local cell spacings.
   real(R8P)                                  :: emin_x, emin_y, emin_z !< Block minimum coordinates.
   real(R8P), pointer                        :: x_cell_gpu(:,:)   !< Cells x coordinates on GPU.
   real(R8P), pointer                        :: y_cell_gpu(:,:)   !< Cells y coordinates on GPU.
   real(R8P), pointer                        :: z_cell_gpu(:,:)   !< Cells z coordinates on GPU.
   real(R8P), pointer                        :: dxyz_gpu(:,:)     !< Cells deltas on GPU.
   integer(I4P), pointer                     :: neighbour_list_gpu(:,:) !< Neighbour list on GPU.

   if (self%particle_number == 0) return

   ni = grid%ni
   nj = grid%nj
   nk = grid%nk
   blocks_number = field_fnl%blocks_number
   x_cell_gpu => field_fnl%x_cell_gpu
   y_cell_gpu => field_fnl%y_cell_gpu
   z_cell_gpu => field_fnl%z_cell_gpu
   dxyz_gpu   => field_fnl%dxyz_gpu
   neighbour_list_gpu => self%neighbour_list_gpu

   !$acc parallel loop independent DEVICEVAR(q_pic_gpu, x_cell_gpu, y_cell_gpu, z_cell_gpu, dxyz_gpu, neighbour_list_gpu)&
   !$acc& private(x_p, y_p, z_p, dx, dy, dz, emin_x, emin_y, emin_z, i_p, j_p, k_p, block_p)
   !$omp OMPLOOP DEVICEPTR(q_pic_gpu, x_cell_gpu, y_cell_gpu, z_cell_gpu, dxyz_gpu, neighbour_list_gpu) &
   !$omp& private(x_p, y_p, z_p, dx, dy, dz, emin_x, emin_y, emin_z, i_p, j_p, k_p, block_p)
   do n = 1, self%particle_number
      x_p = q_pic_gpu(n,1)
      y_p = q_pic_gpu(n,2)
      z_p = q_pic_gpu(n,3)

      block_p = 0_I4P
      i_p = 0_I4P ; j_p = 0_I4P ; k_p = 0_I4P

      !$acc loop seq
      do b = 1, blocks_number
         dx = dxyz_gpu(b,1)
         dy = dxyz_gpu(b,2)
         dz = dxyz_gpu(b,3)
         emin_x = x_cell_gpu(b,1) - 0.5_R8P * dx
         emin_y = y_cell_gpu(b,1) - 0.5_R8P * dy
         emin_z = z_cell_gpu(b,1) - 0.5_R8P * dz
         i_p = ceiling((x_p - emin_x) / dx, kind=I4P)
         j_p = ceiling((y_p - emin_y) / dy, kind=I4P)
         k_p = ceiling((z_p - emin_z) / dz, kind=I4P)
         if (i_p >= 1_I4P .and. i_p <= ni .and. &
             j_p >= 1_I4P .and. j_p <= nj .and. &
             k_p >= 1_I4P .and. k_p <= nk) then
            block_p = b
            exit
         endif
      enddo
      if (block_p == 0_I4P) then
         i_p = 0_I4P
         j_p = 0_I4P
         k_p = 0_I4P
      endif

      neighbour_list_gpu(n,1) = block_p
      neighbour_list_gpu(n,2) = i_p
      neighbour_list_gpu(n,3) = j_p
      neighbour_list_gpu(n,4) = k_p
   enddo
   endsubroutine particle_cartesian_grid_index_dev_impl

   subroutine NGP_charge_weighting_dev(self, field_fnl, field, grid, q_gpu, q_pic_gpu, nv)
   !< Nearest Grid Point weighting of particle charge density to the grid.
   class(prism_fnl_pic_object), intent(inout) :: self                     !< PIC object on device.
   type(field_fnl_object),      intent(in)    :: field_fnl                !< Device field helper.
   type(field_object),          intent(in)    :: field                    !< Host field helper.
   type(grid_object),           intent(in)    :: grid                     !< Host grid helper.
   real(R8P),                   intent(inout) :: q_gpu(1:,               &
                                                      1-grid%ngc:,       &
                                                      1-grid%ngc:,       &
                                                      1-grid%ngc:,1:)    !< Field variables on device.
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)         !< PIC variables on device.
   integer(I4P),                intent(in)    :: nv                       !< Charge variable index.

   call deposit_bspline_charge_dev(self=self, field_fnl=field_fnl, grid=grid, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, nv=nv, order=0_I4P)
   endsubroutine NGP_charge_weighting_dev

   subroutine CIC_charge_weighting_dev(self, field_fnl, field, grid, q_gpu, q_pic_gpu, nv)
   !< Cloud-in-Cell weighting of particle charge density to the grid.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call deposit_bspline_charge_dev(self=self, field_fnl=field_fnl, grid=grid, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, nv=nv, order=1_I4P)
   endsubroutine CIC_charge_weighting_dev

   subroutine TSC_charge_weighting_dev(self, field_fnl, field, grid, q_gpu, q_pic_gpu, nv)
   !< Triangular-Shaped-Cloud weighting of particle charge density to the grid.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call deposit_bspline_charge_dev(self=self, field_fnl=field_fnl, grid=grid, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, nv=nv, order=2_I4P)
   endsubroutine TSC_charge_weighting_dev

   subroutine cubic_charge_weighting_dev(self, field_fnl, field, grid, q_gpu, q_pic_gpu, nv)
   !< Cubic B-spline weighting of particle charge density to the grid.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call deposit_bspline_charge_dev(self=self, field_fnl=field_fnl, grid=grid, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, nv=nv, order=3_I4P)
   endsubroutine cubic_charge_weighting_dev

   subroutine quartic_charge_weighting_dev(self, field_fnl, field, grid, q_gpu, q_pic_gpu, nv)
   !< Quartic B-spline weighting of particle charge density to the grid.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call deposit_bspline_charge_dev(self=self, field_fnl=field_fnl, grid=grid, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, nv=nv, order=4_I4P)
   endsubroutine quartic_charge_weighting_dev

   subroutine quintic_charge_weighting_dev(self, field_fnl, field, grid, q_gpu, q_pic_gpu, nv)
   !< Quintic B-spline weighting of particle charge density to the grid.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call deposit_bspline_charge_dev(self=self, field_fnl=field_fnl, grid=grid, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, nv=nv, order=5_I4P)
   endsubroutine quintic_charge_weighting_dev

   subroutine Gaussian_charge_weighting_dev(self, field_fnl, field, grid, q_gpu, q_pic_gpu, nv)
   !< Gaussian weighting of particle charge density to the grid.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call deposit_gaussian_charge_dev(self=self, field_fnl=field_fnl, grid=grid, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, nv=nv)
   endsubroutine Gaussian_charge_weighting_dev

   subroutine NGP_current_weighting_dev(self, field_fnl, field, grid, q_gpu, q_pic_gpu, nv)
   !< Nearest Grid Point weighting of particle current density to the grid.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call deposit_bspline_current_dev(self=self, field_fnl=field_fnl, grid=grid, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, nv=nv, order=0_I4P)
   endsubroutine NGP_current_weighting_dev

   subroutine CIC_current_weighting_dev(self, field_fnl, field, grid, q_gpu, q_pic_gpu, nv)
   !< Cloud-in-Cell weighting of particle current density to the grid.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call deposit_bspline_current_dev(self=self, field_fnl=field_fnl, grid=grid, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, nv=nv, order=1_I4P)
   endsubroutine CIC_current_weighting_dev

   subroutine TSC_current_weighting_dev(self, field_fnl, field, grid, q_gpu, q_pic_gpu, nv)
   !< Triangular-Shaped-Cloud weighting of particle current density to the grid.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call deposit_bspline_current_dev(self=self, field_fnl=field_fnl, grid=grid, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, nv=nv, order=2_I4P)
   endsubroutine TSC_current_weighting_dev

   subroutine cubic_current_weighting_dev(self, field_fnl, field, grid, q_gpu, q_pic_gpu, nv)
   !< Cubic B-spline weighting of particle current density to the grid.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call deposit_bspline_current_dev(self=self, field_fnl=field_fnl, grid=grid, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, nv=nv, order=3_I4P)
   endsubroutine cubic_current_weighting_dev

   subroutine quartic_current_weighting_dev(self, field_fnl, field, grid, q_gpu, q_pic_gpu, nv)
   !< Quartic B-spline weighting of particle current density to the grid.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call deposit_bspline_current_dev(self=self, field_fnl=field_fnl, grid=grid, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, nv=nv, order=4_I4P)
   endsubroutine quartic_current_weighting_dev

   subroutine quintic_current_weighting_dev(self, field_fnl, field, grid, q_gpu, q_pic_gpu, nv)
   !< Quintic B-spline weighting of particle current density to the grid.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call deposit_bspline_current_dev(self=self, field_fnl=field_fnl, grid=grid, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, nv=nv, order=5_I4P)
   endsubroutine quintic_current_weighting_dev

   subroutine Gaussian_current_weighting_dev(self, field_fnl, field, grid, q_gpu, q_pic_gpu, nv)
   !< Gaussian weighting of particle current density to the grid.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call deposit_gaussian_current_dev(self=self, field_fnl=field_fnl, grid=grid, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, nv=nv)
   endsubroutine Gaussian_current_weighting_dev

   subroutine zeroD_field_weighting_dev(self, field_fnl, field, grid, pic_fields_gpu, q_gpu, q_pic_gpu, nv)
   !< Zeroth-order spatial interpolation of cell-centered fields to particle locations.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: pic_fields_gpu(1:,1:)
   real(R8P),                   intent(in)    :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call gather_bspline_fields_dev(self=self, field_fnl=field_fnl, grid=grid, pic_fields_gpu=pic_fields_gpu, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, order=0_I4P)
   endsubroutine zeroD_field_weighting_dev

   subroutine oneD_field_weighting_dev(self, field_fnl, field, grid, pic_fields_gpu, q_gpu, q_pic_gpu, nv)
   !< First-order spatial interpolation of cell-centered fields to particle locations.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: pic_fields_gpu(1:,1:)
   real(R8P),                   intent(in)    :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call gather_bspline_fields_dev(self=self, field_fnl=field_fnl, grid=grid, pic_fields_gpu=pic_fields_gpu, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, order=1_I4P)
   endsubroutine oneD_field_weighting_dev

   subroutine twoD_field_weighting_dev(self, field_fnl, field, grid, pic_fields_gpu, q_gpu, q_pic_gpu, nv)
   !< Second-order spatial interpolation of cell-centered fields to particle locations.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: pic_fields_gpu(1:,1:)
   real(R8P),                   intent(in)    :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call gather_bspline_fields_dev(self=self, field_fnl=field_fnl, grid=grid, pic_fields_gpu=pic_fields_gpu, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, order=2_I4P)
   endsubroutine twoD_field_weighting_dev

   subroutine threeD_field_weighting_dev(self, field_fnl, field, grid, pic_fields_gpu, q_gpu, q_pic_gpu, nv)
   !< Third-order spatial interpolation of cell-centered fields to particle locations.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: pic_fields_gpu(1:,1:)
   real(R8P),                   intent(in)    :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call gather_bspline_fields_dev(self=self, field_fnl=field_fnl, grid=grid, pic_fields_gpu=pic_fields_gpu, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, order=3_I4P)
   endsubroutine threeD_field_weighting_dev

   subroutine fourD_field_weighting_dev(self, field_fnl, field, grid, pic_fields_gpu, q_gpu, q_pic_gpu, nv)
   !< Fourth-order spatial interpolation of cell-centered fields to particle locations.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: pic_fields_gpu(1:,1:)
   real(R8P),                   intent(in)    :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call gather_bspline_fields_dev(self=self, field_fnl=field_fnl, grid=grid, pic_fields_gpu=pic_fields_gpu, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, order=4_I4P)
   endsubroutine fourD_field_weighting_dev

   subroutine fiveD_field_weighting_dev(self, field_fnl, field, grid, pic_fields_gpu, q_gpu, q_pic_gpu, nv)
   !< Fifth-order spatial interpolation of cell-centered fields to particle locations.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: pic_fields_gpu(1:,1:)
   real(R8P),                   intent(in)    :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call gather_bspline_fields_dev(self=self, field_fnl=field_fnl, grid=grid, pic_fields_gpu=pic_fields_gpu, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu, order=5_I4P)
   endsubroutine fiveD_field_weighting_dev

   subroutine Gaussian_field_weighting_dev(self, field_fnl, field, grid, pic_fields_gpu, q_gpu, q_pic_gpu, nv)
   !< Gaussian interpolation of cell-centered fields to particle locations.
   class(prism_fnl_pic_object), intent(inout) :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(field_object),          intent(in)    :: field
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: pic_fields_gpu(1:,1:)
   real(R8P),                   intent(in)    :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv

   call gather_gaussian_fields_dev(self=self, field_fnl=field_fnl, grid=grid, pic_fields_gpu=pic_fields_gpu, q_gpu=q_gpu, q_pic_gpu=q_pic_gpu)
   endsubroutine Gaussian_field_weighting_dev

   subroutine deposit_bspline_charge_dev(self, field_fnl, grid, q_gpu, q_pic_gpu, nv, order)
   !< Deposit particle charge density with centered cardinal B-splines.
   class(prism_fnl_pic_object), intent(in)    :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv
   integer(I4P),                intent(in)    :: order
   integer(I4P)                               :: b, i, j, k, n
   integer(I4P)                               :: ngc, ni, nj, nk, blocks_number
   integer(I4P)                               :: block_p, i_p, j_p, k_p
   integer(I4P)                               :: i_min, i_max, j_min, j_max, k_min, k_max
   real(R8P)                                  :: dx, dy, dz, charge_density
   real(R8P)                                  :: wx, wy, wz, weight
   real(R8P),    pointer                      :: x_cell_gpu(:,:), y_cell_gpu(:,:), z_cell_gpu(:,:), dxyz_gpu(:,:)
   integer(I4P), pointer                      :: neighbour_list_gpu(:,:)

   if (self%particle_number == 0) return

   ngc = grid%ngc ; ni = grid%ni ; nj = grid%nj ; nk = grid%nk
   blocks_number = field_fnl%blocks_number
   x_cell_gpu => field_fnl%x_cell_gpu
   y_cell_gpu => field_fnl%y_cell_gpu
   z_cell_gpu => field_fnl%z_cell_gpu
   dxyz_gpu   => field_fnl%dxyz_gpu
   neighbour_list_gpu => self%neighbour_list_gpu

   !$acc parallel loop collapse(4) independent DEVICEVAR(q_gpu)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu)
   do k = 1-ngc, nk+ngc
      do j = 1-ngc, nj+ngc
         do i = 1-ngc, ni+ngc
            do b = 1, blocks_number
               q_gpu(b,i,j,k,nv) = 0.0_R8P
            enddo
         enddo
      enddo
   enddo

   !$acc parallel loop independent DEVICEVAR(q_gpu, q_pic_gpu, x_cell_gpu, y_cell_gpu, z_cell_gpu, dxyz_gpu, neighbour_list_gpu)&
   !$acc& private(block_p, i_p, j_p, k_p, i_min, i_max, j_min, j_max, k_min, k_max, dx, dy, dz, charge_density, wx, wy, wz, weight)
   !$omp OMPLOOP DEVICEPTR(q_gpu, q_pic_gpu, x_cell_gpu, y_cell_gpu, z_cell_gpu, dxyz_gpu, neighbour_list_gpu) &
   !$omp& private(block_p, i_p, j_p, k_p, i_min, i_max, j_min, j_max, k_min, k_max, dx, dy, dz, charge_density, wx, wy, wz, weight)
   do n = 1, self%particle_number
      block_p = neighbour_list_gpu(n,1)
      if (block_p <= 0_I4P) cycle
      i_p = neighbour_list_gpu(n,2)
      j_p = neighbour_list_gpu(n,3)
      k_p = neighbour_list_gpu(n,4)

      dx = dxyz_gpu(block_p,1)
      dy = dxyz_gpu(block_p,2)
      dz = dxyz_gpu(block_p,3)
      charge_density = q_pic_gpu(n,7) / (dx * dy * dz)

      call set_bspline_stencil_dev(order=order, x_p=q_pic_gpu(n,1), x_c=x_cell_gpu(block_p,i_p), i_p=i_p, i_min=i_min, i_max=i_max)
      call set_bspline_stencil_dev(order=order, x_p=q_pic_gpu(n,2), x_c=y_cell_gpu(block_p,j_p), i_p=j_p, i_min=j_min, i_max=j_max)
      call set_bspline_stencil_dev(order=order, x_p=q_pic_gpu(n,3), x_c=z_cell_gpu(block_p,k_p), i_p=k_p, i_min=k_min, i_max=k_max)

      i_min = max(i_min, 1-ngc) ; i_max = min(i_max, ni+ngc)
      j_min = max(j_min, 1-ngc) ; j_max = min(j_max, nj+ngc)
      k_min = max(k_min, 1-ngc) ; k_max = min(k_max, nk+ngc)

      !$acc loop seq
      do k = k_min, k_max
         wz = bspline_weight_dev(order=order, r=(q_pic_gpu(n,3) - z_cell_gpu(block_p,k)) / dz)
         !$acc loop seq
         do j = j_min, j_max
            wy = bspline_weight_dev(order=order, r=(q_pic_gpu(n,2) - y_cell_gpu(block_p,j)) / dy)
            !$acc loop seq
            do i = i_min, i_max
               wx = bspline_weight_dev(order=order, r=(q_pic_gpu(n,1) - x_cell_gpu(block_p,i)) / dx)
               weight = wx * wy * wz
               !$acc atomic update
               !$omp atomic update
               q_gpu(block_p,i,j,k,nv) = q_gpu(block_p,i,j,k,nv) + charge_density * weight
            enddo
         enddo
      enddo
   enddo
   endsubroutine deposit_bspline_charge_dev

   subroutine deposit_gaussian_charge_dev(self, field_fnl, grid, q_gpu, q_pic_gpu, nv)
   !< Deposit particle charge density with Gaussian support.
   class(prism_fnl_pic_object), intent(in)    :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv
   integer(I4P)                               :: b, i, j, k, n
   integer(I4P)                               :: ngc, ni, nj, nk, blocks_number
   integer(I4P)                               :: block_p, i_p, j_p, k_p
   integer(I4P)                               :: i_min, i_max, j_min, j_max, k_min, k_max
   integer(I4P)                               :: ni_sigma, nj_sigma, nk_sigma
   real(R8P)                                  :: dx, dy, dz, sigma_x, sigma_y, sigma_z
   real(R8P)                                  :: inverse_cell_volume, charge_prefactor
   real(R8P)                                  :: rx, ry, rz, wx, wy, wz, weight, weight_sum
   real(R8P),    pointer                      :: x_cell_gpu(:,:), y_cell_gpu(:,:), z_cell_gpu(:,:), dxyz_gpu(:,:)
   integer(I4P), pointer                      :: neighbour_list_gpu(:,:)

   if (self%particle_number == 0) return

   ngc = grid%ngc ; ni = grid%ni ; nj = grid%nj ; nk = grid%nk
   blocks_number = field_fnl%blocks_number
   x_cell_gpu => field_fnl%x_cell_gpu
   y_cell_gpu => field_fnl%y_cell_gpu
   z_cell_gpu => field_fnl%z_cell_gpu
   dxyz_gpu   => field_fnl%dxyz_gpu
   neighbour_list_gpu => self%neighbour_list_gpu

   !$acc parallel loop collapse(4) independent DEVICEVAR(q_gpu)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu)
   do k = 1-ngc, nk+ngc
      do j = 1-ngc, nj+ngc
         do i = 1-ngc, ni+ngc
            do b = 1, blocks_number
               q_gpu(b,i,j,k,nv) = 0.0_R8P
            enddo
         enddo
      enddo
   enddo

   !$acc parallel loop independent DEVICEVAR(q_gpu, q_pic_gpu, x_cell_gpu, y_cell_gpu, z_cell_gpu, dxyz_gpu, neighbour_list_gpu)&
   !$acc& private(block_p, i_p, j_p, k_p, i_min, i_max, j_min, j_max, k_min, k_max, ni_sigma, nj_sigma, nk_sigma, &
   !$acc&         dx, dy, dz, sigma_x, sigma_y, sigma_z, inverse_cell_volume, charge_prefactor, rx, ry, rz, wx, wy, wz, weight, weight_sum)
   !$omp OMPLOOP DEVICEPTR(q_gpu, q_pic_gpu, x_cell_gpu, y_cell_gpu, z_cell_gpu, dxyz_gpu, neighbour_list_gpu) &
   !$omp& private(block_p, i_p, j_p, k_p, i_min, i_max, j_min, j_max, k_min, k_max, ni_sigma, nj_sigma, nk_sigma, &
   !$omp&         dx, dy, dz, sigma_x, sigma_y, sigma_z, inverse_cell_volume, charge_prefactor, rx, ry, rz, wx, wy, wz, weight, weight_sum)
   do n = 1, self%particle_number
      block_p = neighbour_list_gpu(n,1)
      if (block_p <= 0_I4P) cycle
      i_p = neighbour_list_gpu(n,2)
      j_p = neighbour_list_gpu(n,3)
      k_p = neighbour_list_gpu(n,4)

      dx = dxyz_gpu(block_p,1)
      dy = dxyz_gpu(block_p,2)
      dz = dxyz_gpu(block_p,3)
      sigma_x = self%sigma
      sigma_y = self%sigma
      sigma_z = self%sigma
      inverse_cell_volume = 1.0_R8P / (dx * dy * dz)
      charge_prefactor = q_pic_gpu(n,7) * inverse_cell_volume

      ni_sigma = ceiling(self%cutoff_sigma * sigma_x / dx, kind=I4P)
      nj_sigma = ceiling(self%cutoff_sigma * sigma_y / dy, kind=I4P)
      nk_sigma = ceiling(self%cutoff_sigma * sigma_z / dz, kind=I4P)

      i_min = max(i_p - ni_sigma, 1-ngc) ; i_max = min(i_p + ni_sigma, ni+ngc)
      j_min = max(j_p - nj_sigma, 1-ngc) ; j_max = min(j_p + nj_sigma, nj+ngc)
      k_min = max(k_p - nk_sigma, 1-ngc) ; k_max = min(k_p + nk_sigma, nk+ngc)

      weight_sum = 0.0_R8P
      !$acc loop seq
      do k = k_min, k_max
         rz = (q_pic_gpu(n,3) - z_cell_gpu(block_p,k)) / sigma_z
         wz = exp(-0.5_R8P * rz * rz)
         !$acc loop seq
         do j = j_min, j_max
            ry = (q_pic_gpu(n,2) - y_cell_gpu(block_p,j)) / sigma_y
            wy = exp(-0.5_R8P * ry * ry)
            !$acc loop seq
            do i = i_min, i_max
               rx = (q_pic_gpu(n,1) - x_cell_gpu(block_p,i)) / sigma_x
               wx = exp(-0.5_R8P * rx * rx)
               weight_sum = weight_sum + wx * wy * wz
            enddo
         enddo
      enddo

      if (weight_sum > tiny(1.0_R8P)) then
         !$acc loop seq
         do k = k_min, k_max
            rz = (q_pic_gpu(n,3) - z_cell_gpu(block_p,k)) / sigma_z
            wz = exp(-0.5_R8P * rz * rz)
            !$acc loop seq
            do j = j_min, j_max
               ry = (q_pic_gpu(n,2) - y_cell_gpu(block_p,j)) / sigma_y
               wy = exp(-0.5_R8P * ry * ry)
               !$acc loop seq
               do i = i_min, i_max
                  rx = (q_pic_gpu(n,1) - x_cell_gpu(block_p,i)) / sigma_x
                  wx = exp(-0.5_R8P * rx * rx)
                  weight = wx * wy * wz / weight_sum
                  !$acc atomic update
                  !$omp atomic update
                  q_gpu(block_p,i,j,k,nv) = q_gpu(block_p,i,j,k,nv) + charge_prefactor * weight
               enddo
            enddo
         enddo
      endif
   enddo
   endsubroutine deposit_gaussian_charge_dev

   subroutine deposit_bspline_current_dev(self, field_fnl, grid, q_gpu, q_pic_gpu, nv, order)
   !< Deposit particle current density with centered cardinal B-splines.
   class(prism_fnl_pic_object), intent(in)    :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv
   integer(I4P),                intent(in)    :: order
   integer(I4P)                               :: b, i, j, k, n
   integer(I4P)                               :: ngc, ni, nj, nk, blocks_number
   integer(I4P)                               :: block_p, i_p, j_p, k_p
   integer(I4P)                               :: i_min, i_max, j_min, j_max, k_min, k_max
   real(R8P)                                  :: dx, dy, dz, prefactor
   real(R8P)                                  :: wx, wy, wz, weight
   real(R8P)                                  :: jx, jy, jz
   real(R8P),    pointer                      :: x_cell_gpu(:,:), y_cell_gpu(:,:), z_cell_gpu(:,:), dxyz_gpu(:,:)
   integer(I4P), pointer                      :: neighbour_list_gpu(:,:)

   if (self%particle_number == 0) return

   ngc = grid%ngc ; ni = grid%ni ; nj = grid%nj ; nk = grid%nk
   blocks_number = field_fnl%blocks_number
   x_cell_gpu => field_fnl%x_cell_gpu
   y_cell_gpu => field_fnl%y_cell_gpu
   z_cell_gpu => field_fnl%z_cell_gpu
   dxyz_gpu   => field_fnl%dxyz_gpu
   neighbour_list_gpu => self%neighbour_list_gpu

   !$acc parallel loop collapse(4) independent DEVICEVAR(q_gpu)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu)
   do k = 1-ngc, nk+ngc
      do j = 1-ngc, nj+ngc
         do i = 1-ngc, ni+ngc
            do b = 1, blocks_number
               q_gpu(b,i,j,k,nv-3) = 0.0_R8P
               q_gpu(b,i,j,k,nv-2) = 0.0_R8P
               q_gpu(b,i,j,k,nv-1) = 0.0_R8P
            enddo
         enddo
      enddo
   enddo

   !$acc parallel loop independent DEVICEVAR(q_gpu, q_pic_gpu, x_cell_gpu, y_cell_gpu, z_cell_gpu, dxyz_gpu, neighbour_list_gpu)&
   !$acc& private(block_p, i_p, j_p, k_p, i_min, i_max, j_min, j_max, k_min, k_max, dx, dy, dz, prefactor, wx, wy, wz, weight, jx, jy, jz)
   !$omp OMPLOOP DEVICEPTR(q_gpu, q_pic_gpu, x_cell_gpu, y_cell_gpu, z_cell_gpu, dxyz_gpu, neighbour_list_gpu) &
   !$omp& private(block_p, i_p, j_p, k_p, i_min, i_max, j_min, j_max, k_min, k_max, dx, dy, dz, prefactor, wx, wy, wz, weight, jx, jy, jz)
   do n = 1, self%particle_number
      block_p = neighbour_list_gpu(n,1)
      if (block_p <= 0_I4P) cycle
      i_p = neighbour_list_gpu(n,2)
      j_p = neighbour_list_gpu(n,3)
      k_p = neighbour_list_gpu(n,4)

      dx = dxyz_gpu(block_p,1)
      dy = dxyz_gpu(block_p,2)
      dz = dxyz_gpu(block_p,3)
      prefactor = q_pic_gpu(n,7) / (dx * dy * dz)
      jx = prefactor * q_pic_gpu(n,4)
      jy = prefactor * q_pic_gpu(n,5)
      jz = prefactor * q_pic_gpu(n,6)

      call set_bspline_stencil_dev(order=order, x_p=q_pic_gpu(n,1), x_c=x_cell_gpu(block_p,i_p), i_p=i_p, i_min=i_min, i_max=i_max)
      call set_bspline_stencil_dev(order=order, x_p=q_pic_gpu(n,2), x_c=y_cell_gpu(block_p,j_p), i_p=j_p, i_min=j_min, i_max=j_max)
      call set_bspline_stencil_dev(order=order, x_p=q_pic_gpu(n,3), x_c=z_cell_gpu(block_p,k_p), i_p=k_p, i_min=k_min, i_max=k_max)

      i_min = max(i_min, 1-ngc) ; i_max = min(i_max, ni+ngc)
      j_min = max(j_min, 1-ngc) ; j_max = min(j_max, nj+ngc)
      k_min = max(k_min, 1-ngc) ; k_max = min(k_max, nk+ngc)

      !$acc loop seq
      do k = k_min, k_max
         wz = bspline_weight_dev(order=order, r=(q_pic_gpu(n,3) - z_cell_gpu(block_p,k)) / dz)
         !$acc loop seq
         do j = j_min, j_max
            wy = bspline_weight_dev(order=order, r=(q_pic_gpu(n,2) - y_cell_gpu(block_p,j)) / dy)
            !$acc loop seq
            do i = i_min, i_max
               wx = bspline_weight_dev(order=order, r=(q_pic_gpu(n,1) - x_cell_gpu(block_p,i)) / dx)
               weight = wx * wy * wz
               !$acc atomic update
               !$omp atomic update
               q_gpu(block_p,i,j,k,nv-3) = q_gpu(block_p,i,j,k,nv-3) + jx * weight
               !$acc atomic update
               !$omp atomic update
               q_gpu(block_p,i,j,k,nv-2) = q_gpu(block_p,i,j,k,nv-2) + jy * weight
               !$acc atomic update
               !$omp atomic update
               q_gpu(block_p,i,j,k,nv-1) = q_gpu(block_p,i,j,k,nv-1) + jz * weight
            enddo
         enddo
      enddo
   enddo
   endsubroutine deposit_bspline_current_dev

   subroutine deposit_gaussian_current_dev(self, field_fnl, grid, q_gpu, q_pic_gpu, nv)
   !< Deposit particle current density with Gaussian support.
   class(prism_fnl_pic_object), intent(in)    :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: nv
   integer(I4P)                               :: b, i, j, k, n
   integer(I4P)                               :: ngc, ni, nj, nk, blocks_number
   integer(I4P)                               :: block_p, i_p, j_p, k_p
   integer(I4P)                               :: i_min, i_max, j_min, j_max, k_min, k_max
   integer(I4P)                               :: ni_sigma, nj_sigma, nk_sigma
   real(R8P)                                  :: dx, dy, dz, sigma_x, sigma_y, sigma_z, inverse_cell_volume
   real(R8P)                                  :: rx, ry, rz, wx, wy, wz, weight, weight_sum
   real(R8P)                                  :: jx, jy, jz
   real(R8P),    pointer                      :: x_cell_gpu(:,:), y_cell_gpu(:,:), z_cell_gpu(:,:), dxyz_gpu(:,:)
   integer(I4P), pointer                      :: neighbour_list_gpu(:,:)

   if (self%particle_number == 0) return

   ngc = grid%ngc ; ni = grid%ni ; nj = grid%nj ; nk = grid%nk
   blocks_number = field_fnl%blocks_number
   x_cell_gpu => field_fnl%x_cell_gpu
   y_cell_gpu => field_fnl%y_cell_gpu
   z_cell_gpu => field_fnl%z_cell_gpu
   dxyz_gpu   => field_fnl%dxyz_gpu
   neighbour_list_gpu => self%neighbour_list_gpu

   !$acc parallel loop collapse(4) independent DEVICEVAR(q_gpu)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu)
   do k = 1-ngc, nk+ngc
      do j = 1-ngc, nj+ngc
         do i = 1-ngc, ni+ngc
            do b = 1, blocks_number
               q_gpu(b,i,j,k,nv-3) = 0.0_R8P
               q_gpu(b,i,j,k,nv-2) = 0.0_R8P
               q_gpu(b,i,j,k,nv-1) = 0.0_R8P
            enddo
         enddo
      enddo
   enddo

   !$acc parallel loop independent DEVICEVAR(q_gpu, q_pic_gpu, x_cell_gpu, y_cell_gpu, z_cell_gpu, dxyz_gpu, neighbour_list_gpu)&
   !$acc& private(block_p, i_p, j_p, k_p, i_min, i_max, j_min, j_max, k_min, k_max, ni_sigma, nj_sigma, nk_sigma, &
   !$acc&         dx, dy, dz, sigma_x, sigma_y, sigma_z, inverse_cell_volume, rx, ry, rz, wx, wy, wz, weight, weight_sum, jx, jy, jz)
   !$omp OMPLOOP DEVICEPTR(q_gpu, q_pic_gpu, x_cell_gpu, y_cell_gpu, z_cell_gpu, dxyz_gpu, neighbour_list_gpu) &
   !$omp& private(block_p, i_p, j_p, k_p, i_min, i_max, j_min, j_max, k_min, k_max, ni_sigma, nj_sigma, nk_sigma, &
   !$omp&         dx, dy, dz, sigma_x, sigma_y, sigma_z, inverse_cell_volume, rx, ry, rz, wx, wy, wz, weight, weight_sum, jx, jy, jz)
   do n = 1, self%particle_number
      block_p = neighbour_list_gpu(n,1)
      if (block_p <= 0_I4P) cycle
      i_p = neighbour_list_gpu(n,2)
      j_p = neighbour_list_gpu(n,3)
      k_p = neighbour_list_gpu(n,4)

      dx = dxyz_gpu(block_p,1)
      dy = dxyz_gpu(block_p,2)
      dz = dxyz_gpu(block_p,3)
      sigma_x = self%sigma
      sigma_y = self%sigma
      sigma_z = self%sigma
      inverse_cell_volume = 1.0_R8P / (dx * dy * dz)
      jx = q_pic_gpu(n,7) * q_pic_gpu(n,4) * inverse_cell_volume
      jy = q_pic_gpu(n,7) * q_pic_gpu(n,5) * inverse_cell_volume
      jz = q_pic_gpu(n,7) * q_pic_gpu(n,6) * inverse_cell_volume

      ni_sigma = ceiling(self%cutoff_sigma * sigma_x / dx, kind=I4P)
      nj_sigma = ceiling(self%cutoff_sigma * sigma_y / dy, kind=I4P)
      nk_sigma = ceiling(self%cutoff_sigma * sigma_z / dz, kind=I4P)

      i_min = max(i_p - ni_sigma, 1-ngc) ; i_max = min(i_p + ni_sigma, ni+ngc)
      j_min = max(j_p - nj_sigma, 1-ngc) ; j_max = min(j_p + nj_sigma, nj+ngc)
      k_min = max(k_p - nk_sigma, 1-ngc) ; k_max = min(k_p + nk_sigma, nk+ngc)

      weight_sum = 0.0_R8P
      !$acc loop seq
      do k = k_min, k_max
         rz = (q_pic_gpu(n,3) - z_cell_gpu(block_p,k)) / sigma_z
         wz = exp(-0.5_R8P * rz * rz)
         !$acc loop seq
         do j = j_min, j_max
            ry = (q_pic_gpu(n,2) - y_cell_gpu(block_p,j)) / sigma_y
            wy = exp(-0.5_R8P * ry * ry)
            !$acc loop seq
            do i = i_min, i_max
               rx = (q_pic_gpu(n,1) - x_cell_gpu(block_p,i)) / sigma_x
               wx = exp(-0.5_R8P * rx * rx)
               weight_sum = weight_sum + wx * wy * wz
            enddo
         enddo
      enddo

      if (weight_sum > tiny(1.0_R8P)) then
         !$acc loop seq
         do k = k_min, k_max
            rz = (q_pic_gpu(n,3) - z_cell_gpu(block_p,k)) / sigma_z
            wz = exp(-0.5_R8P * rz * rz)
            !$acc loop seq
            do j = j_min, j_max
               ry = (q_pic_gpu(n,2) - y_cell_gpu(block_p,j)) / sigma_y
               wy = exp(-0.5_R8P * ry * ry)
               !$acc loop seq
               do i = i_min, i_max
                  rx = (q_pic_gpu(n,1) - x_cell_gpu(block_p,i)) / sigma_x
                  wx = exp(-0.5_R8P * rx * rx)
                  weight = wx * wy * wz / weight_sum
                  !$acc atomic update
                  !$omp atomic update
                  q_gpu(block_p,i,j,k,nv-3) = q_gpu(block_p,i,j,k,nv-3) + jx * weight
                  !$acc atomic update
                  !$omp atomic update
                  q_gpu(block_p,i,j,k,nv-2) = q_gpu(block_p,i,j,k,nv-2) + jy * weight
                  !$acc atomic update
                  !$omp atomic update
                  q_gpu(block_p,i,j,k,nv-1) = q_gpu(block_p,i,j,k,nv-1) + jz * weight
               enddo
            enddo
         enddo
      endif
   enddo
   endsubroutine deposit_gaussian_current_dev

   subroutine gather_bspline_fields_dev(self, field_fnl, grid, pic_fields_gpu, q_gpu, q_pic_gpu, order)
   !< Gather cell-centered fields at particle locations with centered cardinal B-splines.
   class(prism_fnl_pic_object), intent(in)    :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: pic_fields_gpu(1:,1:)
   real(R8P),                   intent(in)    :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P),                intent(in)    :: order
   integer(I4P)                               :: i, j, k, n
   integer(I4P)                               :: ngc, ni, nj, nk
   integer(I4P)                               :: block_p, i_p, j_p, k_p
   integer(I4P)                               :: i_min, i_max, j_min, j_max, k_min, k_max
   real(R8P)                                  :: dx, dy, dz, wx, wy, wz, weight
   real(R8P)                                  :: f1, f2, f3, f4, f5, f6
   real(R8P),    pointer                      :: x_cell_gpu(:,:), y_cell_gpu(:,:), z_cell_gpu(:,:), dxyz_gpu(:,:)
   integer(I4P), pointer                      :: neighbour_list_gpu(:,:)

   if (self%particle_number == 0) return

   ngc = grid%ngc ; ni = grid%ni ; nj = grid%nj ; nk = grid%nk
   x_cell_gpu => field_fnl%x_cell_gpu
   y_cell_gpu => field_fnl%y_cell_gpu
   z_cell_gpu => field_fnl%z_cell_gpu
   dxyz_gpu   => field_fnl%dxyz_gpu
   neighbour_list_gpu => self%neighbour_list_gpu

   !$acc parallel loop independent DEVICEVAR(pic_fields_gpu, q_gpu, q_pic_gpu, x_cell_gpu, y_cell_gpu, z_cell_gpu, dxyz_gpu, neighbour_list_gpu)&
   !$acc& private(block_p, i_p, j_p, k_p, i_min, i_max, j_min, j_max, k_min, k_max, dx, dy, dz, wx, wy, wz, weight, f1, f2, f3, f4, f5, f6)
   !$omp OMPLOOP DEVICEPTR(pic_fields_gpu, q_gpu, q_pic_gpu, x_cell_gpu, y_cell_gpu, z_cell_gpu, dxyz_gpu, neighbour_list_gpu) &
   !$omp& private(block_p, i_p, j_p, k_p, i_min, i_max, j_min, j_max, k_min, k_max, dx, dy, dz, wx, wy, wz, weight, f1, f2, f3, f4, f5, f6)
   do n = 1, self%particle_number
      block_p = neighbour_list_gpu(n,1)
      if (block_p <= 0_I4P) then
         pic_fields_gpu(n,1) = 0.0_R8P
         pic_fields_gpu(n,2) = 0.0_R8P
         pic_fields_gpu(n,3) = 0.0_R8P
         pic_fields_gpu(n,4) = 0.0_R8P
         pic_fields_gpu(n,5) = 0.0_R8P
         pic_fields_gpu(n,6) = 0.0_R8P
         cycle
      endif

      i_p = neighbour_list_gpu(n,2)
      j_p = neighbour_list_gpu(n,3)
      k_p = neighbour_list_gpu(n,4)

      dx = dxyz_gpu(block_p,1)
      dy = dxyz_gpu(block_p,2)
      dz = dxyz_gpu(block_p,3)

      call set_bspline_stencil_dev(order=order, x_p=q_pic_gpu(n,1), x_c=x_cell_gpu(block_p,i_p), i_p=i_p, i_min=i_min, i_max=i_max)
      call set_bspline_stencil_dev(order=order, x_p=q_pic_gpu(n,2), x_c=y_cell_gpu(block_p,j_p), i_p=j_p, i_min=j_min, i_max=j_max)
      call set_bspline_stencil_dev(order=order, x_p=q_pic_gpu(n,3), x_c=z_cell_gpu(block_p,k_p), i_p=k_p, i_min=k_min, i_max=k_max)

      i_min = max(i_min, 1-ngc) ; i_max = min(i_max, ni+ngc)
      j_min = max(j_min, 1-ngc) ; j_max = min(j_max, nj+ngc)
      k_min = max(k_min, 1-ngc) ; k_max = min(k_max, nk+ngc)

      f1 = 0.0_R8P ; f2 = 0.0_R8P ; f3 = 0.0_R8P
      f4 = 0.0_R8P ; f5 = 0.0_R8P ; f6 = 0.0_R8P

      !$acc loop seq
      do k = k_min, k_max
         wz = bspline_weight_dev(order=order, r=(q_pic_gpu(n,3) - z_cell_gpu(block_p,k)) / dz)
         !$acc loop seq
         do j = j_min, j_max
            wy = bspline_weight_dev(order=order, r=(q_pic_gpu(n,2) - y_cell_gpu(block_p,j)) / dy)
            !$acc loop seq
            do i = i_min, i_max
               wx = bspline_weight_dev(order=order, r=(q_pic_gpu(n,1) - x_cell_gpu(block_p,i)) / dx)
               weight = wx * wy * wz
               f1 = f1 + weight * q_gpu(block_p,i,j,k,1)
               f2 = f2 + weight * q_gpu(block_p,i,j,k,2)
               f3 = f3 + weight * q_gpu(block_p,i,j,k,3)
               f4 = f4 + weight * q_gpu(block_p,i,j,k,4)
               f5 = f5 + weight * q_gpu(block_p,i,j,k,5)
               f6 = f6 + weight * q_gpu(block_p,i,j,k,6)
            enddo
         enddo
      enddo

      pic_fields_gpu(n,1) = f1
      pic_fields_gpu(n,2) = f2
      pic_fields_gpu(n,3) = f3
      pic_fields_gpu(n,4) = f4
      pic_fields_gpu(n,5) = f5
      pic_fields_gpu(n,6) = f6
   enddo
   endsubroutine gather_bspline_fields_dev

   subroutine gather_gaussian_fields_dev(self, field_fnl, grid, pic_fields_gpu, q_gpu, q_pic_gpu)
   !< Gather cell-centered fields at particle locations with Gaussian support.
   class(prism_fnl_pic_object), intent(in)    :: self
   type(field_fnl_object),      intent(in)    :: field_fnl
   type(grid_object),           intent(in)    :: grid
   real(R8P),                   intent(inout) :: pic_fields_gpu(1:,1:)
   real(R8P),                   intent(in)    :: q_gpu(1:,1-grid%ngc:,1-grid%ngc:,1-grid%ngc:,1:)
   real(R8P),                   intent(in)    :: q_pic_gpu(1:,1:)
   integer(I4P)                               :: i, j, k, n
   integer(I4P)                               :: ngc, ni, nj, nk
   integer(I4P)                               :: block_p, i_p, j_p, k_p
   integer(I4P)                               :: i_min, i_max, j_min, j_max, k_min, k_max
   integer(I4P)                               :: ni_sigma, nj_sigma, nk_sigma
   real(R8P)                                  :: dx, dy, dz, sigma_x, sigma_y, sigma_z
   real(R8P)                                  :: rx, ry, rz, wx, wy, wz, weight, weight_sum
   real(R8P)                                  :: f1, f2, f3, f4, f5, f6
   real(R8P),    pointer                      :: x_cell_gpu(:,:), y_cell_gpu(:,:), z_cell_gpu(:,:), dxyz_gpu(:,:)
   integer(I4P), pointer                      :: neighbour_list_gpu(:,:)

   if (self%particle_number == 0) return

   ngc = grid%ngc ; ni = grid%ni ; nj = grid%nj ; nk = grid%nk
   x_cell_gpu => field_fnl%x_cell_gpu
   y_cell_gpu => field_fnl%y_cell_gpu
   z_cell_gpu => field_fnl%z_cell_gpu
   dxyz_gpu   => field_fnl%dxyz_gpu
   neighbour_list_gpu => self%neighbour_list_gpu

   !$acc parallel loop independent DEVICEVAR(pic_fields_gpu, q_gpu, q_pic_gpu, x_cell_gpu, y_cell_gpu, z_cell_gpu, dxyz_gpu, neighbour_list_gpu)&
   !$acc& private(block_p, i_p, j_p, k_p, i_min, i_max, j_min, j_max, k_min, k_max, ni_sigma, nj_sigma, nk_sigma, &
   !$acc&         dx, dy, dz, sigma_x, sigma_y, sigma_z, rx, ry, rz, wx, wy, wz, weight, weight_sum, f1, f2, f3, f4, f5, f6)
   !$omp OMPLOOP DEVICEPTR(pic_fields_gpu, q_gpu, q_pic_gpu, x_cell_gpu, y_cell_gpu, z_cell_gpu, dxyz_gpu, neighbour_list_gpu) &
   !$omp& private(block_p, i_p, j_p, k_p, i_min, i_max, j_min, j_max, k_min, k_max, ni_sigma, nj_sigma, nk_sigma, &
   !$omp&         dx, dy, dz, sigma_x, sigma_y, sigma_z, rx, ry, rz, wx, wy, wz, weight, weight_sum, f1, f2, f3, f4, f5, f6)
   do n = 1, self%particle_number
      block_p = neighbour_list_gpu(n,1)
      if (block_p <= 0_I4P) then
         pic_fields_gpu(n,1) = 0.0_R8P
         pic_fields_gpu(n,2) = 0.0_R8P
         pic_fields_gpu(n,3) = 0.0_R8P
         pic_fields_gpu(n,4) = 0.0_R8P
         pic_fields_gpu(n,5) = 0.0_R8P
         pic_fields_gpu(n,6) = 0.0_R8P
         cycle
      endif

      i_p = neighbour_list_gpu(n,2)
      j_p = neighbour_list_gpu(n,3)
      k_p = neighbour_list_gpu(n,4)

      dx = dxyz_gpu(block_p,1)
      dy = dxyz_gpu(block_p,2)
      dz = dxyz_gpu(block_p,3)
      sigma_x = self%sigma
      sigma_y = self%sigma
      sigma_z = self%sigma

      ni_sigma = ceiling(self%cutoff_sigma * sigma_x / dx, kind=I4P)
      nj_sigma = ceiling(self%cutoff_sigma * sigma_y / dy, kind=I4P)
      nk_sigma = ceiling(self%cutoff_sigma * sigma_z / dz, kind=I4P)

      i_min = max(i_p - ni_sigma, 1-ngc) ; i_max = min(i_p + ni_sigma, ni+ngc)
      j_min = max(j_p - nj_sigma, 1-ngc) ; j_max = min(j_p + nj_sigma, nj+ngc)
      k_min = max(k_p - nk_sigma, 1-ngc) ; k_max = min(k_p + nk_sigma, nk+ngc)

      weight_sum = 0.0_R8P
      !$acc loop seq
      do k = k_min, k_max
         rz = (q_pic_gpu(n,3) - z_cell_gpu(block_p,k)) / sigma_z
         wz = exp(-0.5_R8P * rz * rz)
         !$acc loop seq
         do j = j_min, j_max
            ry = (q_pic_gpu(n,2) - y_cell_gpu(block_p,j)) / sigma_y
            wy = exp(-0.5_R8P * ry * ry)
            !$acc loop seq
            do i = i_min, i_max
               rx = (q_pic_gpu(n,1) - x_cell_gpu(block_p,i)) / sigma_x
               wx = exp(-0.5_R8P * rx * rx)
               weight_sum = weight_sum + wx * wy * wz
            enddo
         enddo
      enddo

      f1 = 0.0_R8P ; f2 = 0.0_R8P ; f3 = 0.0_R8P
      f4 = 0.0_R8P ; f5 = 0.0_R8P ; f6 = 0.0_R8P

      if (weight_sum > tiny(1.0_R8P)) then
         !$acc loop seq
         do k = k_min, k_max
            rz = (q_pic_gpu(n,3) - z_cell_gpu(block_p,k)) / sigma_z
            wz = exp(-0.5_R8P * rz * rz)
            !$acc loop seq
            do j = j_min, j_max
               ry = (q_pic_gpu(n,2) - y_cell_gpu(block_p,j)) / sigma_y
               wy = exp(-0.5_R8P * ry * ry)
               !$acc loop seq
               do i = i_min, i_max
                  rx = (q_pic_gpu(n,1) - x_cell_gpu(block_p,i)) / sigma_x
                  wx = exp(-0.5_R8P * rx * rx)
                  weight = wx * wy * wz / weight_sum
                  f1 = f1 + weight * q_gpu(block_p,i,j,k,1)
                  f2 = f2 + weight * q_gpu(block_p,i,j,k,2)
                  f3 = f3 + weight * q_gpu(block_p,i,j,k,3)
                  f4 = f4 + weight * q_gpu(block_p,i,j,k,4)
                  f5 = f5 + weight * q_gpu(block_p,i,j,k,5)
                  f6 = f6 + weight * q_gpu(block_p,i,j,k,6)
               enddo
            enddo
         enddo
      endif

      pic_fields_gpu(n,1) = f1
      pic_fields_gpu(n,2) = f2
      pic_fields_gpu(n,3) = f3
      pic_fields_gpu(n,4) = f4
      pic_fields_gpu(n,5) = f5
      pic_fields_gpu(n,6) = f6
   enddo
   endsubroutine gather_gaussian_fields_dev

   pure subroutine set_bspline_stencil_dev(order, x_p, x_c, i_p, i_min, i_max)
   !< Compute the one-dimensional B-spline stencil.
   !$acc routine seq
   !$omp declare target
   integer(I4P), intent(in)  :: order
   real(R8P),    intent(in)  :: x_p, x_c
   integer(I4P), intent(in)  :: i_p
   integer(I4P), intent(out) :: i_min, i_max

   select case(order)
   case(0_I4P)
      i_min = i_p
      i_max = i_p
   case(1_I4P)
      if (x_p <= x_c) then
         i_min = i_p - 1_I4P
         i_max = i_p
      else
         i_min = i_p
         i_max = i_p + 1_I4P
      endif
   case(2_I4P)
      i_min = i_p - 1_I4P
      i_max = i_p + 1_I4P
   case(3_I4P)
      if (x_p <= x_c) then
         i_min = i_p - 2_I4P
         i_max = i_p + 1_I4P
      else
         i_min = i_p - 1_I4P
         i_max = i_p + 2_I4P
      endif
   case(4_I4P)
      i_min = i_p - 2_I4P
      i_max = i_p + 2_I4P
   case(5_I4P)
      if (x_p <= x_c) then
         i_min = i_p - 3_I4P
         i_max = i_p + 2_I4P
      else
         i_min = i_p - 2_I4P
         i_max = i_p + 3_I4P
      endif
   case default
      i_min = i_p
      i_max = i_p
   endselect
   endsubroutine set_bspline_stencil_dev

   pure function bspline_weight_dev(order, r) result(weight)
   !< Return the centered cardinal B-spline weight.
   !$acc routine seq
   !$omp declare target
   integer(I4P), intent(in) :: order
   real(R8P),    intent(in) :: r
   real(R8P)                :: weight
   real(R8P)                :: a

   a = abs(r)

   select case(order)
   case(0_I4P)
      weight = 1.0_R8P
   case(1_I4P)
      if (a <= 1.0_R8P) then
         weight = 1.0_R8P - a
      else
         weight = 0.0_R8P
      endif
   case(2_I4P)
      if (a <= 0.5_R8P) then
         weight = 0.75_R8P - a * a
      elseif (a <= 1.5_R8P) then
         weight = 0.5_R8P * (1.5_R8P - a)**2
      else
         weight = 0.0_R8P
      endif
   case(3_I4P)
      if (a <= 1.0_R8P) then
         weight = (4.0_R8P - 6.0_R8P * a * a + 3.0_R8P * a**3) / 6.0_R8P
      elseif (a <= 2.0_R8P) then
         weight = (2.0_R8P - a)**3 / 6.0_R8P
      else
         weight = 0.0_R8P
      endif
   case(4_I4P)
      if (a <= 0.5_R8P) then
         weight = ((2.5_R8P - a)**4                                      &
                  - 5.0_R8P  * (1.5_R8P - a)**4                          &
                  + 10.0_R8P * (0.5_R8P - a)**4) / 24.0_R8P
      elseif (a <= 1.5_R8P) then
         weight = ((2.5_R8P - a)**4 - 5.0_R8P * (1.5_R8P - a)**4) / 24.0_R8P
      elseif (a <= 2.5_R8P) then
         weight = (2.5_R8P - a)**4 / 24.0_R8P
      else
         weight = 0.0_R8P
      endif
   case(5_I4P)
      if (a <= 1.0_R8P) then
         weight = ((3.0_R8P - a)**5                                      &
                  - 6.0_R8P  * (2.0_R8P - a)**5                          &
                  + 15.0_R8P * (1.0_R8P - a)**5) / 120.0_R8P
      elseif (a <= 2.0_R8P) then
         weight = ((3.0_R8P - a)**5 - 6.0_R8P * (2.0_R8P - a)**5) / 120.0_R8P
      elseif (a <= 3.0_R8P) then
         weight = (3.0_R8P - a)**5 / 120.0_R8P
      else
         weight = 0.0_R8P
      endif
   case default
      weight = 0.0_R8P
   endselect
   endfunction bspline_weight_dev
endmodule adam_prism_fnl_pic_object
