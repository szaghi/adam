!< ADAM, PRISM perfectly matched layer class definition, FNL backend.

#include "fundal.H"

module adam_prism_fnl_pml_object
!< ADAM, PRISM perfectly matched layer class definition, FNL backend.

use :: adam_prism_common_library, only : prism_pml_object
use :: adam_grid_object,         only : grid_object
use :: adam_field_object,        only : field_object
use :: adam_fnl_mpih_global, only : mpih_fnl
use :: fundal
use :: penf

implicit none
private
public :: prism_fnl_pml_object

integer(I4P), parameter :: PML_FACE_X_M = 1_I4P
integer(I4P), parameter :: PML_FACE_X_P = 2_I4P
integer(I4P), parameter :: PML_FACE_Y_M = 3_I4P
integer(I4P), parameter :: PML_FACE_Y_P = 4_I4P
integer(I4P), parameter :: PML_FACE_Z_M = 5_I4P
integer(I4P), parameter :: PML_FACE_Z_P = 6_I4P
integer(I4P), parameter :: PML_VARS_PER_FACE = 4_I4P
integer(I4P), parameter :: CFS_PROFILE_EXPONENT = 2_I4P
real(R8P),    parameter :: BERMUDEZ_EPS = 1.e-6_R8P
character(len=8), parameter :: PML_TYPE_NONE     = 'NONE'
character(len=8), parameter :: PML_TYPE_CLASSIC  = 'CLASSIC'
character(len=8), parameter :: PML_TYPE_CFS      = 'CFS'
character(len=8), parameter :: PML_TYPE_BERMUDEZ = 'BERMUDEZ'

type :: prism_fnl_pml_object
   logical                   :: enabled = .false.
   character(len=99)         :: pml_type = PML_TYPE_NONE
   logical                   :: layer(6) = .false.
   real(R8P)                 :: width = 0._R8P
   real(R8P)                 :: gamma_max = 0._R8P
   real(R8P)                 :: gamma_exponent = 0._R8P
   real(R8P)                 :: alpha_max = 0._R8P
   real(R8P)                 :: k_max = 1._R8P
   real(R8P)                 :: beta = 0._R8P
   integer(I4P)              :: active_blocks(6) = 0_I4P
   integer(I4P)              :: max_cells(6) = 0_I4P
   integer(I4P), pointer     :: blocks_x_m_gpu(:) => null()
   integer(I4P), pointer     :: blocks_x_p_gpu(:) => null()
   integer(I4P), pointer     :: blocks_y_m_gpu(:) => null()
   integer(I4P), pointer     :: blocks_y_p_gpu(:) => null()
   integer(I4P), pointer     :: blocks_z_m_gpu(:) => null()
   integer(I4P), pointer     :: blocks_z_p_gpu(:) => null()
   integer(I4P), pointer     :: start_x_m_gpu(:)  => null()
   integer(I4P), pointer     :: start_x_p_gpu(:)  => null()
   integer(I4P), pointer     :: start_y_m_gpu(:)  => null()
   integer(I4P), pointer     :: start_y_p_gpu(:)  => null()
   integer(I4P), pointer     :: start_z_m_gpu(:)  => null()
   integer(I4P), pointer     :: start_z_p_gpu(:)  => null()
   integer(I4P), pointer     :: cells_x_m_gpu(:)  => null()
   integer(I4P), pointer     :: cells_x_p_gpu(:)  => null()
   integer(I4P), pointer     :: cells_y_m_gpu(:)  => null()
   integer(I4P), pointer     :: cells_y_p_gpu(:)  => null()
   integer(I4P), pointer     :: cells_z_m_gpu(:)  => null()
   integer(I4P), pointer     :: cells_z_p_gpu(:)  => null()
   real(R8P),    pointer     :: gamma_x_m_gpu(:,:) => null()
   real(R8P),    pointer     :: gamma_x_p_gpu(:,:) => null()
   real(R8P),    pointer     :: gamma_y_m_gpu(:,:) => null()
   real(R8P),    pointer     :: gamma_y_p_gpu(:,:) => null()
   real(R8P),    pointer     :: gamma_z_m_gpu(:,:) => null()
   real(R8P),    pointer     :: gamma_z_p_gpu(:,:) => null()
   real(R8P),    pointer     :: alpha_x_m_gpu(:,:) => null()
   real(R8P),    pointer     :: alpha_x_p_gpu(:,:) => null()
   real(R8P),    pointer     :: alpha_y_m_gpu(:,:) => null()
   real(R8P),    pointer     :: alpha_y_p_gpu(:,:) => null()
   real(R8P),    pointer     :: alpha_z_m_gpu(:,:) => null()
   real(R8P),    pointer     :: alpha_z_p_gpu(:,:) => null()
   real(R8P),    pointer     :: kappa_x_m_gpu(:,:) => null()
   real(R8P),    pointer     :: kappa_x_p_gpu(:,:) => null()
   real(R8P),    pointer     :: kappa_y_m_gpu(:,:) => null()
   real(R8P),    pointer     :: kappa_y_p_gpu(:,:) => null()
   real(R8P),    pointer     :: kappa_z_m_gpu(:,:) => null()
   real(R8P),    pointer     :: kappa_z_p_gpu(:,:) => null()
   real(R8P),    pointer     :: q_pml_x_m_gpu(:,:,:,:,:) => null()
   real(R8P),    pointer     :: q_pml_x_p_gpu(:,:,:,:,:) => null()
   real(R8P),    pointer     :: q_pml_y_m_gpu(:,:,:,:,:) => null()
   real(R8P),    pointer     :: q_pml_y_p_gpu(:,:,:,:,:) => null()
   real(R8P),    pointer     :: q_pml_z_m_gpu(:,:,:,:,:) => null()
   real(R8P),    pointer     :: q_pml_z_p_gpu(:,:,:,:,:) => null()
contains
   procedure, pass(self) :: initialize
endtype prism_fnl_pml_object

contains

   subroutine initialize(self, pml, grid, field)
   !< Initialize the GPU-side PML support from the host compact PML object.
   class(prism_fnl_pml_object), intent(inout) :: self
   type(prism_pml_object),      intent(in)    :: pml
   type(grid_object),           intent(in)    :: grid
   type(field_object),          intent(in)    :: field

   self%enabled         = pml%enabled
   self%pml_type        = pml%pml_type
   self%layer           = pml%layer
   self%width           = pml%width
   self%gamma_max       = pml%gamma_max
   self%gamma_exponent  = pml%gamma_exponent
   self%alpha_max       = pml%alpha_max
   self%k_max           = pml%k_max
   self%beta            = pml%beta
   self%active_blocks   = pml%active_blocks
   self%max_cells       = pml%max_cells

   if (.not. self%enabled) return

   call initialize_x_face(pml_type=self%pml_type, width=self%width, gamma_max=self%gamma_max,                       &
                          gamma_exponent=self%gamma_exponent, alpha_max=self%alpha_max, k_max=self%k_max,          &
                          beta=self%beta, face=PML_FACE_X_M, q_face_cpu=pml%q_pml_x_m, blocks_cpu=pml%blocks_x_m,  &
                          range_cpu=pml%ni_pml, profile_span=pml%profile_span(PML_FACE_X_M), grid=grid, field=field, &
                          blocks_gpu=self%blocks_x_m_gpu,            &
                          start_gpu=self%start_x_m_gpu, cells_gpu=self%cells_x_m_gpu, gamma_gpu=self%gamma_x_m_gpu,&
                          alpha_gpu=self%alpha_x_m_gpu, kappa_gpu=self%kappa_x_m_gpu, q_face_gpu=self%q_pml_x_m_gpu)
   call initialize_x_face(pml_type=self%pml_type, width=self%width, gamma_max=self%gamma_max,                       &
                          gamma_exponent=self%gamma_exponent, alpha_max=self%alpha_max, k_max=self%k_max,          &
                          beta=self%beta, face=PML_FACE_X_P, q_face_cpu=pml%q_pml_x_p, blocks_cpu=pml%blocks_x_p,  &
                          range_cpu=pml%ni_pml, profile_span=pml%profile_span(PML_FACE_X_P), grid=grid, field=field, &
                          blocks_gpu=self%blocks_x_p_gpu,            &
                          start_gpu=self%start_x_p_gpu, cells_gpu=self%cells_x_p_gpu, gamma_gpu=self%gamma_x_p_gpu,&
                          alpha_gpu=self%alpha_x_p_gpu, kappa_gpu=self%kappa_x_p_gpu, q_face_gpu=self%q_pml_x_p_gpu)
   call initialize_y_face(pml_type=self%pml_type, width=self%width, gamma_max=self%gamma_max,                       &
                          gamma_exponent=self%gamma_exponent, alpha_max=self%alpha_max, k_max=self%k_max,          &
                          beta=self%beta, face=PML_FACE_Y_M, q_face_cpu=pml%q_pml_y_m, blocks_cpu=pml%blocks_y_m,  &
                          range_cpu=pml%nj_pml, profile_span=pml%profile_span(PML_FACE_Y_M), grid=grid, field=field, &
                          blocks_gpu=self%blocks_y_m_gpu,            &
                          start_gpu=self%start_y_m_gpu, cells_gpu=self%cells_y_m_gpu, gamma_gpu=self%gamma_y_m_gpu,&
                          alpha_gpu=self%alpha_y_m_gpu, kappa_gpu=self%kappa_y_m_gpu, q_face_gpu=self%q_pml_y_m_gpu)
   call initialize_y_face(pml_type=self%pml_type, width=self%width, gamma_max=self%gamma_max,                       &
                          gamma_exponent=self%gamma_exponent, alpha_max=self%alpha_max, k_max=self%k_max,          &
                          beta=self%beta, face=PML_FACE_Y_P, q_face_cpu=pml%q_pml_y_p, blocks_cpu=pml%blocks_y_p,  &
                          range_cpu=pml%nj_pml, profile_span=pml%profile_span(PML_FACE_Y_P), grid=grid, field=field, &
                          blocks_gpu=self%blocks_y_p_gpu,            &
                          start_gpu=self%start_y_p_gpu, cells_gpu=self%cells_y_p_gpu, gamma_gpu=self%gamma_y_p_gpu,&
                          alpha_gpu=self%alpha_y_p_gpu, kappa_gpu=self%kappa_y_p_gpu, q_face_gpu=self%q_pml_y_p_gpu)
   call initialize_z_face(pml_type=self%pml_type, width=self%width, gamma_max=self%gamma_max,                       &
                          gamma_exponent=self%gamma_exponent, alpha_max=self%alpha_max, k_max=self%k_max,          &
                          beta=self%beta, face=PML_FACE_Z_M, q_face_cpu=pml%q_pml_z_m, blocks_cpu=pml%blocks_z_m,  &
                          range_cpu=pml%nk_pml, profile_span=pml%profile_span(PML_FACE_Z_M), grid=grid, field=field, &
                          blocks_gpu=self%blocks_z_m_gpu,            &
                          start_gpu=self%start_z_m_gpu, cells_gpu=self%cells_z_m_gpu, gamma_gpu=self%gamma_z_m_gpu,&
                          alpha_gpu=self%alpha_z_m_gpu, kappa_gpu=self%kappa_z_m_gpu, q_face_gpu=self%q_pml_z_m_gpu)
   call initialize_z_face(pml_type=self%pml_type, width=self%width, gamma_max=self%gamma_max,                       &
                          gamma_exponent=self%gamma_exponent, alpha_max=self%alpha_max, k_max=self%k_max,          &
                          beta=self%beta, face=PML_FACE_Z_P, q_face_cpu=pml%q_pml_z_p, blocks_cpu=pml%blocks_z_p,  &
                          range_cpu=pml%nk_pml, profile_span=pml%profile_span(PML_FACE_Z_P), grid=grid, field=field, &
                          blocks_gpu=self%blocks_z_p_gpu,            &
                          start_gpu=self%start_z_p_gpu, cells_gpu=self%cells_z_p_gpu, gamma_gpu=self%gamma_z_p_gpu,&
                          alpha_gpu=self%alpha_z_p_gpu, kappa_gpu=self%kappa_z_p_gpu, q_face_gpu=self%q_pml_z_p_gpu)
   endsubroutine initialize

   subroutine initialize_x_face(pml_type, width, gamma_max, gamma_exponent, alpha_max, k_max, beta, face, q_face_cpu, &
                                blocks_cpu, range_cpu, profile_span, grid, field, blocks_gpu, start_gpu, cells_gpu, gamma_gpu, alpha_gpu, &
                                kappa_gpu, q_face_gpu)
   character(len=*),       intent(in)    :: pml_type
   real(R8P),              intent(in)    :: width, gamma_max, gamma_exponent, alpha_max, k_max, beta
   integer(I4P),            intent(in)    :: face
   real(R8P), allocatable,  intent(in)    :: q_face_cpu(:,:,:,:,:)
   integer(I4P), allocatable, intent(in)  :: blocks_cpu(:)
   integer(I4P),            intent(in)    :: range_cpu(:,:,:)
   real(R8P),               intent(in)    :: profile_span
   type(grid_object),       intent(in)    :: grid
   type(field_object),      intent(in)    :: field
   integer(I4P), pointer,   intent(inout) :: blocks_gpu(:)
   integer(I4P), pointer,   intent(inout) :: start_gpu(:)
   integer(I4P), pointer,   intent(inout) :: cells_gpu(:)
   real(R8P),    pointer,   intent(inout) :: gamma_gpu(:,:)
   real(R8P),    pointer,   intent(inout) :: alpha_gpu(:,:)
   real(R8P),    pointer,   intent(inout) :: kappa_gpu(:,:)
   real(R8P),    pointer,   intent(inout) :: q_face_gpu(:,:,:,:,:)
   integer(I4P)                           :: b, cmax, ierr, lid, nblocks
   integer(I4P), allocatable              :: start_host(:), cells_host(:)
   real(R8P), allocatable                 :: gamma_host(:,:), alpha_host(:,:), kappa_host(:,:)

   if (.not. allocated(q_face_cpu)) return
   nblocks = size(q_face_cpu, dim=5)
   cmax = size(q_face_cpu, dim=2)
   allocate(start_host(1:nblocks), cells_host(1:nblocks))
   allocate(gamma_host(1:nblocks,1:cmax), alpha_host(1:nblocks,1:cmax), kappa_host(1:nblocks,1:cmax))
   gamma_host = 0._R8P ; alpha_host = 0._R8P ; kappa_host = 1._R8P
   do lid = 1, nblocks
      b = blocks_cpu(lid)
      start_host(lid) = range_cpu(1,b,face)
      cells_host(lid) = range_cpu(2,b,face) - start_host(lid) + 1_I4P
      call fill_face_coefficients(pml_type=pml_type, width=width, gamma_max=gamma_max,                                    &
                                  gamma_exponent=gamma_exponent, alpha_max=alpha_max, k_max=k_max, beta=beta,            &
                                  face=face, profile_span=profile_span, grid=grid, field=field, block_id=b,             &
                                  start_idx=start_host(lid), cells=cells_host(lid), gamma=gamma_host(lid,:),             &
                                  alpha=alpha_host(lid,:), kappa=kappa_host(lid,:))
   enddo
   call dev_assign_to_device(src=blocks_cpu,  dst=blocks_gpu)
   call dev_assign_to_device(src=start_host,  dst=start_gpu)
   call dev_assign_to_device(src=cells_host,  dst=cells_gpu)
   call dev_assign_to_device(src=gamma_host,  dst=gamma_gpu)
   call dev_assign_to_device(src=alpha_host,  dst=alpha_gpu)
   call dev_assign_to_device(src=kappa_host,  dst=kappa_gpu)
   call dev_alloc(fptr_dev=q_face_gpu, ubounds=[nblocks,cmax,grid%nj,grid%nk,PML_VARS_PER_FACE], lbounds=[1,1,1,1,1], &
                  init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate x-face GPU PML state')
   call copy_face_state_cpu_gpu(q_cpu=q_face_cpu, q_gpu=q_face_gpu)
   endsubroutine initialize_x_face

   subroutine initialize_y_face(pml_type, width, gamma_max, gamma_exponent, alpha_max, k_max, beta, face, q_face_cpu, &
                                blocks_cpu, range_cpu, profile_span, grid, field, blocks_gpu, start_gpu, cells_gpu, gamma_gpu, alpha_gpu, &
                                kappa_gpu, q_face_gpu)
   character(len=*),       intent(in)    :: pml_type
   real(R8P),              intent(in)    :: width, gamma_max, gamma_exponent, alpha_max, k_max, beta
   integer(I4P),            intent(in)    :: face
   real(R8P), allocatable,  intent(in)    :: q_face_cpu(:,:,:,:,:)
   integer(I4P), allocatable, intent(in)  :: blocks_cpu(:)
   integer(I4P),            intent(in)    :: range_cpu(:,:,:)
   real(R8P),               intent(in)    :: profile_span
   type(grid_object),       intent(in)    :: grid
   type(field_object),      intent(in)    :: field
   integer(I4P), pointer,   intent(inout) :: blocks_gpu(:)
   integer(I4P), pointer,   intent(inout) :: start_gpu(:)
   integer(I4P), pointer,   intent(inout) :: cells_gpu(:)
   real(R8P),    pointer,   intent(inout) :: gamma_gpu(:,:)
   real(R8P),    pointer,   intent(inout) :: alpha_gpu(:,:)
   real(R8P),    pointer,   intent(inout) :: kappa_gpu(:,:)
   real(R8P),    pointer,   intent(inout) :: q_face_gpu(:,:,:,:,:)
   integer(I4P)                           :: b, cmax, ierr, lid, nblocks
   integer(I4P), allocatable              :: start_host(:), cells_host(:)
   real(R8P), allocatable                 :: gamma_host(:,:), alpha_host(:,:), kappa_host(:,:)

   if (.not. allocated(q_face_cpu)) return
   nblocks = size(q_face_cpu, dim=5)
   cmax = size(q_face_cpu, dim=3)
   allocate(start_host(1:nblocks), cells_host(1:nblocks))
   allocate(gamma_host(1:nblocks,1:cmax), alpha_host(1:nblocks,1:cmax), kappa_host(1:nblocks,1:cmax))
   gamma_host = 0._R8P ; alpha_host = 0._R8P ; kappa_host = 1._R8P
   do lid = 1, nblocks
      b = blocks_cpu(lid)
      start_host(lid) = range_cpu(1,b,face)
      cells_host(lid) = range_cpu(2,b,face) - start_host(lid) + 1_I4P
      call fill_face_coefficients(pml_type=pml_type, width=width, gamma_max=gamma_max,                                    &
                                  gamma_exponent=gamma_exponent, alpha_max=alpha_max, k_max=k_max, beta=beta,            &
                                  face=face, profile_span=profile_span, grid=grid, field=field, block_id=b,             &
                                  start_idx=start_host(lid), cells=cells_host(lid), gamma=gamma_host(lid,:),             &
                                  alpha=alpha_host(lid,:), kappa=kappa_host(lid,:))
   enddo
   call dev_assign_to_device(src=blocks_cpu,  dst=blocks_gpu)
   call dev_assign_to_device(src=start_host,  dst=start_gpu)
   call dev_assign_to_device(src=cells_host,  dst=cells_gpu)
   call dev_assign_to_device(src=gamma_host,  dst=gamma_gpu)
   call dev_assign_to_device(src=alpha_host,  dst=alpha_gpu)
   call dev_assign_to_device(src=kappa_host,  dst=kappa_gpu)
   call dev_alloc(fptr_dev=q_face_gpu, ubounds=[nblocks,grid%ni,cmax,grid%nk,PML_VARS_PER_FACE], lbounds=[1,1,1,1,1], &
                  init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate y-face GPU PML state')
   call copy_face_state_cpu_gpu(q_cpu=q_face_cpu, q_gpu=q_face_gpu)
   endsubroutine initialize_y_face

   subroutine initialize_z_face(pml_type, width, gamma_max, gamma_exponent, alpha_max, k_max, beta, face, q_face_cpu, &
                                blocks_cpu, range_cpu, profile_span, grid, field, blocks_gpu, start_gpu, cells_gpu, gamma_gpu, alpha_gpu, &
                                kappa_gpu, q_face_gpu)
   character(len=*),       intent(in)    :: pml_type
   real(R8P),              intent(in)    :: width, gamma_max, gamma_exponent, alpha_max, k_max, beta
   integer(I4P),            intent(in)    :: face
   real(R8P), allocatable,  intent(in)    :: q_face_cpu(:,:,:,:,:)
   integer(I4P), allocatable, intent(in)  :: blocks_cpu(:)
   integer(I4P),            intent(in)    :: range_cpu(:,:,:)
   real(R8P),               intent(in)    :: profile_span
   type(grid_object),       intent(in)    :: grid
   type(field_object),      intent(in)    :: field
   integer(I4P), pointer,   intent(inout) :: blocks_gpu(:)
   integer(I4P), pointer,   intent(inout) :: start_gpu(:)
   integer(I4P), pointer,   intent(inout) :: cells_gpu(:)
   real(R8P),    pointer,   intent(inout) :: gamma_gpu(:,:)
   real(R8P),    pointer,   intent(inout) :: alpha_gpu(:,:)
   real(R8P),    pointer,   intent(inout) :: kappa_gpu(:,:)
   real(R8P),    pointer,   intent(inout) :: q_face_gpu(:,:,:,:,:)
   integer(I4P)                           :: b, cmax, ierr, lid, nblocks
   integer(I4P), allocatable              :: start_host(:), cells_host(:)
   real(R8P), allocatable                 :: gamma_host(:,:), alpha_host(:,:), kappa_host(:,:)

   if (.not. allocated(q_face_cpu)) return
   nblocks = size(q_face_cpu, dim=5)
   cmax = size(q_face_cpu, dim=4)
   allocate(start_host(1:nblocks), cells_host(1:nblocks))
   allocate(gamma_host(1:nblocks,1:cmax), alpha_host(1:nblocks,1:cmax), kappa_host(1:nblocks,1:cmax))
   gamma_host = 0._R8P ; alpha_host = 0._R8P ; kappa_host = 1._R8P
   do lid = 1, nblocks
      b = blocks_cpu(lid)
      start_host(lid) = range_cpu(1,b,face)
      cells_host(lid) = range_cpu(2,b,face) - start_host(lid) + 1_I4P
      call fill_face_coefficients(pml_type=pml_type, width=width, gamma_max=gamma_max,                                    &
                                  gamma_exponent=gamma_exponent, alpha_max=alpha_max, k_max=k_max, beta=beta,            &
                                  face=face, profile_span=profile_span, grid=grid, field=field, block_id=b,             &
                                  start_idx=start_host(lid), cells=cells_host(lid), gamma=gamma_host(lid,:),             &
                                  alpha=alpha_host(lid,:), kappa=kappa_host(lid,:))
   enddo
   call dev_assign_to_device(src=blocks_cpu,  dst=blocks_gpu)
   call dev_assign_to_device(src=start_host,  dst=start_gpu)
   call dev_assign_to_device(src=cells_host,  dst=cells_gpu)
   call dev_assign_to_device(src=gamma_host,  dst=gamma_gpu)
   call dev_assign_to_device(src=alpha_host,  dst=alpha_gpu)
   call dev_assign_to_device(src=kappa_host,  dst=kappa_gpu)
   call dev_alloc(fptr_dev=q_face_gpu, ubounds=[nblocks,grid%ni,grid%nj,cmax,PML_VARS_PER_FACE], lbounds=[1,1,1,1,1], &
                  init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate z-face GPU PML state')
   call copy_face_state_cpu_gpu(q_cpu=q_face_cpu, q_gpu=q_face_gpu)
   endsubroutine initialize_z_face

   subroutine copy_face_state_cpu_gpu(q_cpu, q_gpu)
   real(R8P), intent(in)    :: q_cpu(:,:,:,:,:)
   real(R8P), intent(inout) :: q_gpu(:,:,:,:,:)
   integer(I4P)             :: i1, i2, i3, i4, i5
   real(R8P), allocatable   :: q_t(:,:,:,:,:)

   allocate(q_t(1:size(q_cpu,5), 1:size(q_cpu,2), 1:size(q_cpu,3), 1:size(q_cpu,4), 1:size(q_cpu,1)))
   do i5 = 1, size(q_cpu,5)
      do i4 = 1, size(q_cpu,4)
         do i3 = 1, size(q_cpu,3)
            do i2 = 1, size(q_cpu,2)
               do i1 = 1, size(q_cpu,1)
                  q_t(i5,i2,i3,i4,i1) = q_cpu(i1,i2,i3,i4,i5)
               enddo
            enddo
         enddo
      enddo
   enddo
   call dev_memcpy_to_device(dst=q_gpu, src=q_t)
   deallocate(q_t)
   endsubroutine copy_face_state_cpu_gpu

   subroutine fill_face_coefficients(pml_type, width, gamma_max, gamma_exponent, alpha_max, k_max, beta, face, profile_span, grid, &
                                     field, block_id, start_idx, cells, &
                                     gamma, alpha, kappa)
   character(len=*), intent(in) :: pml_type
   real(R8P),        intent(in) :: width
   real(R8P),        intent(in) :: gamma_max
   real(R8P),        intent(in) :: gamma_exponent
   real(R8P),        intent(in) :: alpha_max
   real(R8P),        intent(in) :: k_max
   real(R8P),        intent(in) :: beta
    integer(I4P),    intent(in) :: face
   real(R8P),        intent(in) :: profile_span
   type(grid_object), intent(in) :: grid
   type(field_object), intent(in) :: field
   integer(I4P),     intent(in) :: block_id
   integer(I4P),     intent(in) :: start_idx
   integer(I4P),     intent(in) :: cells
   real(R8P),    intent(out) :: gamma(1:)
   real(R8P),    intent(out) :: alpha(1:)
   real(R8P),    intent(out) :: kappa(1:)
   real(R8P)                 :: center_distance
   real(R8P)                 :: depth
   real(R8P)                 :: distance_to_outer
   real(R8P)                 :: ds
   real(R8P)                 :: distance0
   integer(I4P)              :: idx
   integer(I4P)              :: li

   gamma = 0._R8P
   alpha = 0._R8P
   kappa = 1._R8P
   if (cells <= 0_I4P) return

   do li = 1, min(cells, size(gamma))
      select case (face)
      case (PML_FACE_X_M)
         ds = field%dxyz(1,block_id)
         distance0 = field%emin(1,block_id) - grid%domain_emin(1)
         idx = start_idx + li - 1_I4P
         center_distance = distance0 + real(idx - 1_I4P, R8P) * ds
      case (PML_FACE_X_P)
         ds = field%dxyz(1,block_id)
         distance0 = grid%domain_emax(1) - field%emax(1,block_id)
         idx = start_idx + li - 1_I4P
         center_distance = distance0 + real(grid%ni - idx, R8P) * ds
      case (PML_FACE_Y_M)
         ds = field%dxyz(2,block_id)
         distance0 = field%emin(2,block_id) - grid%domain_emin(2)
         idx = start_idx + li - 1_I4P
         center_distance = distance0 + real(idx - 1_I4P, R8P) * ds
      case (PML_FACE_Y_P)
         ds = field%dxyz(2,block_id)
         distance0 = grid%domain_emax(2) - field%emax(2,block_id)
         idx = start_idx + li - 1_I4P
         center_distance = distance0 + real(grid%nj - idx, R8P) * ds
      case (PML_FACE_Z_M)
         ds = field%dxyz(3,block_id)
         distance0 = field%emin(3,block_id) - grid%domain_emin(3)
         idx = start_idx + li - 1_I4P
         center_distance = distance0 + real(idx - 1_I4P, R8P) * ds
      case default
         ds = field%dxyz(3,block_id)
         distance0 = grid%domain_emax(3) - field%emax(3,block_id)
         idx = start_idx + li - 1_I4P
         center_distance = distance0 + real(grid%nk - idx, R8P) * ds
      end select
      select case (trim(pml_type))
      case (PML_TYPE_CLASSIC)
         if (profile_span > 0._R8P) then
            depth = 1._R8P - center_distance / profile_span
         else
            depth = 1._R8P
         endif
         depth = max(0._R8P, min(1._R8P, depth))
         gamma(li) = gamma_max * depth**gamma_exponent
      case (PML_TYPE_BERMUDEZ)
         if (profile_span > 0._R8P) then
            depth = 1._R8P - center_distance / profile_span
            distance_to_outer = width * center_distance / profile_span + BERMUDEZ_EPS
         else
            depth = 1._R8P
            distance_to_outer = BERMUDEZ_EPS
         endif
         depth = max(0._R8P, min(1._R8P, depth))
         gamma(li) = beta / distance_to_outer**gamma_exponent
      case (PML_TYPE_CFS)
         if (profile_span > 0._R8P) then
            depth = 1._R8P - center_distance / profile_span
         else
            depth = 1._R8P
         endif
         depth = max(0._R8P, min(1._R8P, depth))
         gamma(li) = gamma_max * depth**CFS_PROFILE_EXPONENT
         alpha(li) = alpha_max * (1._R8P - depth)**CFS_PROFILE_EXPONENT
         kappa(li) = 1._R8P + (k_max - 1._R8P) * depth**CFS_PROFILE_EXPONENT
      end select
   enddo
   endsubroutine fill_face_coefficients

endmodule adam_prism_fnl_pml_object
