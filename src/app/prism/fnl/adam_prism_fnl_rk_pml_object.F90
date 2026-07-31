!< PRISM, device-side RK support for PML auxiliary variables.

#include "fundal.H"

module adam_prism_fnl_rk_pml_object

use :: adam_rk_object,            only : rk_object
use :: adam_prism_fnl_pml_object, only : prism_fnl_pml_object
use :: adam_fnl_mpih_global,      only : mpih_fnl
use :: fundal
use :: penf

implicit none
private
public :: prism_fnl_rk_pml_object

type :: prism_fnl_rk_pml_object
   logical                :: enabled = .false.
   integer(I4P)           :: nrk = 0_I4P
   real(R8P), pointer     :: alph_gpu(:,:) => null()
   real(R8P), pointer     :: beta_gpu(:)   => null()
   real(R8P), pointer     :: gamm_gpu(:)   => null()
   real(R8P), pointer     :: q_pml_x_m_rk_gpu(:,:,:,:,:,:) => null()
   real(R8P), pointer     :: q_pml_x_p_rk_gpu(:,:,:,:,:,:) => null()
   real(R8P), pointer     :: q_pml_y_m_rk_gpu(:,:,:,:,:,:) => null()
   real(R8P), pointer     :: q_pml_y_p_rk_gpu(:,:,:,:,:,:) => null()
   real(R8P), pointer     :: q_pml_z_m_rk_gpu(:,:,:,:,:,:) => null()
   real(R8P), pointer     :: q_pml_z_p_rk_gpu(:,:,:,:,:,:) => null()
   real(R8P), pointer     :: dq_pml_x_m_gpu(:,:,:,:,:) => null()
   real(R8P), pointer     :: dq_pml_x_p_gpu(:,:,:,:,:) => null()
   real(R8P), pointer     :: dq_pml_y_m_gpu(:,:,:,:,:) => null()
   real(R8P), pointer     :: dq_pml_y_p_gpu(:,:,:,:,:) => null()
   real(R8P), pointer     :: dq_pml_z_m_gpu(:,:,:,:,:) => null()
   real(R8P), pointer     :: dq_pml_z_p_gpu(:,:,:,:,:) => null()
contains
   procedure, pass(self) :: initialize
   procedure, pass(self) :: initialize_stages
   procedure, pass(self) :: reset_rhs
   procedure, pass(self) :: compute_stage
   procedure, pass(self) :: assign_stage
   procedure, pass(self) :: update_q_pml
endtype prism_fnl_rk_pml_object

contains

   subroutine initialize(self, rk, pml_fnl)
   class(prism_fnl_rk_pml_object), intent(inout) :: self
   type(rk_object),                intent(in)    :: rk
   type(prism_fnl_pml_object),     intent(in)    :: pml_fnl
   integer(I4P)                                  :: ierr

   if (.not. pml_fnl%enabled) return
   self%enabled = .true.
   self%nrk = rk%nrk
   if (allocated(rk%alph)) call dev_assign_to_device(src=rk%alph, dst=self%alph_gpu)
   if (allocated(rk%beta)) call dev_assign_to_device(src=rk%beta, dst=self%beta_gpu)
   if (allocated(rk%gamm)) call dev_assign_to_device(src=rk%gamm, dst=self%gamm_gpu)

   call allocate_face_buffers(q_face_gpu=pml_fnl%q_pml_x_m_gpu, nrk=self%nrk, q_face_rk_gpu=self%q_pml_x_m_rk_gpu, &
                              dq_face_gpu=self%dq_pml_x_m_gpu, label='x-')
   call allocate_face_buffers(q_face_gpu=pml_fnl%q_pml_x_p_gpu, nrk=self%nrk, q_face_rk_gpu=self%q_pml_x_p_rk_gpu, &
                              dq_face_gpu=self%dq_pml_x_p_gpu, label='x+')
   call allocate_face_buffers(q_face_gpu=pml_fnl%q_pml_y_m_gpu, nrk=self%nrk, q_face_rk_gpu=self%q_pml_y_m_rk_gpu, &
                              dq_face_gpu=self%dq_pml_y_m_gpu, label='y-')
   call allocate_face_buffers(q_face_gpu=pml_fnl%q_pml_y_p_gpu, nrk=self%nrk, q_face_rk_gpu=self%q_pml_y_p_rk_gpu, &
                              dq_face_gpu=self%dq_pml_y_p_gpu, label='y+')
   call allocate_face_buffers(q_face_gpu=pml_fnl%q_pml_z_m_gpu, nrk=self%nrk, q_face_rk_gpu=self%q_pml_z_m_rk_gpu, &
                              dq_face_gpu=self%dq_pml_z_m_gpu, label='z-')
   call allocate_face_buffers(q_face_gpu=pml_fnl%q_pml_z_p_gpu, nrk=self%nrk, q_face_rk_gpu=self%q_pml_z_p_rk_gpu, &
                              dq_face_gpu=self%dq_pml_z_p_gpu, label='z+')
   if (ierr /= 0_I4P) continue
   endsubroutine initialize

   subroutine initialize_stages(self, pml_fnl)
   class(prism_fnl_rk_pml_object), intent(inout) :: self
   type(prism_fnl_pml_object),     intent(in)    :: pml_fnl

   if (.not. self%enabled) return
   call copy_face_to_stages(q_face_gpu=pml_fnl%q_pml_x_m_gpu, q_face_rk_gpu=self%q_pml_x_m_rk_gpu, nrk=self%nrk)
   call copy_face_to_stages(q_face_gpu=pml_fnl%q_pml_x_p_gpu, q_face_rk_gpu=self%q_pml_x_p_rk_gpu, nrk=self%nrk)
   call copy_face_to_stages(q_face_gpu=pml_fnl%q_pml_y_m_gpu, q_face_rk_gpu=self%q_pml_y_m_rk_gpu, nrk=self%nrk)
   call copy_face_to_stages(q_face_gpu=pml_fnl%q_pml_y_p_gpu, q_face_rk_gpu=self%q_pml_y_p_rk_gpu, nrk=self%nrk)
   call copy_face_to_stages(q_face_gpu=pml_fnl%q_pml_z_m_gpu, q_face_rk_gpu=self%q_pml_z_m_rk_gpu, nrk=self%nrk)
   call copy_face_to_stages(q_face_gpu=pml_fnl%q_pml_z_p_gpu, q_face_rk_gpu=self%q_pml_z_p_rk_gpu, nrk=self%nrk)
   call self%reset_rhs()
   endsubroutine initialize_stages

   subroutine reset_rhs(self)
   class(prism_fnl_rk_pml_object), intent(inout) :: self

   if (.not. self%enabled) return
   call zero_face_rhs(dq_face_gpu=self%dq_pml_x_m_gpu)
   call zero_face_rhs(dq_face_gpu=self%dq_pml_x_p_gpu)
   call zero_face_rhs(dq_face_gpu=self%dq_pml_y_m_gpu)
   call zero_face_rhs(dq_face_gpu=self%dq_pml_y_p_gpu)
   call zero_face_rhs(dq_face_gpu=self%dq_pml_z_m_gpu)
   call zero_face_rhs(dq_face_gpu=self%dq_pml_z_p_gpu)
   endsubroutine reset_rhs

   subroutine compute_stage(self, s, dt)
   class(prism_fnl_rk_pml_object), intent(inout) :: self
   integer(I4P),                   intent(in)    :: s
   real(R8P),                      intent(in)    :: dt

   if (.not. self%enabled) return
   if (s <= 1_I4P) return
   call compute_face_stage(q_face_rk_gpu=self%q_pml_x_m_rk_gpu, alph_gpu=self%alph_gpu, s=s, dt=dt)
   call compute_face_stage(q_face_rk_gpu=self%q_pml_x_p_rk_gpu, alph_gpu=self%alph_gpu, s=s, dt=dt)
   call compute_face_stage(q_face_rk_gpu=self%q_pml_y_m_rk_gpu, alph_gpu=self%alph_gpu, s=s, dt=dt)
   call compute_face_stage(q_face_rk_gpu=self%q_pml_y_p_rk_gpu, alph_gpu=self%alph_gpu, s=s, dt=dt)
   call compute_face_stage(q_face_rk_gpu=self%q_pml_z_m_rk_gpu, alph_gpu=self%alph_gpu, s=s, dt=dt)
   call compute_face_stage(q_face_rk_gpu=self%q_pml_z_p_rk_gpu, alph_gpu=self%alph_gpu, s=s, dt=dt)
   endsubroutine compute_stage

   subroutine assign_stage(self, s)
   class(prism_fnl_rk_pml_object), intent(inout) :: self
   integer(I4P),                   intent(in)    :: s

   if (.not. self%enabled) return
   call assign_face_stage(q_face_rk_gpu=self%q_pml_x_m_rk_gpu, dq_face_gpu=self%dq_pml_x_m_gpu, s=s)
   call assign_face_stage(q_face_rk_gpu=self%q_pml_x_p_rk_gpu, dq_face_gpu=self%dq_pml_x_p_gpu, s=s)
   call assign_face_stage(q_face_rk_gpu=self%q_pml_y_m_rk_gpu, dq_face_gpu=self%dq_pml_y_m_gpu, s=s)
   call assign_face_stage(q_face_rk_gpu=self%q_pml_y_p_rk_gpu, dq_face_gpu=self%dq_pml_y_p_gpu, s=s)
   call assign_face_stage(q_face_rk_gpu=self%q_pml_z_m_rk_gpu, dq_face_gpu=self%dq_pml_z_m_gpu, s=s)
   call assign_face_stage(q_face_rk_gpu=self%q_pml_z_p_rk_gpu, dq_face_gpu=self%dq_pml_z_p_gpu, s=s)
   endsubroutine assign_stage

   subroutine update_q_pml(self, dt, pml_fnl)
   class(prism_fnl_rk_pml_object), intent(in)    :: self
   real(R8P),                      intent(in)    :: dt
   type(prism_fnl_pml_object),     intent(inout) :: pml_fnl

   if (.not. self%enabled) return
   call update_face_q(beta_gpu=self%beta_gpu, dt=dt, q_face_rk_gpu=self%q_pml_x_m_rk_gpu, q_face_gpu=pml_fnl%q_pml_x_m_gpu)
   call update_face_q(beta_gpu=self%beta_gpu, dt=dt, q_face_rk_gpu=self%q_pml_x_p_rk_gpu, q_face_gpu=pml_fnl%q_pml_x_p_gpu)
   call update_face_q(beta_gpu=self%beta_gpu, dt=dt, q_face_rk_gpu=self%q_pml_y_m_rk_gpu, q_face_gpu=pml_fnl%q_pml_y_m_gpu)
   call update_face_q(beta_gpu=self%beta_gpu, dt=dt, q_face_rk_gpu=self%q_pml_y_p_rk_gpu, q_face_gpu=pml_fnl%q_pml_y_p_gpu)
   call update_face_q(beta_gpu=self%beta_gpu, dt=dt, q_face_rk_gpu=self%q_pml_z_m_rk_gpu, q_face_gpu=pml_fnl%q_pml_z_m_gpu)
   call update_face_q(beta_gpu=self%beta_gpu, dt=dt, q_face_rk_gpu=self%q_pml_z_p_rk_gpu, q_face_gpu=pml_fnl%q_pml_z_p_gpu)
   endsubroutine update_q_pml

   subroutine allocate_face_buffers(q_face_gpu, nrk, q_face_rk_gpu, dq_face_gpu, label)
   real(R8P), pointer,   intent(in)    :: q_face_gpu(:,:,:,:,:)
   integer(I4P),         intent(in)    :: nrk
   real(R8P), pointer,   intent(inout) :: q_face_rk_gpu(:,:,:,:,:,:)
   real(R8P), pointer,   intent(inout) :: dq_face_gpu(:,:,:,:,:)
   character(*),         intent(in)    :: label
   integer(I4P)                        :: ierr

   if (.not. associated(q_face_gpu)) return
   call dev_alloc(fptr_dev=q_face_rk_gpu, ubounds=[ubound(q_face_gpu,1),ubound(q_face_gpu,2),ubound(q_face_gpu,3), &
                                                   ubound(q_face_gpu,4),ubound(q_face_gpu,5),nrk],                  &
                  lbounds=[1,1,1,1,1,1], init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate RK PML stages on face '//trim(label))
   call dev_alloc(fptr_dev=dq_face_gpu, ubounds=[ubound(q_face_gpu,1),ubound(q_face_gpu,2),ubound(q_face_gpu,3), &
                                                 ubound(q_face_gpu,4),ubound(q_face_gpu,5)],                      &
                  lbounds=[1,1,1,1,1], init_value=0._R8P, ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate RK PML RHS on face '//trim(label))
   endsubroutine allocate_face_buffers

   subroutine copy_face_to_stages(q_face_gpu, q_face_rk_gpu, nrk)
   real(R8P), pointer, intent(in)    :: q_face_gpu(:,:,:,:,:)
   real(R8P), pointer, intent(inout) :: q_face_rk_gpu(:,:,:,:,:,:)
   integer(I4P),       intent(in)    :: nrk
   integer(I4P)                      :: i1, i2, i3, i4, i5, s

   if (.not. associated(q_face_rk_gpu)) return
   !$acc parallel loop collapse(6) independent DEVICEVAR(q_face_gpu, q_face_rk_gpu)
   !$omp OMPLOOP collapse(6) DEVICEPTR(q_face_gpu, q_face_rk_gpu)
   do s = 1, nrk
      do i5 = 1, ubound(q_face_rk_gpu,5)
         do i4 = 1, ubound(q_face_rk_gpu,4)
            do i3 = 1, ubound(q_face_rk_gpu,3)
               do i2 = 1, ubound(q_face_rk_gpu,2)
                  do i1 = 1, ubound(q_face_rk_gpu,1)
                     q_face_rk_gpu(i1,i2,i3,i4,i5,s) = q_face_gpu(i1,i2,i3,i4,i5)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine copy_face_to_stages

   subroutine compute_face_stage(q_face_rk_gpu, alph_gpu, s, dt)
   real(R8P), pointer, intent(inout) :: q_face_rk_gpu(:,:,:,:,:,:)
   real(R8P), pointer, intent(in)    :: alph_gpu(:,:)
   integer(I4P),       intent(in)    :: s
   real(R8P),          intent(in)    :: dt
   integer(I4P)                      :: i1, i2, i3, i4, i5, ss

   if (.not. associated(q_face_rk_gpu)) return
   !$acc parallel loop collapse(5) independent DEVICEVAR(q_face_rk_gpu, alph_gpu)
   !$omp OMPLOOP collapse(5) DEVICEPTR(q_face_rk_gpu, alph_gpu)
   do i5 = 1, ubound(q_face_rk_gpu,5)
      do i4 = 1, ubound(q_face_rk_gpu,4)
         do i3 = 1, ubound(q_face_rk_gpu,3)
            do i2 = 1, ubound(q_face_rk_gpu,2)
               do i1 = 1, ubound(q_face_rk_gpu,1)
                  !$acc loop seq
                  do ss = 1, s - 1
                     q_face_rk_gpu(i1,i2,i3,i4,i5,s) = q_face_rk_gpu(i1,i2,i3,i4,i5,s) + &
                                                        dt * alph_gpu(s,ss) * q_face_rk_gpu(i1,i2,i3,i4,i5,ss)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine compute_face_stage

   subroutine assign_face_stage(q_face_rk_gpu, dq_face_gpu, s)
   real(R8P), pointer, intent(inout) :: q_face_rk_gpu(:,:,:,:,:,:)
   real(R8P), pointer, intent(in)    :: dq_face_gpu(:,:,:,:,:)
   integer(I4P),       intent(in)    :: s
   integer(I4P)                      :: i1, i2, i3, i4, i5

   if (.not. associated(q_face_rk_gpu)) return
   !$acc parallel loop collapse(5) independent DEVICEVAR(q_face_rk_gpu, dq_face_gpu)
   !$omp OMPLOOP collapse(5) DEVICEPTR(q_face_rk_gpu, dq_face_gpu)
   do i5 = 1, ubound(q_face_rk_gpu,5)
      do i4 = 1, ubound(q_face_rk_gpu,4)
         do i3 = 1, ubound(q_face_rk_gpu,3)
            do i2 = 1, ubound(q_face_rk_gpu,2)
               do i1 = 1, ubound(q_face_rk_gpu,1)
                  q_face_rk_gpu(i1,i2,i3,i4,i5,s) = dq_face_gpu(i1,i2,i3,i4,i5)
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine assign_face_stage

   subroutine update_face_q(beta_gpu, dt, q_face_rk_gpu, q_face_gpu)
   real(R8P), pointer, intent(in)    :: beta_gpu(:)
   real(R8P),          intent(in)    :: dt
   real(R8P), pointer, intent(in)    :: q_face_rk_gpu(:,:,:,:,:,:)
   real(R8P), pointer, intent(inout) :: q_face_gpu(:,:,:,:,:)
   integer(I4P)                      :: i1, i2, i3, i4, i5, s
   real(R8P)                         :: increment

   if (.not. associated(q_face_rk_gpu)) return
   !$acc parallel loop collapse(5) independent DEVICEVAR(beta_gpu, q_face_rk_gpu, q_face_gpu) private(increment)
   !$omp OMPLOOP collapse(5) DEVICEPTR(beta_gpu, q_face_rk_gpu, q_face_gpu) private(increment)
   do i5 = 1, ubound(q_face_gpu,5)
      do i4 = 1, ubound(q_face_gpu,4)
         do i3 = 1, ubound(q_face_gpu,3)
            do i2 = 1, ubound(q_face_gpu,2)
               do i1 = 1, ubound(q_face_gpu,1)
                  increment = 0._R8P
                  !$acc loop seq
                  do s = 1, ubound(q_face_rk_gpu,6)
                     increment = increment + beta_gpu(s) * q_face_rk_gpu(i1,i2,i3,i4,i5,s)
                  enddo
                  q_face_gpu(i1,i2,i3,i4,i5) = q_face_gpu(i1,i2,i3,i4,i5) + dt * increment
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine update_face_q

   subroutine zero_face_rhs(dq_face_gpu)
   real(R8P), pointer, intent(inout) :: dq_face_gpu(:,:,:,:,:)
   integer(I4P)                      :: i1, i2, i3, i4, i5

   if (.not. associated(dq_face_gpu)) return
   !$acc parallel loop collapse(5) independent DEVICEVAR(dq_face_gpu)
   !$omp OMPLOOP collapse(5) DEVICEPTR(dq_face_gpu)
   do i5 = 1, ubound(dq_face_gpu,5)
      do i4 = 1, ubound(dq_face_gpu,4)
         do i3 = 1, ubound(dq_face_gpu,3)
            do i2 = 1, ubound(dq_face_gpu,2)
               do i1 = 1, ubound(dq_face_gpu,1)
                  dq_face_gpu(i1,i2,i3,i4,i5) = 0._R8P
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine zero_face_rhs

endmodule adam_prism_fnl_rk_pml_object
