!< PRISM, RK support for PML auxiliary variables.
module adam_prism_rk_pml_object

use :: adam_rk_object,        only : rk_object, RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54
use :: adam_mpih_global,      only : mpih
use :: adam_prism_pml_object, only : prism_pml_object
use :: penf

implicit none
private
public :: prism_rk_pml_object

type :: prism_rk_pml_object
   logical                   :: enabled = .false. !< True only when a PML is active.
   character(:), allocatable :: scheme             !< RK scheme mirrored from the field RK object.
   integer(I4P)              :: nrk = 0_I4P        !< Number of SSP stages.
   real(R8P), allocatable    :: alph(:,:)          !< SSP alpha coefficients.
   real(R8P), allocatable    :: beta(:)            !< SSP beta coefficients.
   real(R8P), allocatable    :: gamm(:)            !< SSP gamma coefficients.
   real(R8P), allocatable    :: q_pml_x_m_rk(:,:,:,:,:,:) !< x-minus stage storage.
   real(R8P), allocatable    :: q_pml_x_p_rk(:,:,:,:,:,:) !< x-plus stage storage.
   real(R8P), allocatable    :: q_pml_y_m_rk(:,:,:,:,:,:) !< y-minus stage storage.
   real(R8P), allocatable    :: q_pml_y_p_rk(:,:,:,:,:,:) !< y-plus stage storage.
   real(R8P), allocatable    :: q_pml_z_m_rk(:,:,:,:,:,:) !< z-minus stage storage.
   real(R8P), allocatable    :: q_pml_z_p_rk(:,:,:,:,:,:) !< z-plus stage storage.
   real(R8P), allocatable    :: dq_pml_x_m(:,:,:,:,:)      !< x-minus current RHS.
   real(R8P), allocatable    :: dq_pml_x_p(:,:,:,:,:)      !< x-plus current RHS.
   real(R8P), allocatable    :: dq_pml_y_m(:,:,:,:,:)      !< y-minus current RHS.
   real(R8P), allocatable    :: dq_pml_y_p(:,:,:,:,:)      !< y-plus current RHS.
   real(R8P), allocatable    :: dq_pml_z_m(:,:,:,:,:)      !< z-minus current RHS.
   real(R8P), allocatable    :: dq_pml_z_p(:,:,:,:,:)      !< z-plus current RHS.
contains
   procedure, pass(self) :: assign_stage
   procedure, pass(self) :: compute_stage
   procedure, pass(self) :: initialize
   procedure, pass(self) :: initialize_stages
   procedure, pass(self) :: reset_rhs
   procedure, pass(self) :: update_q_pml
endtype prism_rk_pml_object

contains

   subroutine initialize(self, rk, pml)
   !< Mirror the SSP RK metadata and allocate stage storage for active PML faces.
   class(prism_rk_pml_object), intent(inout) :: self
   type(rk_object),            intent(in)    :: rk
   type(prism_pml_object),     intent(in)    :: pml

   call reset_rk_pml_object(self=self)
   if (.not. pml%enabled) return

   select case (trim(rk%scheme))
   case (RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54)
      continue
   case default
      call mpih%error_stop(msg=': PML auxiliary RK storage is available only for SSP RK schemes, found "'// &
                               trim(rk%scheme)//'"')
   endselect

   self%enabled = .true.
   self%scheme  = rk%scheme
   self%nrk     = rk%nrk

   allocate(self%alph(1:self%nrk,1:self%nrk))
   allocate(self%beta(1:self%nrk))
   allocate(self%gamm(1:self%nrk))
   self%alph = rk%alph
   self%beta = rk%beta
   self%gamm = rk%gamm

   call allocate_face_buffers(q_face=pml%q_pml_x_m, nrk=self%nrk, q_face_rk=self%q_pml_x_m_rk, dq_face=self%dq_pml_x_m)
   call allocate_face_buffers(q_face=pml%q_pml_x_p, nrk=self%nrk, q_face_rk=self%q_pml_x_p_rk, dq_face=self%dq_pml_x_p)
   call allocate_face_buffers(q_face=pml%q_pml_y_m, nrk=self%nrk, q_face_rk=self%q_pml_y_m_rk, dq_face=self%dq_pml_y_m)
   call allocate_face_buffers(q_face=pml%q_pml_y_p, nrk=self%nrk, q_face_rk=self%q_pml_y_p_rk, dq_face=self%dq_pml_y_p)
   call allocate_face_buffers(q_face=pml%q_pml_z_m, nrk=self%nrk, q_face_rk=self%q_pml_z_m_rk, dq_face=self%dq_pml_z_m)
   call allocate_face_buffers(q_face=pml%q_pml_z_p, nrk=self%nrk, q_face_rk=self%q_pml_z_p_rk, dq_face=self%dq_pml_z_p)
   endsubroutine initialize

   subroutine initialize_stages(self, pml)
   !< Replicate the committed PML state into every SSP stage.
   class(prism_rk_pml_object), intent(inout) :: self
   type(prism_pml_object),     intent(in)    :: pml

   if (.not. self%enabled) return
   call copy_face_to_stages(q_face=pml%q_pml_x_m, q_face_rk=self%q_pml_x_m_rk)
   call copy_face_to_stages(q_face=pml%q_pml_x_p, q_face_rk=self%q_pml_x_p_rk)
   call copy_face_to_stages(q_face=pml%q_pml_y_m, q_face_rk=self%q_pml_y_m_rk)
   call copy_face_to_stages(q_face=pml%q_pml_y_p, q_face_rk=self%q_pml_y_p_rk)
   call copy_face_to_stages(q_face=pml%q_pml_z_m, q_face_rk=self%q_pml_z_m_rk)
   call copy_face_to_stages(q_face=pml%q_pml_z_p, q_face_rk=self%q_pml_z_p_rk)
   call self%reset_rhs()
   endsubroutine initialize_stages

   subroutine reset_rhs(self)
   !< Zero the current PML auxiliary residual storage.
   class(prism_rk_pml_object), intent(inout) :: self

   if (.not. self%enabled) return
   call zero_face_rhs(dq_face=self%dq_pml_x_m)
   call zero_face_rhs(dq_face=self%dq_pml_x_p)
   call zero_face_rhs(dq_face=self%dq_pml_y_m)
   call zero_face_rhs(dq_face=self%dq_pml_y_p)
   call zero_face_rhs(dq_face=self%dq_pml_z_m)
   call zero_face_rhs(dq_face=self%dq_pml_z_p)
   endsubroutine reset_rhs

   subroutine compute_stage(self, s, dt)
   !< Assemble the SSP stage state for each active face-local PML buffer.
   class(prism_rk_pml_object), intent(inout) :: self
   integer(I4P),               intent(in)    :: s
   real(R8P),                  intent(in)    :: dt

   if (.not. self%enabled) return
   if (s <= 1_I4P) return

   call compute_face_stage(q_face_rk=self%q_pml_x_m_rk, alph=self%alph, s=s, dt=dt)
   call compute_face_stage(q_face_rk=self%q_pml_x_p_rk, alph=self%alph, s=s, dt=dt)
   call compute_face_stage(q_face_rk=self%q_pml_y_m_rk, alph=self%alph, s=s, dt=dt)
   call compute_face_stage(q_face_rk=self%q_pml_y_p_rk, alph=self%alph, s=s, dt=dt)
   call compute_face_stage(q_face_rk=self%q_pml_z_m_rk, alph=self%alph, s=s, dt=dt)
   call compute_face_stage(q_face_rk=self%q_pml_z_p_rk, alph=self%alph, s=s, dt=dt)
   endsubroutine compute_stage

   subroutine assign_stage(self, s)
   !< Store the current PML auxiliary residuals into SSP stage `s`.
   class(prism_rk_pml_object), intent(inout) :: self
   integer(I4P),               intent(in)    :: s

   if (.not. self%enabled) return

   call assign_face_stage(q_face_rk=self%q_pml_x_m_rk, dq_face=self%dq_pml_x_m, s=s)
   call assign_face_stage(q_face_rk=self%q_pml_x_p_rk, dq_face=self%dq_pml_x_p, s=s)
   call assign_face_stage(q_face_rk=self%q_pml_y_m_rk, dq_face=self%dq_pml_y_m, s=s)
   call assign_face_stage(q_face_rk=self%q_pml_y_p_rk, dq_face=self%dq_pml_y_p, s=s)
   call assign_face_stage(q_face_rk=self%q_pml_z_m_rk, dq_face=self%dq_pml_z_m, s=s)
   call assign_face_stage(q_face_rk=self%q_pml_z_p_rk, dq_face=self%dq_pml_z_p, s=s)
   endsubroutine assign_stage

   subroutine update_q_pml(self, dt, pml)
   !< Commit the SSP update to the reduced PML buffers.
   class(prism_rk_pml_object), intent(in)    :: self
   real(R8P),                  intent(in)    :: dt
   type(prism_pml_object),     intent(inout) :: pml

   if (.not. self%enabled) return

   call update_face_q(beta=self%beta, dt=dt, q_face_rk=self%q_pml_x_m_rk, q_face=pml%q_pml_x_m)
   call update_face_q(beta=self%beta, dt=dt, q_face_rk=self%q_pml_x_p_rk, q_face=pml%q_pml_x_p)
   call update_face_q(beta=self%beta, dt=dt, q_face_rk=self%q_pml_y_m_rk, q_face=pml%q_pml_y_m)
   call update_face_q(beta=self%beta, dt=dt, q_face_rk=self%q_pml_y_p_rk, q_face=pml%q_pml_y_p)
   call update_face_q(beta=self%beta, dt=dt, q_face_rk=self%q_pml_z_m_rk, q_face=pml%q_pml_z_m)
   call update_face_q(beta=self%beta, dt=dt, q_face_rk=self%q_pml_z_p_rk, q_face=pml%q_pml_z_p)
   endsubroutine update_q_pml

   subroutine allocate_face_buffers(q_face, nrk, q_face_rk, dq_face)
   !< Allocate one face-local stage buffer and its current RHS.
   real(R8P), allocatable, intent(in)    :: q_face(:,:,:,:,:)
   integer(I4P),           intent(in)    :: nrk
   real(R8P), allocatable, intent(inout) :: q_face_rk(:,:,:,:,:,:)
   real(R8P), allocatable, intent(inout) :: dq_face(:,:,:,:,:)

   if (.not. allocated(q_face)) return

   allocate(q_face_rk(1:size(q_face,1), 1:size(q_face,2), 1:size(q_face,3), 1:size(q_face,4), 1:size(q_face,5), 1:nrk))
   allocate(dq_face(1:size(q_face,1), 1:size(q_face,2), 1:size(q_face,3), 1:size(q_face,4), 1:size(q_face,5)))
   q_face_rk = 0._R8P
   dq_face   = 0._R8P
   endsubroutine allocate_face_buffers

   subroutine copy_face_to_stages(q_face, q_face_rk)
   !< Copy the committed state into every stage slot.
   real(R8P), allocatable, intent(in)    :: q_face(:,:,:,:,:)
   real(R8P), allocatable, intent(inout) :: q_face_rk(:,:,:,:,:,:)
   integer(I4P)                          :: i1, i2, i3, i4, i5, s

   if (.not. allocated(q_face_rk)) return
   !$omp parallel do collapse(6) default(firstprivate) shared(q_face, q_face_rk)
   do s=1, ubound(q_face_rk, dim=6)
      do i5=1, ubound(q_face_rk, dim=5)
         do i4=1, ubound(q_face_rk, dim=4)
            do i3=1, ubound(q_face_rk, dim=3)
               do i2=1, ubound(q_face_rk, dim=2)
                  do i1=1, ubound(q_face_rk, dim=1)
                     q_face_rk(i1,i2,i3,i4,i5,s) = q_face(i1,i2,i3,i4,i5)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine copy_face_to_stages

   subroutine compute_face_stage(q_face_rk, alph, s, dt)
   !< Apply the SSP linear combination for one face-local stage array.
   real(R8P), allocatable, intent(inout) :: q_face_rk(:,:,:,:,:,:)
   real(R8P),              intent(in)    :: alph(1:,1:)
   integer(I4P),           intent(in)    :: s
   real(R8P),              intent(in)    :: dt
   integer(I4P)                          :: i1, i2, i3, i4, i5, ss

   if (.not. allocated(q_face_rk)) return
   ! The stage accumulation over `ss` is a per-point reduction into
   ! q_face_rk(:,:,:,:,:,s), so `ss` must stay the inner sequential loop.
   !$omp parallel do collapse(5) default(firstprivate) shared(q_face_rk, alph, s, dt)
   do i5=1, ubound(q_face_rk, dim=5)
      do i4=1, ubound(q_face_rk, dim=4)
         do i3=1, ubound(q_face_rk, dim=3)
            do i2=1, ubound(q_face_rk, dim=2)
               do i1=1, ubound(q_face_rk, dim=1)
                  do ss=1, s-1
                     q_face_rk(i1,i2,i3,i4,i5,s) = q_face_rk(i1,i2,i3,i4,i5,s) + &
                                                    dt * alph(s,ss) * q_face_rk(i1,i2,i3,i4,i5,ss)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine compute_face_stage

   subroutine assign_face_stage(q_face_rk, dq_face, s)
   !< Store one face-local RHS into stage slot `s`.
   real(R8P), allocatable, intent(inout) :: q_face_rk(:,:,:,:,:,:)
   real(R8P), allocatable, intent(in)    :: dq_face(:,:,:,:,:)
   integer(I4P),           intent(in)    :: s
   integer(I4P)                          :: i1, i2, i3, i4, i5

   if (.not. allocated(q_face_rk)) return
   !$omp parallel do collapse(5) default(firstprivate) shared(q_face_rk, dq_face, s)
   do i5=1, ubound(q_face_rk, dim=5)
      do i4=1, ubound(q_face_rk, dim=4)
         do i3=1, ubound(q_face_rk, dim=3)
            do i2=1, ubound(q_face_rk, dim=2)
               do i1=1, ubound(q_face_rk, dim=1)
                  q_face_rk(i1,i2,i3,i4,i5,s) = dq_face(i1,i2,i3,i4,i5)
               enddo
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine assign_face_stage

   subroutine update_face_q(beta, dt, q_face_rk, q_face)
   !< Accumulate the SSP combination back into one face-local committed state.
   real(R8P),              intent(in)    :: beta(1:)
   real(R8P),              intent(in)    :: dt
   real(R8P), allocatable, intent(in)    :: q_face_rk(:,:,:,:,:,:)
   real(R8P), allocatable, intent(inout) :: q_face(:,:,:,:,:)
   integer(I4P)                          :: i1, i2, i3, i4, i5, s

   if (.not. allocated(q_face_rk)) return
   ! The final SSP combination over `s` is a per-point reduction into q_face,
   ! so the stage loop must stay the inner sequential loop.
   !$omp parallel do collapse(5) default(firstprivate) shared(q_face, q_face_rk, beta, dt)
   do i5=1, ubound(q_face, dim=5)
      do i4=1, ubound(q_face, dim=4)
         do i3=1, ubound(q_face, dim=3)
            do i2=1, ubound(q_face, dim=2)
               do i1=1, ubound(q_face, dim=1)
                  do s=1, ubound(q_face_rk, dim=6)
                     q_face(i1,i2,i3,i4,i5) = q_face(i1,i2,i3,i4,i5) + dt * beta(s) * q_face_rk(i1,i2,i3,i4,i5,s)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine update_face_q

   subroutine zero_face_rhs(dq_face)
   !< Zero one face-local RHS buffer with OpenMP over the face-local extents.
   real(R8P), allocatable, intent(inout) :: dq_face(:,:,:,:,:)
   integer(I4P)                          :: i1, i2, i3, i4, i5

   if (.not. allocated(dq_face)) return
   !$omp parallel do collapse(5) default(firstprivate) shared(dq_face)
   do i5=1, ubound(dq_face, dim=5)
      do i4=1, ubound(dq_face, dim=4)
         do i3=1, ubound(dq_face, dim=3)
            do i2=1, ubound(dq_face, dim=2)
               do i1=1, ubound(dq_face, dim=1)
                  dq_face(i1,i2,i3,i4,i5) = 0._R8P
               enddo
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine zero_face_rhs

   subroutine reset_rk_pml_object(self)
   !< Release all RK-PML metadata and storage.
   class(prism_rk_pml_object), intent(inout) :: self

   self%enabled = .false.
   self%nrk     = 0_I4P
   if (allocated(self%scheme))      deallocate(self%scheme)
   if (allocated(self%alph))        deallocate(self%alph)
   if (allocated(self%beta))        deallocate(self%beta)
   if (allocated(self%gamm))        deallocate(self%gamm)
   if (allocated(self%q_pml_x_m_rk)) deallocate(self%q_pml_x_m_rk)
   if (allocated(self%q_pml_x_p_rk)) deallocate(self%q_pml_x_p_rk)
   if (allocated(self%q_pml_y_m_rk)) deallocate(self%q_pml_y_m_rk)
   if (allocated(self%q_pml_y_p_rk)) deallocate(self%q_pml_y_p_rk)
   if (allocated(self%q_pml_z_m_rk)) deallocate(self%q_pml_z_m_rk)
   if (allocated(self%q_pml_z_p_rk)) deallocate(self%q_pml_z_p_rk)
   if (allocated(self%dq_pml_x_m))    deallocate(self%dq_pml_x_m)
   if (allocated(self%dq_pml_x_p))    deallocate(self%dq_pml_x_p)
   if (allocated(self%dq_pml_y_m))    deallocate(self%dq_pml_y_m)
   if (allocated(self%dq_pml_y_p))    deallocate(self%dq_pml_y_p)
   if (allocated(self%dq_pml_z_m))    deallocate(self%dq_pml_z_m)
   if (allocated(self%dq_pml_z_p))    deallocate(self%dq_pml_z_p)
   endsubroutine reset_rk_pml_object

endmodule adam_prism_rk_pml_object
