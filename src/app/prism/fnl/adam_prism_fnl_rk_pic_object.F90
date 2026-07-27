!< ADAM, PRISM PIC Runge-Kutta integrator, FNL backend.

#include "fundal.H"

module adam_prism_fnl_rk_pic_object
!< ADAM, PRISM PIC Runge-Kutta integrator, FNL backend.

use :: adam_prism_common_library
use :: adam_fnl_library
use :: fundal
use :: penf

implicit none
private
public :: prism_fnl_rk_pic_object

integer(I4P), parameter :: PIC_VARIABLES_NUMBER = 8_I4P

type :: prism_fnl_rk_pic_object
   real(R8P),    pointer :: alph_gpu(:,:)       => null() !< SSP RK alpha coefficients on device.
   real(R8P),    pointer :: beta_gpu(:)         => null() !< SSP RK beta coefficients on device.
   real(R8P),    pointer :: q_pic_rk_gpu(:,:,:) => null() !< RK stages [particle, variable, stage].
   integer(I4P)          :: particle_number = 0_I4P       !< Total number of particles.
   integer(I4P)          :: nrk             = 0_I4P       !< Number of stages.
contains
   procedure, pass(self) :: initialize
   procedure, pass(self) :: initialize_stages
   procedure, pass(self) :: compute_stage
   procedure, pass(self) :: assign_stage
   procedure, pass(self) :: update_q_pic
endtype prism_fnl_rk_pic_object

contains
   subroutine initialize(self, pic, rk_pic)
   !< Initialize device-side PIC RK state.
   class(prism_fnl_rk_pic_object), intent(inout) :: self   !< FNL PIC RK object.
   type(prism_pic_object),         intent(in)    :: pic    !< Host PIC object.
   type(prism_rk_pic_object),      intent(in)    :: rk_pic !< Host PIC RK object.
   integer(I4P)                                  :: ierr   !< Error status.

   self%particle_number = pic%particle_number
   self%nrk             = rk_pic%nrk

   if (self%particle_number == 0 .or. self%nrk == 0) return

   if (allocated(rk_pic%alph)) call dev_assign_to_device(src=rk_pic%alph, dst=self%alph_gpu)
   if (allocated(rk_pic%beta)) call dev_assign_to_device(src=rk_pic%beta, dst=self%beta_gpu)
   call dev_alloc(fptr_dev=self%q_pic_rk_gpu, ubounds=[self%particle_number, PIC_VARIABLES_NUMBER, self%nrk], &
                  lbounds=[1,1,1], init_value=0._R8P, ierr=ierr)
   endsubroutine initialize

   subroutine initialize_stages(self, q_pic_gpu)
   !< Initialize device SSP RK stages for PIC variables.
   class(prism_fnl_rk_pic_object), intent(inout) :: self             !< FNL PIC RK object.
   real(R8P),                      intent(in)    :: q_pic_gpu(1:,1:) !< PIC variables on device.
   integer(I4P)                                   :: n, s, v          !< Counters.
   real(R8P), pointer                             :: q_pic_rk_gpu(:,:,:) !< PIC stages on device.

   if (.not.associated(self%q_pic_rk_gpu)) return
   q_pic_rk_gpu => self%q_pic_rk_gpu

   !$acc parallel loop collapse(3) independent DEVICEVAR(q_pic_gpu, q_pic_rk_gpu)
   !$omp OMPLOOP collapse(3) DEVICEPTR(q_pic_gpu, q_pic_rk_gpu)
   do s = 1, self%nrk
      do v = 1, PIC_VARIABLES_NUMBER
         do n = 1, self%particle_number
            q_pic_rk_gpu(n,v,s) = q_pic_gpu(n,v)
         enddo
      enddo
   enddo
   endsubroutine initialize_stages

   subroutine compute_stage(self, s, dt)
   !< Compute one PIC SSP RK stage state on device.
   class(prism_fnl_rk_pic_object), intent(inout) :: self !< FNL PIC RK object.
   integer(I4P),                   intent(in)    :: s    !< Current stage number.
   real(R8P),                      intent(in)    :: dt   !< Time step.
   integer(I4P)                                   :: n, ss, v !< Counters.
   real(R8P), pointer                             :: q_pic_rk_gpu(:,:,:) !< PIC stages on device.
   real(R8P), pointer                             :: alph_gpu(:,:)       !< SSP RK alpha coefficients.

   if (.not.associated(self%q_pic_rk_gpu)) return
   q_pic_rk_gpu => self%q_pic_rk_gpu
   alph_gpu => self%alph_gpu

   !$acc parallel loop collapse(2) independent DEVICEVAR(q_pic_rk_gpu, alph_gpu)
   !$omp OMPLOOP collapse(2) DEVICEPTR(q_pic_rk_gpu, alph_gpu)
   do v = 1, 6
      do n = 1, self%particle_number
         !$acc loop seq
         do ss = 1, s - 1
            q_pic_rk_gpu(n,v,s) = q_pic_rk_gpu(n,v,s) + dt * alph_gpu(s,ss) * q_pic_rk_gpu(n,v,ss)
         enddo
      enddo
   enddo
   endsubroutine compute_stage

   subroutine assign_stage(self, s, pic_fields_gpu)
   !< Convert one PIC SSP RK stage state into its right-hand side on device.
   class(prism_fnl_rk_pic_object), intent(inout) :: self                  !< FNL PIC RK object.
   integer(I4P),                   intent(in)    :: s                     !< Current stage number.
   real(R8P),                      intent(in)    :: pic_fields_gpu(1:,1:) !< Fields at particle locations on device.
   integer(I4P)                                   :: n                     !< Counter.
   real(R8P)                                      :: vx, vy, vz            !< Particle velocity.
   real(R8P)                                      :: ex, ey, ez            !< Electric field.
   real(R8P)                                      :: bx, by, bz            !< Magnetic field.
   real(R8P)                                      :: fx, fy, fz            !< Lorentz force.
   real(R8P)                                      :: charge, mass          !< Particle charge and mass.
   real(R8P), pointer                             :: q_pic_rk_gpu(:,:,:)   !< PIC stages on device.

   if (.not.associated(self%q_pic_rk_gpu)) return
   q_pic_rk_gpu => self%q_pic_rk_gpu

   !$acc parallel loop independent DEVICEVAR(q_pic_rk_gpu, pic_fields_gpu)&
   !$acc& private(vx, vy, vz, ex, ey, ez, bx, by, bz, fx, fy, fz, charge, mass)
   !$omp OMPLOOP DEVICEPTR(q_pic_rk_gpu, pic_fields_gpu) &
   !$omp& private(vx, vy, vz, ex, ey, ez, bx, by, bz, fx, fy, fz, charge, mass)
   do n = 1, self%particle_number
      vx = q_pic_rk_gpu(n,4,s)
      vy = q_pic_rk_gpu(n,5,s)
      vz = q_pic_rk_gpu(n,6,s)
      charge = q_pic_rk_gpu(n,7,s)
      mass = q_pic_rk_gpu(n,8,s)

      ex = pic_fields_gpu(n,1) / EPS0
      ey = pic_fields_gpu(n,2) / EPS0
      ez = pic_fields_gpu(n,3) / EPS0
      bx = pic_fields_gpu(n,4)
      by = pic_fields_gpu(n,5)
      bz = pic_fields_gpu(n,6)

      fx = charge * (ex + (vy * bz - vz * by))
      fy = charge * (ey + (vz * bx - vx * bz))
      fz = charge * (ez + (vx * by - vy * bx))

      q_pic_rk_gpu(n,1,s) = vx
      q_pic_rk_gpu(n,2,s) = vy
      q_pic_rk_gpu(n,3,s) = vz
      q_pic_rk_gpu(n,4,s) = fx / mass
      q_pic_rk_gpu(n,5,s) = fy / mass
      q_pic_rk_gpu(n,6,s) = fz / mass
      q_pic_rk_gpu(n,7,s) = 0._R8P
      q_pic_rk_gpu(n,8,s) = 0._R8P
   enddo
   endsubroutine assign_stage

   subroutine update_q_pic(self, dt, q_pic_gpu)
   !< Complete the PIC SSP RK update on device.
   class(prism_fnl_rk_pic_object), intent(inout) :: self             !< FNL PIC RK object.
   real(R8P),                      intent(in)    :: dt               !< Time step.
   real(R8P),                      intent(inout) :: q_pic_gpu(1:,1:) !< PIC variables on device.
   integer(I4P)                                   :: n, s, v          !< Counters.
   real(R8P)                                      :: increment        !< RK increment accumulator.
   real(R8P), pointer                             :: q_pic_rk_gpu(:,:,:) !< PIC stages on device.
   real(R8P), pointer                             :: beta_gpu(:)         !< SSP RK beta coefficients.

   if (.not.associated(self%q_pic_rk_gpu)) return
   q_pic_rk_gpu => self%q_pic_rk_gpu
   beta_gpu => self%beta_gpu

   !$acc parallel loop collapse(2) independent DEVICEVAR(q_pic_gpu, q_pic_rk_gpu, beta_gpu)&
   !$acc& private(increment)
   !$omp OMPLOOP collapse(2) DEVICEPTR(q_pic_gpu, q_pic_rk_gpu, beta_gpu) &
   !$omp& private(increment)
   do v = 1, 6
      do n = 1, self%particle_number
         increment = 0._R8P
         !$acc loop seq
         do s = 1, self%nrk
            increment = increment + beta_gpu(s) * q_pic_rk_gpu(n,v,s)
         enddo
         q_pic_gpu(n,v) = q_pic_gpu(n,v) + dt * increment
      enddo
   enddo
   endsubroutine update_q_pic
endmodule adam_prism_fnl_rk_pic_object
