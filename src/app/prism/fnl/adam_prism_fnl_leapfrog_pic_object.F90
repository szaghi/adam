!< ADAM, PRISM PIC leapfrog integrator, FNL backend.

#include "fundal.H"

module adam_prism_fnl_leapfrog_pic_object
!< ADAM, PRISM PIC leapfrog integrator, FNL backend.

use :: adam_prism_common_library
use :: adam_fnl_library
use :: fundal
use :: penf

implicit none
private
public :: prism_fnl_leapfrog_pic_object

integer(I4P), parameter :: PIC_VARIABLES_NUMBER = 8_I4P

type :: prism_fnl_leapfrog_pic_object
   real(R8P),    pointer :: q_pic_old_gpu(:,:,:) => null() !< Leapfrog history [particle, variable, step].
   real(R8P)             :: nu = 0.01_R8P                    !< RAW filter coefficient.
   real(R8P)             :: alpha = 0.53_R8P                !< RAW-Williams coefficient.
   logical               :: is_filtered = .false.           !< Enable RAW filter.
   integer(I4P)          :: particle_number = 0_I4P         !< Total number of particles.
contains
   procedure, pass(self) :: initialize
   procedure, pass(self) :: assign_step
   procedure, pass(self) :: prime_step
   procedure, pass(self) :: integrate
endtype prism_fnl_leapfrog_pic_object

contains
   subroutine initialize(self, pic, leapfrog_pic)
   !< Initialize device-side PIC leapfrog state.
   class(prism_fnl_leapfrog_pic_object), intent(inout) :: self         !< FNL PIC leapfrog object.
   type(prism_pic_object),               intent(in)    :: pic          !< Host PIC object.
   type(prism_leapfrog_pic_object),      intent(in)    :: leapfrog_pic !< Host PIC leapfrog object.
   integer(I4P)                                        :: ierr         !< Error status.

   self%particle_number = pic%particle_number
   self%nu              = leapfrog_pic%nu
   self%alpha           = leapfrog_pic%alpha
   self%is_filtered     = leapfrog_pic%is_filtered

   if (self%particle_number == 0) return

   call dev_alloc(fptr_dev=self%q_pic_old_gpu, ubounds=[self%particle_number, PIC_VARIABLES_NUMBER, 2_I4P], &
                  lbounds=[1,1,1], init_value=0._R8P, ierr=ierr)
   endsubroutine initialize

   subroutine assign_step(self, s, q_pic_gpu)
   !< Store one PIC state in the device leapfrog history buffer.
   class(prism_fnl_leapfrog_pic_object), intent(inout) :: self             !< FNL PIC leapfrog object.
   integer(I4P),                         intent(in)    :: s                !< Leapfrog history slot.
   real(R8P),                            intent(in)    :: q_pic_gpu(1:,1:) !< PIC variables on device.
   integer(I4P)                                        :: n, v             !< Counters.
   real(R8P), pointer                                  :: q_pic_old_gpu(:,:,:) !< Leapfrog history on device.

   if (.not.associated(self%q_pic_old_gpu)) return
   q_pic_old_gpu => self%q_pic_old_gpu

   !$acc parallel loop collapse(2) independent DEVICEVAR(q_pic_gpu, q_pic_old_gpu)
   !$omp OMPLOOP collapse(2) DEVICEPTR(q_pic_gpu, q_pic_old_gpu)
   do v = 1, PIC_VARIABLES_NUMBER
      do n = 1, self%particle_number
         q_pic_old_gpu(n,v,s) = q_pic_gpu(n,v)
      enddo
   enddo
   endsubroutine assign_step

   subroutine prime_step(self, dt, q_pic_gpu, pic_fields_gpu)
   !< Perform the explicit Euler priming step for PIC leapfrog on device.
   class(prism_fnl_leapfrog_pic_object), intent(inout) :: self                  !< FNL PIC leapfrog object.
   real(R8P),                            intent(in)    :: dt                    !< Time step.
   real(R8P),                            intent(inout) :: q_pic_gpu(1:,1:)      !< PIC variables on device.
   real(R8P),                            intent(in)    :: pic_fields_gpu(1:,1:) !< Fields at particle locations on device.
   integer(I4P)                                        :: n                     !< Counter.
   real(R8P)                                           :: vx, vy, vz            !< Particle velocity.
   real(R8P)                                           :: fx, fy, fz            !< Lorentz force.
   real(R8P)                                           :: charge_mass_ratio     !< Charge-over-mass ratio.

   if (self%particle_number == 0) return

   !$acc parallel loop independent DEVICEVAR(q_pic_gpu, pic_fields_gpu)&
   !$acc& private(vx, vy, vz, fx, fy, fz, charge_mass_ratio)
   !$omp OMPLOOP DEVICEPTR(q_pic_gpu, pic_fields_gpu) &
   !$omp& private(vx, vy, vz, fx, fy, fz, charge_mass_ratio)
   do n = 1, self%particle_number
      vx = q_pic_gpu(n,4)
      vy = q_pic_gpu(n,5)
      vz = q_pic_gpu(n,6)

      q_pic_gpu(n,1) = q_pic_gpu(n,1) + dt * vx
      q_pic_gpu(n,2) = q_pic_gpu(n,2) + dt * vy
      q_pic_gpu(n,3) = q_pic_gpu(n,3) + dt * vz

      fx = (vy * pic_fields_gpu(n,6)) - (vz * pic_fields_gpu(n,5))
      fy = (vz * pic_fields_gpu(n,4)) - (vx * pic_fields_gpu(n,6))
      fz = (vx * pic_fields_gpu(n,5)) - (vy * pic_fields_gpu(n,4))
      charge_mass_ratio = q_pic_gpu(n,7) / q_pic_gpu(n,8)

      q_pic_gpu(n,4) = q_pic_gpu(n,4) + dt * charge_mass_ratio * (pic_fields_gpu(n,1) + fx)
      q_pic_gpu(n,5) = q_pic_gpu(n,5) + dt * charge_mass_ratio * (pic_fields_gpu(n,2) + fy)
      q_pic_gpu(n,6) = q_pic_gpu(n,6) + dt * charge_mass_ratio * (pic_fields_gpu(n,3) + fz)
   enddo
   endsubroutine prime_step

   subroutine integrate(self, dt, q_pic_gpu, pic_fields_gpu)
   !< Integrate PIC particles with leapfrog on device.
   class(prism_fnl_leapfrog_pic_object), intent(inout) :: self                  !< FNL PIC leapfrog object.
   real(R8P),                            intent(in)    :: dt                    !< Time step.
   real(R8P),                            intent(inout) :: q_pic_gpu(1:,1:)      !< PIC variables on device.
   real(R8P),                            intent(in)    :: pic_fields_gpu(1:,1:) !< Fields at particle locations on device.
   integer(I4P)                                        :: n, v                  !< Counters.
   real(R8P)                                           :: filter                !< RAW filter increment.
   real(R8P)                                           :: charge_mass_ratio     !< Charge-over-mass ratio.
   real(R8P)                                           :: t1, t2, t3            !< Buneman-Boris rotation vector.
   real(R8P)                                           :: s1, s2, s3            !< Buneman-Boris auxiliary vector.
   real(R8P)                                           :: w1, w2, w3            !< Buneman-Boris auxiliary vector.
   real(R8P)                                           :: v_star_1, v_star_2, v_star_3
   real(R8P)                                           :: v_star_star_1, v_star_star_2, v_star_star_3
   real(R8P)                                           :: t_norm2               !< Squared norm of t.
   real(R8P), pointer                                  :: q_pic_old_gpu(:,:,:)  !< Leapfrog history on device.

   if (.not.associated(self%q_pic_old_gpu)) return
   q_pic_old_gpu => self%q_pic_old_gpu

   !$acc parallel loop independent DEVICEVAR(q_pic_gpu, pic_fields_gpu, q_pic_old_gpu)&
   !$acc& private(charge_mass_ratio, t1, t2, t3, s1, s2, s3, w1, w2, w3, v_star_1, v_star_2, v_star_3, &
   !$acc&         v_star_star_1, v_star_star_2, v_star_star_3, t_norm2)
   !$omp OMPLOOP DEVICEPTR(q_pic_gpu, pic_fields_gpu, q_pic_old_gpu) &
   !$omp& private(charge_mass_ratio, t1, t2, t3, s1, s2, s3, w1, w2, w3, v_star_1, v_star_2, v_star_3, &
   !$omp&         v_star_star_1, v_star_star_2, v_star_star_3, t_norm2)
   do n = 1, self%particle_number
      charge_mass_ratio = q_pic_gpu(n,7) / q_pic_gpu(n,8)

      v_star_1 = q_pic_old_gpu(n,4,1) + dt * pic_fields_gpu(n,1) * charge_mass_ratio
      v_star_2 = q_pic_old_gpu(n,5,1) + dt * pic_fields_gpu(n,2) * charge_mass_ratio
      v_star_3 = q_pic_old_gpu(n,6,1) + dt * pic_fields_gpu(n,3) * charge_mass_ratio

      t1 = dt * charge_mass_ratio * pic_fields_gpu(n,4)
      t2 = dt * charge_mass_ratio * pic_fields_gpu(n,5)
      t3 = dt * charge_mass_ratio * pic_fields_gpu(n,6)

      w1 = v_star_1 + (v_star_2 * t3 - v_star_3 * t2)
      w2 = v_star_2 + (v_star_3 * t1 - v_star_1 * t3)
      w3 = v_star_3 + (v_star_1 * t2 - v_star_2 * t1)

      t_norm2 = t1 * t1 + t2 * t2 + t3 * t3
      s1 = 2._R8P * t1 / (1._R8P + t_norm2)
      s2 = 2._R8P * t2 / (1._R8P + t_norm2)
      s3 = 2._R8P * t3 / (1._R8P + t_norm2)

      v_star_star_1 = v_star_1 + (w2 * s3 - w3 * s2)
      v_star_star_2 = v_star_2 + (w3 * s1 - w1 * s3)
      v_star_star_3 = v_star_3 + (w1 * s2 - w2 * s1)

      q_pic_old_gpu(n,4,1) = q_pic_gpu(n,4)
      q_pic_old_gpu(n,5,1) = q_pic_gpu(n,5)
      q_pic_old_gpu(n,6,1) = q_pic_gpu(n,6)

      q_pic_old_gpu(n,4,2) = v_star_star_1 + dt * pic_fields_gpu(n,1) * charge_mass_ratio
      q_pic_old_gpu(n,5,2) = v_star_star_2 + dt * pic_fields_gpu(n,2) * charge_mass_ratio
      q_pic_old_gpu(n,6,2) = v_star_star_3 + dt * pic_fields_gpu(n,3) * charge_mass_ratio

      q_pic_gpu(n,4) = q_pic_old_gpu(n,4,2)
      q_pic_gpu(n,5) = q_pic_old_gpu(n,5,2)
      q_pic_gpu(n,6) = q_pic_old_gpu(n,6,2)

      do v = 1, 3
         q_pic_old_gpu(n,v,2) = q_pic_old_gpu(n,v,1) + 2._R8P * dt * q_pic_old_gpu(n,v+3,1)
         q_pic_old_gpu(n,v,1) = q_pic_gpu(n,v)
         q_pic_gpu(n,v) = q_pic_old_gpu(n,v,2)
      enddo
   enddo

   if (self%is_filtered) then
      !$acc parallel loop collapse(2) independent DEVICEVAR(q_pic_gpu, q_pic_old_gpu)&
      !$acc& private(filter)
      !$omp OMPLOOP collapse(2) DEVICEPTR(q_pic_gpu, q_pic_old_gpu) &
      !$omp& private(filter)
      do v = 1, PIC_VARIABLES_NUMBER
         do n = 1, self%particle_number
            filter = (q_pic_old_gpu(n,v,1) - 2._R8P * q_pic_old_gpu(n,v,2) + q_pic_gpu(n,v)) * self%nu * 0.5_R8P
            q_pic_old_gpu(n,v,2) = q_pic_old_gpu(n,v,2) + filter * self%alpha
            q_pic_gpu(n,v) = q_pic_gpu(n,v) + filter * (self%alpha - 1._R8P)
         enddo
      enddo
   endif
   endsubroutine integrate
endmodule adam_prism_fnl_leapfrog_pic_object
