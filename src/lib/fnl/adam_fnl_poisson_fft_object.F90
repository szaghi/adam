!< ADAM, FNL FFT Poisson solver on uniform periodic grids (cuFFT).
module adam_fnl_poisson_fft_object
!< ADAM, FNL FFT Poisson solver on uniform periodic grids (cuFFT).
!<
!< Extends `poisson_fft_object` (issue #81, P2): the gather to rank 0, the placement and the scatter of `phi` with its
!< ghosts are the CPU ones; the spectral solve runs on rank 0's GPU through cuFFT (`-D_CUFFT`, `-cudalib=cufft`). The
!< global source goes to the device, a D2Z transform, the division by the eigenvalues of `L_h` (the same host-computed
!< `lam1..3` as the CPU solver), a Z2D transform, and the potential comes back. cuFFT runs on the stream of the
!< synchronous OpenACC queue, so the transforms and the kernel between them are ordered.
!<
!< Isolated boundaries (issue #81, M8-P4): the transforms are the padded ones (`m`), the kernel's spectrum `khat` is
!< computed once on the device (`transform_kernel`) and stays there (`enter data`, released by `destroy_plans`); the
!< solve multiplies by it instead of dividing by the eigenvalues.
!<
!< cuFFT is NVIDIA only: built without `_CUFFT` (an AMD build), making the plans stops the run.
use :: adam_mpih_global, only : mpih
use :: adam_poisson_fft_object, only : poisson_fft_object
#ifdef _CUFFT
use :: cufft
use :: openacc
#endif
use :: iso_c_binding
use :: penf, only : I4P, R8P, str

implicit none
private
public :: poisson_fft_fnl_object

type, extends(poisson_fft_object) :: poisson_fft_fnl_object
   !< FNL FFT Poisson solver on uniform periodic grids (cuFFT).
   integer(c_int) :: plan_d2z=0_c_int  !< cuFFT forward plan (rank 0).
   integer(c_int) :: plan_z2d=0_c_int  !< cuFFT backward plan (rank 0).
   logical        :: has_plans=.false. !< The cuFFT plans exist.
   logical        :: khat_on_device=.false. !< `khat` is on the device (isolated).
   contains
      ! backend hooks
      procedure, pass(self) :: destroy_plans => destroy_plans_fnl !< Destroy the cuFFT plans.
      procedure, pass(self) :: make_plans    => make_plans_fnl    !< Make the cuFFT plans (rank 0).
      procedure, pass(self) :: transform     => transform_fnl     !< The spectral solve on the GPU (rank 0).
      procedure, pass(self) :: transform_kernel => transform_kernel_fnl !< The kernel's spectrum on the GPU (isolated).
endtype poisson_fft_fnl_object

contains
   ! backend hooks
   subroutine destroy_plans_fnl(self)
   !< Destroy the cuFFT plans.
   class(poisson_fft_fnl_object), intent(inout) :: self !< The solver.
   integer(c_int)                               :: ierr !< cuFFT status.

   ierr = 0_c_int
#ifdef _CUFFT
   if (self%has_plans) then
      ierr = cufftDestroy(self%plan_d2z)
      ierr = cufftDestroy(self%plan_z2d)
   endif
   if (self%khat_on_device) then
      !$acc exit data delete(self%khat)
   endif
#endif
   self%has_plans = .false.
   self%khat_on_device = .false.
   endsubroutine destroy_plans_fnl

   subroutine make_plans_fnl(self)
   !< Make the cuFFT plans (D2Z and Z2D, cuFFT is row-major: the slowest dimension first) on the synchronous queue.
   class(poisson_fft_fnl_object), intent(inout) :: self !< The solver.
   integer(c_int)                               :: ierr !< cuFFT status.

   ierr = 0_c_int
#ifdef _CUFFT
   ierr = cufftPlan3d(self%plan_d2z, int(self%m(3), c_int), int(self%m(2), c_int), int(self%m(1), c_int), CUFFT_D2Z)
   if (ierr /= CUFFT_SUCCESS) call mpih%error_stop(msg='poisson_fft_fnl_object%make_plans: cufftPlan3d D2Z failed, '// &
                                                       'status '//trim(str(int(ierr, I4P), .true.)))
   ierr = cufftPlan3d(self%plan_z2d, int(self%m(3), c_int), int(self%m(2), c_int), int(self%m(1), c_int), CUFFT_Z2D)
   if (ierr /= CUFFT_SUCCESS) call mpih%error_stop(msg='poisson_fft_fnl_object%make_plans: cufftPlan3d Z2D failed, '// &
                                                       'status '//trim(str(int(ierr, I4P), .true.)))
   ierr = cufftSetStream(self%plan_d2z, acc_get_cuda_stream(acc_async_sync))
   ierr = ierr + cufftSetStream(self%plan_z2d, acc_get_cuda_stream(acc_async_sync))
   if (ierr /= CUFFT_SUCCESS) call mpih%error_stop(msg='poisson_fft_fnl_object%make_plans: cufftSetStream failed')
   self%has_plans = .true.
#else
   call mpih%error_stop(msg='poisson_fft_fnl_object: built without cuFFT (build a *-fft mode with -D_CUFFT; cuFFT is '// &
                            'NVIDIA only)')
#endif
   endsubroutine make_plans_fnl

   subroutine transform_fnl(self, factor)
   !< Rank 0: `f_global <- factor * L_h^-1 f_global` on the GPU, the modes of zero eigenvalue (the mean) set to zero.
   class(poisson_fft_fnl_object), intent(inout) :: self   !< The solver.
   real(R8P),                     intent(in)    :: factor !< Source factor.

#ifdef _CUFFT
   if (self%isolated) then
      call spectral_convolve_dev(plan_d2z=self%plan_d2z, plan_z2d=self%plan_z2d, m1=self%m(1), m2=self%m(2), &
                                 m3=self%m(3), factor=factor, khat=self%khat, f=self%f_global, f_hat=self%f_hat)
      return
   endif
   call spectral_solve_dev(plan_d2z=self%plan_d2z, plan_z2d=self%plan_z2d,                     &
                           n1=self%nt(1), n2=self%nt(2), n3=self%nt(3),                         &
                           scale=factor / real(product(self%nt), R8P),                          &
                           lam1=self%lam1, lam2=self%lam2, lam3=self%lam3, f=self%f_global, f_hat=self%f_hat)
#else
   call mpih%error_stop(msg='poisson_fft_fnl_object%transform: built without cuFFT, factor '//trim(str(factor)))
#endif
   endsubroutine transform_fnl

   subroutine transform_kernel_fnl(self)
   !< Rank 0, isolated: `khat <- FFT(f_global) / prod(m)` on the GPU (real: the kernel is even); `khat` stays on the
   !< device for every solve.
   class(poisson_fft_fnl_object), intent(inout) :: self !< The solver.

#ifdef _CUFFT
   call kernel_spectrum_dev(plan_d2z=self%plan_d2z, m1=self%m(1), m2=self%m(2), m3=self%m(3), f=self%f_global, &
                            f_hat=self%f_hat)
   self%khat = real(self%f_hat, R8P) / real(product(self%m), R8P)
   !$acc enter data copyin(self%khat)
   self%khat_on_device = .true.
#else
   call mpih%error_stop(msg='poisson_fft_fnl_object%transform_kernel: built without cuFFT')
#endif
   endsubroutine transform_kernel_fnl

#ifdef _CUFFT
   ! private procedures
   subroutine spectral_solve_dev(plan_d2z, plan_z2d, n1, n2, n3, scale, lam1, lam2, lam3, f, f_hat)
   !< The spectral solve on the device: D2Z, division by the eigenvalues, Z2D; f in, the potential out (in place).
   integer(c_int), intent(in)    :: plan_d2z            !< cuFFT forward plan.
   integer(c_int), intent(in)    :: plan_z2d            !< cuFFT backward plan.
   integer(I4P),   intent(in)    :: n1, n2, n3          !< Global cells per direction.
   real(R8P),      intent(in)    :: scale               !< Factor over the unnormalised transform pair.
   real(R8P),      intent(in)    :: lam1(0:n1/2)        !< Eigenvalue part along x.
   real(R8P),      intent(in)    :: lam2(0:n2-1)        !< Eigenvalue part along y.
   real(R8P),      intent(in)    :: lam3(0:n3-1)        !< Eigenvalue part along z.
   real(R8P),      intent(inout) :: f(n1,n2,n3)         !< Source in, potential out.
   complex(R8P),   intent(inout) :: f_hat(n1/2+1,n2,n3) !< Transform workspace.
   real(R8P)                     :: lam                 !< Eigenvalue of the mode.
   integer(c_int)                :: ierr                !< cuFFT status.
   integer(I4P)                  :: i, j, k             !< Counters.

   !$acc data copyin(lam1, lam2, lam3) copy(f) create(f_hat)
   !$acc host_data use_device(f, f_hat)
   ierr = cufftExecD2Z(plan_d2z, f, f_hat)
   !$acc end host_data
   if (ierr /= CUFFT_SUCCESS) call mpih%error_stop(msg='poisson_fft_fnl_object: cufftExecD2Z failed')
   !$acc parallel loop collapse(3) present(f_hat, lam1, lam2, lam3) private(lam)
   do k=1, n3
      do j=1, n2
         do i=1, n1/2+1
            lam = lam1(i-1) + lam2(j-1) + lam3(k-1)
            if (lam == 0._R8P) then
               f_hat(i,j,k) = (0._R8P, 0._R8P)
            else
               f_hat(i,j,k) = f_hat(i,j,k) * (scale / lam)
            endif
         enddo
      enddo
   enddo
   !$acc host_data use_device(f_hat, f)
   ierr = cufftExecZ2D(plan_z2d, f_hat, f)
   !$acc end host_data
   if (ierr /= CUFFT_SUCCESS) call mpih%error_stop(msg='poisson_fft_fnl_object: cufftExecZ2D failed')
   !$acc end data
   endsubroutine spectral_solve_dev

   subroutine spectral_convolve_dev(plan_d2z, plan_z2d, m1, m2, m3, factor, khat, f, f_hat)
   !< The isolated solve on the device: D2Z, the product with the kernel's spectrum (already on the device), Z2D; the
   !< padded source in, the potential out (in place).
   integer(c_int), intent(in)    :: plan_d2z            !< cuFFT forward plan.
   integer(c_int), intent(in)    :: plan_z2d            !< cuFFT backward plan.
   integer(I4P),   intent(in)    :: m1, m2, m3          !< Transform cells per direction.
   real(R8P),      intent(in)    :: factor              !< Source factor.
   real(R8P),      intent(in)    :: khat(m1/2+1,m2,m3)  !< Kernel's spectrum over prod(m), on the device.
   real(R8P),      intent(inout) :: f(m1,m2,m3)         !< Padded source in, potential out.
   complex(R8P),   intent(inout) :: f_hat(m1/2+1,m2,m3) !< Transform workspace.
   integer(c_int)                :: ierr                !< cuFFT status.
   integer(I4P)                  :: i, j, k             !< Counters.

   !$acc data copy(f) create(f_hat)
   !$acc host_data use_device(f, f_hat)
   ierr = cufftExecD2Z(plan_d2z, f, f_hat)
   !$acc end host_data
   if (ierr /= CUFFT_SUCCESS) call mpih%error_stop(msg='poisson_fft_fnl_object: cufftExecD2Z failed')
   !$acc parallel loop collapse(3) present(f_hat, khat)
   do k=1, m3
      do j=1, m2
         do i=1, m1/2+1
            f_hat(i,j,k) = f_hat(i,j,k) * (factor * khat(i,j,k))
         enddo
      enddo
   enddo
   !$acc host_data use_device(f_hat, f)
   ierr = cufftExecZ2D(plan_z2d, f_hat, f)
   !$acc end host_data
   if (ierr /= CUFFT_SUCCESS) call mpih%error_stop(msg='poisson_fft_fnl_object: cufftExecZ2D failed')
   !$acc end data
   endsubroutine spectral_convolve_dev

   subroutine kernel_spectrum_dev(plan_d2z, m1, m2, m3, f, f_hat)
   !< The kernel's forward transform on the device: f in, its spectrum out to the host.
   integer(c_int), intent(in)    :: plan_d2z            !< cuFFT forward plan.
   integer(I4P),   intent(in)    :: m1, m2, m3          !< Transform cells per direction.
   real(R8P),      intent(in)    :: f(m1,m2,m3)         !< Kernel.
   complex(R8P),   intent(inout) :: f_hat(m1/2+1,m2,m3) !< Its spectrum.
   integer(c_int)                :: ierr                !< cuFFT status.

   !$acc data copyin(f) copyout(f_hat)
   !$acc host_data use_device(f, f_hat)
   ierr = cufftExecD2Z(plan_d2z, f, f_hat)
   !$acc end host_data
   !$acc end data
   if (ierr /= CUFFT_SUCCESS) call mpih%error_stop(msg='poisson_fft_fnl_object: cufftExecD2Z failed (kernel)')
   endsubroutine kernel_spectrum_dev
#endif
endmodule adam_fnl_poisson_fft_object
