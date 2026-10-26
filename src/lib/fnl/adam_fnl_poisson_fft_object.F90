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
!< **The distributed slab transform (issue #84, E5d)** keeps the CPU redistribution and return of the block planes
!< (`adam_poisson_fft_object`, on the host); the transform runs on each rank's GPU, its arrays resident there
!< (`slab_r`, `slab_c`, `pen`, `khat_pen`, the transpose buffers: `enter data` with the plans): `slab_r` goes up, the
!< local planes (one D2Z or Z2D plan for one plane, executed on every plane in turn), the transpose to z-lines (packed
!< and unpacked on the device, `MPI_ALLTOALLV` on the device buffers: GPU-aware MPI, as the ghost exchange), the
!< local z-lines (one Z2Z plan for one row of h1 lines, on every row in turn), the spectral solve, back, and `slab_r`
!< comes down. The arithmetic of a value does not depend on the rank count.
!<
!< cuFFT is NVIDIA only: built without `_CUFFT` (an AMD build), making the plans stops the run.
use :: adam_mpih_global, only : mpih
use :: adam_poisson_fft_object, only : poisson_fft_object
#ifdef _CUFFT
use :: cufft
use :: openacc
#endif
use :: iso_c_binding
use :: mpi
use :: penf, only : I4P, I8P, R8P, str

implicit none
private
public :: poisson_fft_fnl_object

type, extends(poisson_fft_object) :: poisson_fft_fnl_object
   !< FNL FFT Poisson solver on uniform periodic grids (cuFFT).
   integer(c_int) :: plan_d2z=0_c_int  !< cuFFT forward plan (rank 0).
   integer(c_int) :: plan_z2d=0_c_int  !< cuFFT backward plan (rank 0).
   logical        :: has_plans=.false. !< The cuFFT plans exist.
   logical        :: khat_on_device=.false. !< `khat` is on the device (isolated).
   integer(c_int) :: plan_sd2z=0_c_int  !< cuFFT slab plan: one plane, D2Z (#84 E5d).
   integer(c_int) :: plan_sz2d=0_c_int  !< cuFFT slab plan: one plane, Z2D.
   integer(c_int) :: plan_sz2z=0_c_int  !< cuFFT slab plan: one row of z-lines, Z2Z.
   logical        :: has_slab_plans=.false. !< The cuFFT slab plans exist.
   logical        :: slab_on_device=.false. !< The slab arrays are on the device.
   integer(I4P), allocatable :: yown(:)     !< Owner rank of each y-row [0:m(2)-1].
   complex(R8P), allocatable :: t_sb(:)     !< Device: the transpose's send buffer.
   complex(R8P), allocatable :: t_rb(:)     !< Device: the transpose's receive buffer.
   contains
      ! backend hooks
      procedure, pass(self) :: destroy_plans => destroy_plans_fnl !< Destroy the cuFFT plans.
      procedure, pass(self) :: make_plans    => make_plans_fnl    !< Make the cuFFT plans (rank 0).
      procedure, pass(self) :: transform     => transform_fnl     !< The spectral solve on the GPU (rank 0).
      procedure, pass(self) :: transform_kernel => transform_kernel_fnl !< The kernel's spectrum on the GPU (isolated).
      procedure, nopass     :: root_cells_max => root_cells_max_fnl   !< 32^3: the device smooths, the host transforms.
      procedure, pass(self) :: slab_make_plans    => slab_make_plans_fnl    !< The cuFFT slab plans (#84 E5d).
      procedure, pass(self) :: slab_destroy_plans => slab_destroy_plans_fnl !< Destroy them.
      procedure, pass(self) :: slab_planes        => slab_planes_fnl        !< Every local plane on the GPU.
      procedure, pass(self) :: slab_lines         => slab_lines_fnl         !< Every local z-line on the GPU.
      procedure, pass(self) :: slab_transpose     => slab_transpose_fnl     !< The transpose on the GPU.
      procedure, pass(self) :: slab_spectral      => slab_spectral_fnl      !< The spectral solve on the GPU.
      procedure, pass(self) :: slab_keep_kernel   => slab_keep_kernel_fnl   !< khat_pen on the GPU.
      procedure, pass(self) :: slab_update        => slab_update_fnl        !< slab_r to and from the GPU.
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

   pure function root_cells_max_fnl(growth) result(n)
   !< The most cells of a multigrid root on cuFFT (#84 E6, E6b): 32^3 whatever the transform's growth. The multigrid
   !< smooths on the device and the root's transform goes through the host's redistribution, so a small root wins
   !< (medooza FNL, 3.9 M AMR cells, 1 / 2 / 4 GPUs: periodic 0.33 / 0.19 / 0.12 s with a 32^3 root, 0.80 / 0.45 /
   !< 0.37 s with 128^3); for the isolated root 16^3 gains only 3-10% over 32^3 (0.83 / 0.51 / 0.39 against 0.85 /
   !< 0.52 / 0.43 s), not worth a uniform 32^3 grid losing its exact one-V-cycle solve.
   integer(I4P), intent(in) :: growth !< Transform cells per root cell (unused: the cap is in cells).
   integer(I4P)             :: n      !< Cells.

   n = 32768_I4P
   endfunction root_cells_max_fnl

   subroutine slab_make_plans_fnl(self)
   !< The cuFFT slab plans (#84 E5d): one plane D2Z and Z2D, one row of z-lines Z2Z (h1 lines of m(3) cells,
   !< contiguous), on the synchronous OpenACC queue's stream; the slab arrays and the transpose buffers on the device.
   class(poisson_fft_fnl_object), intent(inout) :: self !< The solver.
   integer(c_int)                               :: ierr !< cuFFT status.
   integer(I4P)                                 :: r, me, nb !< Counter, rank, buffer size.

   ierr = 0_c_int
#ifdef _CUFFT
   ierr = cufftPlan2d(self%plan_sd2z, int(self%m(2), c_int), int(self%m(1), c_int), CUFFT_D2Z)
   ierr = ierr + cufftPlan2d(self%plan_sz2d, int(self%m(2), c_int), int(self%m(1), c_int), CUFFT_Z2D)
   ierr = ierr + cufftPlanMany(self%plan_sz2z, 1_c_int, [int(self%m(3), c_int)], [int(self%m(3), c_int)], 1_c_int, &
                               int(self%m(3), c_int), [int(self%m(3), c_int)], 1_c_int, int(self%m(3), c_int),    &
                               CUFFT_Z2Z, int(self%h1, c_int))
   if (ierr /= CUFFT_SUCCESS) call mpih%error_stop(msg='poisson_fft_fnl_object%slab_make_plans: cuFFT plans failed')
   ierr = cufftSetStream(self%plan_sd2z, acc_get_cuda_stream(acc_async_sync))
   ierr = ierr + cufftSetStream(self%plan_sz2d, acc_get_cuda_stream(acc_async_sync))
   ierr = ierr + cufftSetStream(self%plan_sz2z, acc_get_cuda_stream(acc_async_sync))
   if (ierr /= CUFFT_SUCCESS) call mpih%error_stop(msg='poisson_fft_fnl_object%slab_make_plans: cufftSetStream failed')
   self%has_slab_plans = .true.
   me = mpih%myrank
   allocate(self%yown(0:self%m(2)-1))
   do r=0, mpih%procs_number-1
      if (self%yr(r+1) > self%yr(r)) self%yown(self%yr(r):self%yr(r+1)-1) = r
   enddo
   nb = max(self%h1 * self%m(2) * (self%zr(me+1) - self%zr(me)), self%h1 * (self%yr(me+1) - self%yr(me)) * self%m(3), 1)
   allocate(self%t_sb(nb), self%t_rb(nb))
   !$acc enter data create(self%slab_r, self%slab_c, self%pen, self%t_sb, self%t_rb)
   self%slab_on_device = .true.
#else
   r = 0 ; me = 0 ; nb = 0
   call mpih%error_stop(msg='poisson_fft_fnl_object: built without cuFFT (build a *-fft mode with -D_CUFFT; cuFFT is '// &
                            'NVIDIA only)')
#endif
   endsubroutine slab_make_plans_fnl

   subroutine slab_destroy_plans_fnl(self)
   !< Destroy the cuFFT slab plans.
   class(poisson_fft_fnl_object), intent(inout) :: self !< The solver.
   integer(c_int)                               :: ierr !< cuFFT status.

#ifdef _CUFFT
   if (self%has_slab_plans) then
      ierr = cufftDestroy(self%plan_sd2z)
      ierr = cufftDestroy(self%plan_sz2d)
      ierr = cufftDestroy(self%plan_sz2z)
   endif
   if (self%slab_on_device) then
      !$acc exit data delete(self%slab_r, self%slab_c, self%pen, self%t_sb, self%t_rb)
      if (allocated(self%khat_pen)) then
         !$acc exit data delete(self%khat_pen)
      endif
   endif
#endif
   ierr = 0_c_int
   self%has_slab_plans = .false.
   self%slab_on_device = .false.
   if (allocated(self%yown)) deallocate(self%yown)
   if (allocated(self%t_sb)) deallocate(self%t_sb)
   if (allocated(self%t_rb)) deallocate(self%t_rb)
   endsubroutine slab_destroy_plans_fnl

   subroutine slab_planes_fnl(self, forward)
   !< Every local plane on the GPU: forward slab_r -> slab_c (D2Z), backward slab_c -> slab_r (Z2D), on the device.
   class(poisson_fft_fnl_object), intent(inout) :: self    !< The solver.
   logical,                       intent(in)    :: forward !< Direction.
   integer(I4P)                                 :: nzl     !< Local planes.

   nzl = self%zr(mpih%myrank+1) - self%zr(mpih%myrank)
   if (nzl == 0) return
#ifdef _CUFFT
   call planes_dev(plan=merge(self%plan_sd2z, self%plan_sz2d, forward), m1=self%m(1), m2=self%m(2), h1=self%h1, &
                   nz=nzl, forward=forward, r=self%slab_r, c=self%slab_c)
#else
   call mpih%error_stop(msg='poisson_fft_fnl_object%slab_planes: built without cuFFT')
#endif
   endsubroutine slab_planes_fnl

   subroutine slab_lines_fnl(self, forward)
   !< Every local z-line on the GPU, in place in pen (Z2Z, forward or inverse), on the device.
   class(poisson_fft_fnl_object), intent(inout) :: self    !< The solver.
   logical,                       intent(in)    :: forward !< Direction.
   integer(I4P)                                 :: nyl     !< Local rows.

   nyl = self%yr(mpih%myrank+1) - self%yr(mpih%myrank)
   if (nyl == 0) return
#ifdef _CUFFT
   call lines_dev(plan=self%plan_sz2z, m3=self%m(3), h1=self%h1, ny=nyl, forward=forward, p=self%pen)
#else
   call mpih%error_stop(msg='poisson_fft_fnl_object%slab_lines: built without cuFFT')
#endif
   endsubroutine slab_lines_fnl

   subroutine slab_transpose_fnl(self, forward)
   !< The transpose on the device: forward slab_c(x, y, z local) -> pen(z, x, y local), backward the reverse; the blocks
   !< packed by destination rank (the CPU layout), MPI_ALLTOALLV on the device buffers (skipped on one rank).
   class(poisson_fft_fnl_object), intent(inout) :: self    !< The solver.
   logical,                       intent(in)    :: forward !< Direction.
   integer(I4P), allocatable :: sc(:), sd(:), rc(:), rd(:) !< Counts and displacements (values).
   integer(I4P)              :: nr, me, r, nzl, nyl, nb    !< Counters, sizes.
   integer(I4P)              :: ierr                       !< MPI error status.

   nr = mpih%procs_number ; me = mpih%myrank
   nzl = self%zr(me+1) - self%zr(me) ; nyl = self%yr(me+1) - self%yr(me)
   nb = size(self%t_sb)
#ifdef _CUFFT
   if (forward) then
      if (nzl > 0) call pack_planes_dev(h1=self%h1, m2=self%m(2), nz=nzl, nr=nr, yr=self%yr, yown=self%yown, &
                                        c=self%slab_c, nb=nb, sb=self%t_sb)
   else
      if (nyl > 0) call pack_lines_dev(m3=self%m(3), h1=self%h1, ny=nyl, p=self%pen, nb=nb, sb=self%t_sb)
   endif
   if (nr > 1) then
      allocate(sc(0:nr-1), sd(0:nr-1), rc(0:nr-1), rd(0:nr-1))
      do r=0, nr-1
         if (forward) then
            sc(r) = self%h1 * (self%yr(r+1) - self%yr(r)) * nzl
            rc(r) = self%h1 * nyl * (self%zr(r+1) - self%zr(r))
         else
            sc(r) = self%h1 * nyl * (self%zr(r+1) - self%zr(r))
            rc(r) = self%h1 * (self%yr(r+1) - self%yr(r)) * nzl
         endif
      enddo
      sd(0) = 0 ; rd(0) = 0
      do r=1, nr-1
         sd(r) = sd(r-1) + sc(r-1) ; rd(r) = rd(r-1) + rc(r-1)
      enddo
      call alltoallv_dev(nb=nb, sb=self%t_sb, sc=sc, sd=sd, rb=self%t_rb, rc=rc, rd=rd, ierr=ierr)
      if (forward) then
         if (nyl > 0) call unpack_lines_dev(m3=self%m(3), h1=self%h1, ny=nyl, nb=nb, rb=self%t_rb, p=self%pen)
      else
         if (nzl > 0) call unpack_planes_dev(h1=self%h1, m2=self%m(2), nz=nzl, nr=nr, yr=self%yr, yown=self%yown, &
                                             nb=nb, rb=self%t_rb, c=self%slab_c)
      endif
   else ! one rank: the packed buffer is the received one
      if (forward) then
         call unpack_lines_dev(m3=self%m(3), h1=self%h1, ny=nyl, nb=nb, rb=self%t_sb, p=self%pen)
      else
         call unpack_planes_dev(h1=self%h1, m2=self%m(2), nz=nzl, nr=nr, yr=self%yr, yown=self%yown, nb=nb, &
                                rb=self%t_sb, c=self%slab_c)
      endif
   endif
#else
   r = 0 ; ierr = 0
   call mpih%error_stop(msg='poisson_fft_fnl_object%slab_transpose: built without cuFFT')
#endif
   endsubroutine slab_transpose_fnl

   subroutine slab_spectral_fnl(self, factor)
   !< The spectral solve on the device: the division by the eigenvalues (the zero modes dropped) or, isolated, the
   !< product with the kernel's spectrum.
   class(poisson_fft_fnl_object), intent(inout) :: self   !< The solver.
   real(R8P),                     intent(in)    :: factor !< Source factor.
   integer(I4P)                                 :: nyl    !< Local rows.

   nyl = self%yr(mpih%myrank+1) - self%yr(mpih%myrank)
   if (nyl == 0) return
#ifdef _CUFFT
   if (self%isolated) then
      call convolve_lines_dev(m3=self%m(3), h1=self%h1, ny=nyl, factor=factor, k=self%khat_pen, p=self%pen)
   else
      call divide_lines_dev(m3=self%m(3), h1=self%h1, ny=nyl, y0=self%yr(mpih%myrank),                       &
                            n1=size(self%lam1), n2=size(self%lam2), n3=size(self%lam3),                       &
                            scale=factor / real(product(int(self%nt, I8P)), R8P),                             &
                            l1=self%lam1, l2=self%lam2, l3=self%lam3, p=self%pen)
   endif
#else
   call mpih%error_stop(msg='poisson_fft_fnl_object%slab_spectral: built without cuFFT, factor '//trim(str(factor)))
#endif
   endsubroutine slab_spectral_fnl

   subroutine slab_keep_kernel_fnl(self)
   !< Isolated, at setup: `khat_pen <- real(pen) / prod(m)` on the device, where it stays.
   class(poisson_fft_fnl_object), intent(inout) :: self !< The solver.

#ifdef _CUFFT
   !$acc enter data create(self%khat_pen)
   call keep_kernel_dev(n=size(self%pen), d=real(product(int(self%m, I8P)), R8P), p=self%pen, k=self%khat_pen)
#else
   call mpih%error_stop(msg='poisson_fft_fnl_object%slab_keep_kernel: built without cuFFT')
#endif
   endsubroutine slab_keep_kernel_fnl

   subroutine slab_update_fnl(self, to_device)
   !< `slab_r` to (`to_device`) or from the device.
   class(poisson_fft_fnl_object), intent(inout) :: self      !< The solver.
   logical,                       intent(in)    :: to_device !< Direction.

   if (self%zr(mpih%myrank+1) == self%zr(mpih%myrank)) return
   if (to_device) then
      !$acc update device(self%slab_r)
   else
      !$acc update self(self%slab_r)
   endif
   endsubroutine slab_update_fnl

#ifdef _CUFFT
   ! private procedures

   subroutine planes_dev(plan, m1, m2, h1, nz, forward, r, c)
   !< The 2-D transform of every plane on the device (r and c present), one plan executed on each plane in turn
   !< (#84 E5d).
   integer(c_int), intent(in)    :: plan            !< cuFFT plan: one plane, D2Z or Z2D.
   integer(I4P),   intent(in)    :: m1, m2, h1, nz  !< Plane sizes, planes.
   logical,        intent(in)    :: forward         !< D2Z (r -> c) or Z2D (c -> r).
   real(R8P),      intent(inout) :: r(m1,m2,nz)     !< Real planes.
   complex(R8P),   intent(inout) :: c(h1,m2,nz)     !< Their transforms.
   integer(c_int)                :: ierr            !< cuFFT status.
   integer(I4P)                  :: k               !< Counter.

   ierr = CUFFT_SUCCESS
   if (forward) then
      !$acc host_data use_device(r, c)
      do k=1, nz
         ierr = ierr + cufftExecD2Z(plan, r(1,1,k), c(1,1,k))
      enddo
      !$acc end host_data
   else
      !$acc host_data use_device(c, r)
      do k=1, nz
         ierr = ierr + cufftExecZ2D(plan, c(1,1,k), r(1,1,k))
      enddo
      !$acc end host_data
   endif
   if (ierr /= CUFFT_SUCCESS) call mpih%error_stop(msg='poisson_fft_fnl_object: slab plane transform failed')
   endsubroutine planes_dev

   subroutine lines_dev(plan, m3, h1, ny, forward, p)
   !< The 1-D transform of every z-line on the device (p present), in place, one plan (a row of h1 lines) on each row
   !< in turn.
   integer(c_int), intent(in)    :: plan         !< cuFFT plan: one row of lines, Z2Z.
   integer(I4P),   intent(in)    :: m3, h1, ny   !< Line length, lines per row, rows.
   logical,        intent(in)    :: forward      !< CUFFT_FORWARD or CUFFT_INVERSE.
   complex(R8P),   intent(inout) :: p(m3,h1,ny)  !< Lines.
   integer(c_int)                :: ierr, dir    !< cuFFT status, direction.
   integer(I4P)                  :: j            !< Counter.

   ierr = CUFFT_SUCCESS
   dir = merge(CUFFT_FORWARD, CUFFT_INVERSE, forward)
   !$acc host_data use_device(p)
   do j=1, ny
      ierr = ierr + cufftExecZ2Z(plan, p(1,1,j), p(1,1,j), dir)
   enddo
   !$acc end host_data
   if (ierr /= CUFFT_SUCCESS) call mpih%error_stop(msg='poisson_fft_fnl_object: slab line transform failed')
   endsubroutine lines_dev

   subroutine pack_planes_dev(h1, m2, nz, nr, yr, yown, c, nb, sb)
   !< Forward transpose, on the device: the planes' rows by destination rank r (the rows [yr(r), yr(r+1)) of every
   !< local plane, x fastest, then y, then z), at r's displacement h1 nz yr(r).
   integer(I4P), intent(in)    :: h1, m2, nz, nr, nb !< Sizes, ranks, buffer size.
   integer(I4P), intent(in)    :: yr(0:nr)           !< First y-row of each rank.
   integer(I4P), intent(in)    :: yown(0:m2-1)       !< Owner rank of each y-row.
   complex(R8P), intent(in)    :: c(h1,m2,nz)        !< Transformed planes (present).
   complex(R8P), intent(inout) :: sb(nb)             !< Send buffer (present).
   integer(I4P)                :: x, y, z, r         !< Counters, rank.

   !$acc parallel loop collapse(3) present(c, sb) copyin(yr, yown) private(r)
   do z=1, nz
      do y=1, m2
         do x=1, h1
            r = yown(y-1)
            sb(h1*nz*yr(r) + x + h1*((y-1-yr(r)) + (yr(r+1)-yr(r))*(z-1))) = c(x,y,z)
         enddo
      enddo
   enddo
   endsubroutine pack_planes_dev

   subroutine unpack_planes_dev(h1, m2, nz, nr, yr, yown, nb, rb, c)
   !< Backward transpose, on the device: the planes' rows from each source rank (the inverse of pack_planes_dev).
   integer(I4P), intent(in)    :: h1, m2, nz, nr, nb !< Sizes, ranks, buffer size.
   integer(I4P), intent(in)    :: yr(0:nr)           !< First y-row of each rank.
   integer(I4P), intent(in)    :: yown(0:m2-1)       !< Owner rank of each y-row.
   complex(R8P), intent(in)    :: rb(nb)             !< Receive buffer (present).
   complex(R8P), intent(inout) :: c(h1,m2,nz)        !< Transformed planes (present).
   integer(I4P)                :: x, y, z, r         !< Counters, rank.

   !$acc parallel loop collapse(3) present(c, rb) copyin(yr, yown) private(r)
   do z=1, nz
      do y=1, m2
         do x=1, h1
            r = yown(y-1)
            c(x,y,z) = rb(h1*nz*yr(r) + x + h1*((y-1-yr(r)) + (yr(r+1)-yr(r))*(z-1)))
         enddo
      enddo
   enddo
   endsubroutine unpack_planes_dev

   subroutine pack_lines_dev(m3, h1, ny, p, nb, sb)
   !< Backward transpose, on the device: the z-lines as (x, y local, z), so the planes of rank r start at h1 ny zr(r).
   integer(I4P), intent(in)    :: m3, h1, ny, nb !< Sizes, buffer size.
   complex(R8P), intent(in)    :: p(m3,h1,ny)    !< Lines (present).
   complex(R8P), intent(inout) :: sb(nb)         !< Send buffer (present).
   integer(I4P)                :: x, y, z        !< Counters.

   !$acc parallel loop collapse(3) present(p, sb)
   do y=1, ny
      do x=1, h1
         do z=1, m3
            sb(x + h1*(y-1) + h1*ny*(z-1)) = p(z,x,y)
         enddo
      enddo
   enddo
   endsubroutine pack_lines_dev

   subroutine unpack_lines_dev(m3, h1, ny, nb, rb, p)
   !< Forward transpose, on the device: the z-lines from (x, y local, z) (the inverse of pack_lines_dev).
   integer(I4P), intent(in)    :: m3, h1, ny, nb !< Sizes, buffer size.
   complex(R8P), intent(in)    :: rb(nb)         !< Receive buffer (present).
   complex(R8P), intent(inout) :: p(m3,h1,ny)    !< Lines (present).
   integer(I4P)                :: x, y, z        !< Counters.

   !$acc parallel loop collapse(3) present(p, rb)
   do y=1, ny
      do x=1, h1
         do z=1, m3
            p(z,x,y) = rb(x + h1*(y-1) + h1*ny*(z-1))
         enddo
      enddo
   enddo
   endsubroutine unpack_lines_dev

   subroutine alltoallv_dev(nb, sb, sc, sd, rb, rc, rd, ierr)
   !< MPI_ALLTOALLV of the device buffers (GPU-aware MPI).
   integer(I4P), intent(in)    :: nb                             !< Buffer size.
   complex(R8P), intent(in)    :: sb(nb)                         !< Send buffer (present).
   integer(I4P), intent(in)    :: sc(0:), sd(0:), rc(0:), rd(0:) !< Counts and displacements (values).
   complex(R8P), intent(inout) :: rb(nb)                         !< Receive buffer (present).
   integer(I4P), intent(out)   :: ierr                           !< MPI error status.

   !$acc host_data use_device(sb, rb)
   call MPI_ALLTOALLV(sb, sc, sd, MPI_DOUBLE_COMPLEX, rb, rc, rd, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
   !$acc end host_data
   endsubroutine alltoallv_dev

   subroutine divide_lines_dev(m3, h1, ny, y0, n1, n2, n3, scale, l1, l2, l3, p)
   !< The periodic spectral solve on the device: p <- p scale / lambda, the zero modes dropped (the CPU slab_spectral).
   integer(I4P), intent(in)    :: m3, h1, ny, y0 !< Sizes, first local row.
   integer(I4P), intent(in)    :: n1, n2, n3     !< Eigenvalue parts' sizes.
   real(R8P),    intent(in)    :: scale          !< Factor over the unnormalised transform pair.
   real(R8P),    intent(in)    :: l1(0:n1-1)     !< Eigenvalue part along x.
   real(R8P),    intent(in)    :: l2(0:n2-1)     !< Eigenvalue part along y.
   real(R8P),    intent(in)    :: l3(0:n3-1)     !< Eigenvalue part along z.
   complex(R8P), intent(inout) :: p(m3,h1,ny)    !< Lines (present).
   real(R8P)                   :: lam            !< Eigenvalue of the mode.
   integer(I4P)                :: i, j, k        !< Counters.

   !$acc parallel loop collapse(3) present(p) copyin(l1, l2, l3) private(lam)
   do j=1, ny
      do i=1, h1
         do k=1, m3
            lam = l1(i-1) + l2(y0+j-1) + l3(k-1)
            if (lam == 0._R8P) then
               p(k,i,j) = (0._R8P, 0._R8P)
            else
               p(k,i,j) = p(k,i,j) * (scale / lam)
            endif
         enddo
      enddo
   enddo
   endsubroutine divide_lines_dev

   subroutine convolve_lines_dev(m3, h1, ny, factor, k, p)
   !< The isolated spectral solve on the device: p <- p factor khat_pen.
   integer(I4P), intent(in)    :: m3, h1, ny  !< Sizes.
   real(R8P),    intent(in)    :: factor      !< Source factor.
   real(R8P),    intent(in)    :: k(m3,h1,ny) !< Kernel's spectrum over prod(m) (present).
   complex(R8P), intent(inout) :: p(m3,h1,ny) !< Lines (present).
   integer(I4P)                :: i, j, z     !< Counters.

   !$acc parallel loop collapse(3) present(p, k)
   do j=1, ny
      do i=1, h1
         do z=1, m3
            p(z,i,j) = p(z,i,j) * (factor * k(z,i,j))
         enddo
      enddo
   enddo
   endsubroutine convolve_lines_dev

   subroutine keep_kernel_dev(n, d, p, k)
   !< k <- real(p) / d on the device (the kernel's spectrum, real: the kernel is even).
   integer(I4P), intent(in)    :: n    !< Values.
   real(R8P),    intent(in)    :: d    !< Divisor, prod(m).
   complex(R8P), intent(in)    :: p(n) !< Transformed kernel (present).
   real(R8P),    intent(inout) :: k(n) !< Its spectrum (present).
   integer(I4P)                :: i    !< Counter.

   !$acc parallel loop present(p, k)
   do i=1, n
      k(i) = real(p(i), R8P) / d
   enddo
   endsubroutine keep_kernel_dev
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
