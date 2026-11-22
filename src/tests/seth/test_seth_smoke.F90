!< Smoke test of the SETH dependency: the WENO-JS kernels on the host and on the device give the same faces (ADAM #139, S1).
program test_seth_smoke
!< Smoke test of the SETH dependency: the WENO-JS kernels on the host and on the device give the same faces (ADAM #139, S1).
!<
!< **Why this test exists.** ADAM compiles SETH from source with its own flags: the host build in the CPU modes, the
!< OpenACC build in the FNL modes (`-DDEV_OAC`), the OpenMP-offload build in the AMD modes (`-DDEV_OMP`). This program is
!< the first ADAM target that reaches SETH: it proves that the dependency is fetched, compiled and linked on each backend,
!< and that its kernels run on the device, before any ADAM object uses them.
!<
!< **What it does.** N random faces (smooth and rough stencils, magnitudes 1e-3..1e3) for every WENO-JS width (S = 1..5),
!< each reconstructed by `seth_weno<2S-1>_js_wr` with the scale-invariant descaler `seth_weno_descaling` (sigma = 1), on
!< the host and in one device loop (OpenACC or OpenMP target; a host loop in the CPU build). Check: the two sets of states
!< and weights ratios agree bitwise in the CPU build (two host loops of one compiler); in the device builds the device
!< compiler contracts other products into FMAs than the host one (nvfortran -fast), so the check there is a difference
!< of at most RTOL = 8 epsilon, the states relative to the largest |value| of their stencil, the weights ratios (in
!< [0, 1]) absolute.
use penf, only : I4P, R8P
use seth, only : R_P, SETH_WENO_JS_EPS, seth_weno_descaling,                                      &
                 seth_weno1_js_wr, seth_weno3_js_wr, seth_weno5_js_wr, seth_weno7_js_wr, seth_weno9_js_wr

implicit none

integer(I4P), parameter :: N=4096           !< Faces per width.
integer(I4P), parameter :: SMAX=5           !< Widest stencil (WENO-9).
real(R8P),    allocatable :: vp(:,:), vm(:,:) !< Stencils, (1-SMAX:SMAX, N).
real(R8P),    allocatable :: hst(:,:), dev(:,:) !< Results (rp, rm, wrp, wrm, N), host and device.
#if defined DEV_OAC || defined DEV_OMP
real(R8P),    parameter   :: RTOL=8._R8P*epsilon(1._R8P) !< Device builds: largest relative difference allowed.
#else
real(R8P),    parameter   :: RTOL=0._R8P                 !< CPU build: bitwise.
#endif
real(R8P)                 :: dmax, rmax     !< Largest absolute and relative differences.
integer(I4P)              :: S, f, nbad, nfail !< Width, face, counters.

if (R_P /= R8P) error stop 'test_seth_smoke: SETH must be built with R8P in ADAM'
#if defined DEV_OAC
print '(A)', 'FNL build (OpenACC): the device loop runs on the GPU'
#elif defined DEV_OMP
print '(A)', 'AMD build (OpenMP offload): the device loop runs on the GPU'
#else
print '(A)', 'CPU build: the "device" loop is a second host loop'
#endif
allocate(vp(1-SMAX:SMAX,N), vm(1-SMAX:SMAX,N), hst(4,N), dev(4,N))
call fill_stencils(vp=vp, vm=vm)
nfail = 0
do S=1, SMAX
   do f=1, N
      call face(S=S, vp=vp(:,f), vm=vm(:,f), r=hst(:,f))
   enddo
   dev = 0._R8P
#if defined DEV_OAC
   !$acc parallel loop copyin(vp, vm) copyout(dev)
#elif defined DEV_OMP
   !$omp target teams distribute parallel do map(to:vp, vm) map(from:dev)
#endif
   do f=1, N
      call face(S=S, vp=vp(:,f), vm=vm(:,f), r=dev(:,f))
   enddo
   nbad = count(hst /= dev)
   dmax = maxval(abs(hst - dev))
   rmax = 0._R8P
   do f=1, N ! the states relative to the magnitude of the face's stencil, the weights ratios (in [0, 1]) absolute
      rmax = max(rmax, maxval(abs(hst(1:2,f) - dev(1:2,f)))/maxval(abs(vp(:,f))), maxval(abs(hst(3:4,f) - dev(3:4,f))))
   enddo
   print '(A,I1,A,I0,A,I0,A,ES10.3,A,ES10.3,A,F5.2,A)', 'WENO-', 2*S-1, ': ', nbad, ' of ', 4*N, &
      ' values differ, largest ', dmax, ', relative ', rmax, ' (', rmax/epsilon(1._R8P), ' eps)'
   if (rmax > RTOL) nfail = nfail + 1
enddo
if (nfail > 0) error stop 'test_seth_smoke: host and device differ beyond the tolerance'
print '(A)', 'test_seth_smoke: ok'

contains
   pure subroutine face(S, vp, vm, r)
   !< One face of width S: the descaler, then the WENO-JS kernel of that width with its weights ratios.
   integer(I4P), intent(in)  :: S                 !< Candidates of the scheme.
   real(R8P),    intent(in)  :: vp(1-SMAX:SMAX)   !< Values with a left bias.
   real(R8P),    intent(in)  :: vm(1-SMAX:SMAX)   !< Values with a right bias.
   real(R8P),    intent(out) :: r(4)              !< rp, rm, wrp, wrm.
   real(R8P)                 :: rmu               !< Inverse descaler.
   !$acc routine seq
   !$omp declare target

   rmu = seth_weno_descaling(S=S, sigma=1._R8P, vp=vp(1-S:S), vm=vm(1-S:S))
   select case(S)
   case(1) ; call seth_weno1_js_wr(vp(1-S:S), vm(1-S:S), SETH_WENO_JS_EPS, rmu, r(1), r(2), r(3), r(4))
   case(2) ; call seth_weno3_js_wr(vp(1-S:S), vm(1-S:S), SETH_WENO_JS_EPS, rmu, r(1), r(2), r(3), r(4))
   case(3) ; call seth_weno5_js_wr(vp(1-S:S), vm(1-S:S), SETH_WENO_JS_EPS, rmu, r(1), r(2), r(3), r(4))
   case(4) ; call seth_weno7_js_wr(vp(1-S:S), vm(1-S:S), SETH_WENO_JS_EPS, rmu, r(1), r(2), r(3), r(4))
   case(5) ; call seth_weno9_js_wr(vp(1-S:S), vm(1-S:S), SETH_WENO_JS_EPS, rmu, r(1), r(2), r(3), r(4))
   endselect
   endsubroutine face

   subroutine fill_stencils(vp, vm)
   !< Random stencils, a fixed seed: half smooth (a sampled sine plus 1e-3 noise), half rough (a jump in the middle), each
   !< face scaled by 10**k, k uniform in [-3, 3].
   real(R8P), intent(out) :: vp(1-SMAX:,:), vm(1-SMAX:,:) !< Stencils.
   real(R8P)              :: u(2*SMAX), a, ph, sc         !< Random numbers, amplitude, phase, scale.
   integer, allocatable   :: seed(:)                      !< Seed.
   integer                :: ns, f, m                     !< Seed size, counters.

   call random_seed(size=ns) ; allocate(seed(ns)) ; seed = 139 ; call random_seed(put=seed)
   do f=1, size(vp, dim=2)
      call random_number(u) ; call random_number(a) ; call random_number(ph) ; call random_number(sc)
      sc = 10._R8P**(6._R8P*sc - 3._R8P)
      do m=1-SMAX, SMAX
         if (mod(f, 2) == 0) then
            vp(m,f) = sc*(1._R8P + a*sin(0.7_R8P*m + 6.28_R8P*ph) + 1.e-3_R8P*u(m+SMAX))
         else
            vp(m,f) = sc*merge(1._R8P, 0.1_R8P, m <= 0) + sc*1.e-2_R8P*u(m+SMAX)
         endif
      enddo
      vm(:,f) = vp(:,f)*(1._R8P + 1.e-3_R8P*(u(1) - 0.5_R8P))
   enddo
   endsubroutine fill_stencils
endprogram test_seth_smoke
