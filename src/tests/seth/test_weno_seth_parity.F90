!< Device-callable faces of the SETH parity test (ADAM #139, S2): ADAM's FNL per-width kernels and SETH's, one face each.
module test_weno_seth_parity_kernels
!< Device-callable faces of the SETH parity test (ADAM #139, S2): ADAM's FNL per-width kernels and SETH's, one face each.
!<
!< `seth_face` reconstructs one face through the adapter's run-time dispatchers (`weno_reconstruct`, `_wr`, `_interp_wr`
!< of `adam_weno_object`) with SETH's descaler; `adam_dev_face` (device builds only) through ADAM's FNL kernels of the
!< same width (`weno_reconstruct_upwind[_wratio]_s<N>_dev`) on the tables of a `weno_object`. `kind`: 1 reconstruction,
!< 2 reconstruction with the weights ratios, 3 interpolation with the weights ratios. `r = [rp, rm, wrp, wrm]`.
use :: adam_weno_object, only : weno_reconstruct, weno_reconstruct_wr, weno_reconstruct_interp_wr
#if defined DEV_OAC || defined DEV_OMP
use :: adam_fnl_weno_width_kernels, only : weno_reconstruct_upwind_s1_dev, weno_reconstruct_upwind_s2_dev,                &
                                           weno_reconstruct_upwind_s3_dev, weno_reconstruct_upwind_s4_dev,                &
                                           weno_reconstruct_upwind_s5_dev,                                                &
                                           weno_reconstruct_upwind_wratio_s1_dev, weno_reconstruct_upwind_wratio_s2_dev,  &
                                           weno_reconstruct_upwind_wratio_s3_dev, weno_reconstruct_upwind_wratio_s4_dev,  &
                                           weno_reconstruct_upwind_wratio_s5_dev
#endif
use :: penf, only : I4P, R8P
use :: seth, only : seth_weno_descaling

implicit none
private
public :: SMAX
public :: seth_face
#if defined DEV_OAC || defined DEV_OMP
public :: adam_dev_face
#endif

integer(I4P), parameter :: SMAX=5 !< Widest stencil (WENO-9).

contains
   pure subroutine seth_face(S, kind, sigma, eps, vp, vm, r)
   !< One face by the adapter's dispatchers (SETH's kernels) and SETH's descaler.
   integer(I4P), intent(in)  :: S               !< Candidates of the scheme.
   integer(I4P), intent(in)  :: kind            !< 1 reconstruction, 2 with the weights ratios, 3 interpolation.
   real(R8P),    intent(in)  :: sigma           !< Descaler switch.
   real(R8P),    intent(in)  :: eps             !< Regulariser.
   real(R8P),    intent(in)  :: vp(1-SMAX:SMAX) !< Left-biased values.
   real(R8P),    intent(in)  :: vm(1-SMAX:SMAX) !< Right-biased values.
   real(R8P),    intent(out) :: r(4)            !< rp, rm, wrp, wrm.
   real(R8P)                 :: rmu             !< Inverse descaler.
   !$acc routine seq
   !$omp declare target

   rmu = seth_weno_descaling(S=S, sigma=sigma, vp=vp(1-S:S), vm=vm(1-S:S))
   r = 0._R8P
   select case(kind)
   case(1) ; call weno_reconstruct(S=S, vp=vp(1-S:S), vm=vm(1-S:S), eps=eps, rmu=rmu, rp=r(1), rm=r(2))
   case(2) ; call weno_reconstruct_wr(S=S, vp=vp(1-S:S), vm=vm(1-S:S), eps=eps, rmu=rmu, rp=r(1), rm=r(2), wrp=r(3), wrm=r(4))
   case(3) ; call weno_reconstruct_interp_wr(S=S, vp=vp(1-S:S), vm=vm(1-S:S), eps=eps, rmu=rmu, &
                                             rp=r(1), rm=r(2), wrp=r(3), wrm=r(4))
   endselect
   endsubroutine seth_face

#if defined DEV_OAC || defined DEV_OMP
   subroutine adam_dev_face(S, kind, rmu, zeps, a, p, d, v, r)
   !< One face by ADAM's FNL kernel of width S, ADAM's rows `v(1:2,1-S:S-1)` (row 2 = vp, row 1 = vm shifted one cell).
   integer(I4P), intent(in)  :: S                   !< Candidates of the scheme.
   integer(I4P), intent(in)  :: kind                !< 1 reconstruction, 2 and 3 with the weights ratios.
   real(R8P),    intent(in)  :: rmu                 !< Inverse descaler.
   real(R8P),    intent(in)  :: zeps                !< Regulariser.
   real(R8P),    intent(in)  :: a(1:,0:,1:)         !< Optimal weights (reconstruction or interpolation).
   real(R8P),    intent(in)  :: p(1:,0:,0:,1:)      !< Polynomial coefficients (reconstruction or interpolation).
   real(R8P),    intent(in)  :: d(0:,0:,0:,1:)      !< Smoothness indicator coefficients.
   real(R8P),    intent(in)  :: v(1:2,1-SMAX:SMAX-1) !< ADAM's rows.
   real(R8P),    intent(out) :: r(4)                !< rp, rm, wrp, wrm.
   real(R8P)                 :: vr(2), wr(2)        !< ADAM's states and ratios.
   !$acc routine seq
   !$omp declare target

   vr = 0._R8P ; wr = 0._R8P
   if (kind == 1) then
      select case(S)
      case(1) ; call weno_reconstruct_upwind_s1_dev(a, p, d, zeps, rmu, v(:,1-1:1-1), vr)
      case(2) ; call weno_reconstruct_upwind_s2_dev(a, p, d, zeps, rmu, v(:,1-2:2-1), vr)
      case(3) ; call weno_reconstruct_upwind_s3_dev(a, p, d, zeps, rmu, v(:,1-3:3-1), vr)
      case(4) ; call weno_reconstruct_upwind_s4_dev(a, p, d, zeps, rmu, v(:,1-4:4-1), vr)
      case(5) ; call weno_reconstruct_upwind_s5_dev(a, p, d, zeps, rmu, v(:,1-5:5-1), vr)
      endselect
   else
      select case(S)
      case(1) ; call weno_reconstruct_upwind_wratio_s1_dev(a, p, d, zeps, rmu, v(:,1-1:1-1), vr, wr)
      case(2) ; call weno_reconstruct_upwind_wratio_s2_dev(a, p, d, zeps, rmu, v(:,1-2:2-1), vr, wr)
      case(3) ; call weno_reconstruct_upwind_wratio_s3_dev(a, p, d, zeps, rmu, v(:,1-3:3-1), vr, wr)
      case(4) ; call weno_reconstruct_upwind_wratio_s4_dev(a, p, d, zeps, rmu, v(:,1-4:4-1), vr, wr)
      case(5) ; call weno_reconstruct_upwind_wratio_s5_dev(a, p, d, zeps, rmu, v(:,1-5:5-1), vr, wr)
      endselect
   endif
   r = [vr(2), vr(1), wr(2), wr(1)]
   endsubroutine adam_dev_face
#endif
endmodule test_weno_seth_parity_kernels

program test_weno_seth_parity
!< Parity of the SETH adapter with ADAM's WENO, hybrid and seam sensor (ADAM #139, S2; a transition test, deleted in S6).
!<
!< **Why this test exists.** #139 replaces ADAM's WENO kernels, its KEEP-PE flux, its Ducros sensor and its seam WENO-5
!< sensor by SETH's, under the rule that the results do not change: CPU runs bitwise, FNL runs bitwise expected (D7).
!< Before any application calls SETH, this program measures the substitution on random data, so that a later A/B
!< difference has one place to be explained.
!<
!< **What it checks.** Random stencils for S = 1..5 (half smooth, half rough, each scaled by 10^k, k uniform in
!< [-12, 12]; `vp` and `vm` independent, as a flux split's two halves are):
!< 1. host, gated (bitwise in the CPU build, 8 eps in the device builds, whose host compiler contracts FMAs): the
!<    adapter's dispatchers (`weno_reconstruct`, `_wr`, `_interp_wr`) against ADAM's host routines
!<    (`weno_reconstruct_upwind[_wratio]` on the reconstruction and on the interpolation tables), for the Jiang-Shu and
!<    the scale-invariant weights (`seth_weno_descaling` against `weno_descaling`, under the same gate: nvfortran -fast
!<    reassociates the host sum);
!< 2. host, gated bitwise in the CPU build: `seth_keep_pe6` against `compute_kep_flux` (20000 states, each direction),
!<    `seth_ducros` plus FLUME's density-jump override against `compute_ducros`;
!< 3. host: `min(wrp, wrm)` of `seth_weno5_js_wr` against FLUME's `seam_weno5_sensor`, the largest difference reported
!<    (the two sum the indicators' products in another order, and on smooth data the sums cancel: 3.3e-12 measured), the
!<    faces whose `reflux_policy = troubled` decision (sensor >= 0.2) flips gated at none;
!< 4. device builds: the dispatchers on the device against the host (8 eps), and SETH's kernels on the device against
!<    ADAM's FNL per-width kernels on the device (reported: the predictor of the FNL A/B of #139 S3b, D7).
use :: adam_globals,                  only : mpih
use :: adam_flume_hybrid_library,     only : compute_ducros, compute_kep_flux
use :: adam_flume_seam_sync_object,   only : seam_weno5_sensor
use :: adam_weno_object,              only : weno_descaling, weno_object, weno_reconstruct_upwind,                        &
                                             weno_reconstruct_upwind_wratio
use :: penf,                          only : I4P, R8P
use :: seth,                          only : SETH_WENO_JS_EPS, seth_ducros, seth_keep_pe6, seth_weno_descaling,           &
                                             seth_weno5_js_wr
use :: test_weno_seth_parity_kernels, only : SMAX, seth_face
#if defined DEV_OAC || defined DEV_OMP
use :: test_weno_seth_parity_kernels, only : adam_dev_face
#endif

implicit none

integer(I4P),     parameter :: N=20000                       !< Faces per width and weights.
character(len=8), parameter :: SCHEMES(5)=['weno-u-1', 'weno-u-3', 'weno-u-5', 'weno-u-7', 'weno-u-9'] !< One per S.
character(len=6), parameter :: KINDS(3)=['recon ', 'recwr ', 'intwr ']  !< Kinds of face.
#if defined DEV_OAC || defined DEV_OMP
real(R8P),        parameter :: RTOL_HOST=8._R8P*epsilon(1._R8P) !< Host gate of the device builds.
#else
real(R8P),        parameter :: RTOL_HOST=0._R8P                 !< Host gate of the CPU build: bitwise.
#endif
real(R8P),        parameter :: RTOL_DEV=8._R8P*epsilon(1._R8P)  !< Device against host.
type(weno_object)           :: weno(5)                         !< ADAM's WENO objects, one per S (the tables).
real(R8P),        allocatable :: vp(:,:), vm(:,:)              !< Stencils (1-SMAX:SMAX, N).
real(R8P),        allocatable :: v(:,:,:)                      !< ADAM's rows (1:2, 1-SMAX:SMAX-1, N).
real(R8P),        allocatable :: ra(:,:), rs(:,:), rd(:,:)     !< Results (4, N): ADAM host, SETH host, device.
real(R8P)                   :: sigma, zeps, rmu_a, rmu_s, dmax, rel !< Descaler switch, regulariser, descalers, diffs.
integer(I4P)                :: S, k, f, iw, nbad, nfail, nrmu  !< Width, kind, face, weights, counters.
logical                     :: dev_built                       !< Device build flag.

call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)
dev_built = .false.
#if defined DEV_OAC || defined DEV_OMP
dev_built = .true.
#endif
nfail = 0
allocate(vp(1-SMAX:SMAX,N), vm(1-SMAX:SMAX,N), v(1:2,1-SMAX:SMAX-1,N), ra(4,N), rs(4,N), rd(4,N))
do S=1, 5
   call weno(S)%initialize(scheme=trim(SCHEMES(S)), nb=1, ngc=5, ni=1, nj=1, nk=1)
   call weno(S)%initialize_interpolation
enddo
zeps = weno(1)%zeps
if (zeps /= SETH_WENO_JS_EPS .or. weno(1)%eps /= SETH_WENO_JS_EPS) then
   print '(A)', 'FAIL: weno_object zeps/eps differ from SETH_WENO_JS_EPS' ; nfail = nfail + 1
endif
call fill_stencils(vp=vp, vm=vm)
do f=1, N ! ADAM's rows: row 2 = vp, row 1 = vm shifted one cell
   v(2,1-SMAX:SMAX-1,f) = vp(1-SMAX:SMAX-1,f)
   v(1,1-SMAX:SMAX-1,f) = vm(2-SMAX:SMAX,f)
enddo

! 1. the WENO dispatchers against ADAM's host routines (and the device, in the device builds)
do iw=0, 1
   sigma = real(iw, R8P)
   do S=1, 5
      nrmu = 0
      do k=1, 3
         do f=1, N
            rmu_a = weno_descaling(S=S, weno_sigma=sigma, v=v(:,1-S:S-1,f))
            rmu_s = seth_weno_descaling(S=S, sigma=sigma, vp=vp(1-S:S,f), vm=vm(1-S:S,f))
            if (abs(rmu_a - rmu_s) > RTOL_HOST*abs(rmu_a)) nrmu = nrmu + 1 ! nvfortran -fast reassociates the host sum
            call adam_host_face(S=S, kind=k, rmu=rmu_a, f=f, r=ra(:,f))
            call seth_face(S=S, kind=k, sigma=sigma, eps=SETH_WENO_JS_EPS, vp=vp(:,f), vm=vm(:,f), r=rs(:,f))
         enddo
         call compare(label='host '//trim(KINDS(k)), a=ra, b=rs, rtol=RTOL_HOST, gate=.true.)
         if (dev_built) call device_checks(S=S, k=k)
      enddo
      if (nrmu > 0) then
         print '(A,I0,A,I0)', 'FAIL: descaler differs on ', nrmu, ' faces, S = ', S ; nfail = nfail + 1
      endif
   enddo
enddo

! 2. the hybrid's flux and sensor
call hybrid_checks
! 3. the seam sensor
call seam_sensor_check

call mpih%finalize
if (nfail > 0) then
   print '(A,I0,A)', 'test_weno_seth_parity: FAILED (', nfail, ' checks)'
   error stop 1
endif
print '(A)', 'test_weno_seth_parity: ok'

contains
   subroutine adam_host_face(S, kind, rmu, f, r)
   !< One face by ADAM's host routines on the tables of weno(S).
   integer(I4P), intent(in)  :: S, kind, f !< Width, kind, face.
   real(R8P),    intent(in)  :: rmu        !< Inverse descaler.
   real(R8P),    intent(out) :: r(4)       !< rp, rm, wrp, wrm.
   real(R8P)                 :: vr(2), wr(2) !< States and ratios, ADAM's rows.

   vr = 0._R8P ; wr = 0._R8P
   associate(w => weno(S))
      select case(kind)
      case(1) ; call weno_reconstruct_upwind(S=S, weno_a=w%a, weno_p=w%p, weno_d=w%d, weno_zeps=w%zeps, weno_rmu=rmu, &
                                             v=v(:,1-S:S-1,f), vr=vr)
      case(2) ; call weno_reconstruct_upwind_wratio(S=S, weno_a=w%a, weno_p=w%p, weno_d=w%d, weno_zeps=w%zeps,        &
                                                    weno_rmu=rmu, v=v(:,1-S:S-1,f), vr=vr, wr=wr)
      case(3) ; call weno_reconstruct_upwind_wratio(S=S, weno_a=w%a_interp, weno_p=w%p_interp, weno_d=w%d,            &
                                                    weno_zeps=w%zeps, weno_rmu=rmu, v=v(:,1-S:S-1,f), vr=vr, wr=wr)
      endselect
   endassociate
   r = [vr(2), vr(1), wr(2), wr(1)]
   endsubroutine adam_host_face

   subroutine compare(label, a, b, rtol, gate)
   !< Count the faces that differ, the largest difference (states relative to the stencil magnitude, ratios absolute).
   character(*), intent(in) :: label    !< Report label.
   real(R8P),    intent(in) :: a(:,:)   !< Reference (4, N).
   real(R8P),    intent(in) :: b(:,:)   !< Candidate (4, N).
   real(R8P),    intent(in) :: rtol     !< Largest difference passed.
   logical,      intent(in) :: gate     !< Count a miss as a failure.
   integer(I4P)             :: ff       !< Face.

   nbad = count(a /= b) ; dmax = 0._R8P
   do ff=1, N
      rel = max(maxval(abs(a(1:2,ff) - b(1:2,ff)))/max(maxval(abs(vp(:,ff))), maxval(abs(vm(:,ff)))), &
                maxval(abs(a(3:4,ff) - b(3:4,ff))))
      dmax = max(dmax, rel)
   enddo
   print '(A,A,I0,A,F3.1,A,A,A,I0,A,I0,A,ES10.3,A)', label, ' S=', S, ' sigma=', sigma, ' ', &
      merge('gated   ', 'reported', gate), ': ', nbad, ' of ', 4*N, ' values differ, largest ', dmax, ' (scaled)'
   if (gate .and. dmax > rtol) then
      print '(A)', '   FAIL' ; nfail = nfail + 1
   endif
   endsubroutine compare

   subroutine device_checks(S, k)
   !< The dispatchers on the device against the host; SETH on the device against ADAM's FNL kernels on the device.
   integer(I4P), intent(in) :: S, k !< Width, kind.
#if defined DEV_OAC || defined DEV_OMP
   real(R8P), allocatable   :: a(:,:,:), p(:,:,:,:), d(:,:,:,:) !< ADAM's tables of this kind.
   real(R8P)                :: s_sigma, s_eps, rmu              !< Copies for the device.
   integer(I4P)             :: ff, s_S, s_k                     !< Face, copies.

   s_sigma = sigma ; s_eps = SETH_WENO_JS_EPS ; s_S = S ; s_k = k
   rd = 0._R8P
#if defined DEV_OAC
   !$acc parallel loop copyin(vp, vm) copyout(rd)
#else
   !$omp target teams distribute parallel do map(to:vp, vm) map(from:rd)
#endif
   do ff=1, N
      call seth_face(S=s_S, kind=s_k, sigma=s_sigma, eps=s_eps, vp=vp(:,ff), vm=vm(:,ff), r=rd(:,ff))
   enddo
   call compare(label='dev SETH vs host SETH '//trim(KINDS(k)), a=rs, b=rd, rtol=RTOL_DEV, gate=.true.)
   ra = rd ! SETH on the device, the reference of the next comparison
   if (k == 3) then
      a = weno(S)%a_interp ; p = weno(S)%p_interp
   else
      a = weno(S)%a ; p = weno(S)%p
   endif
   d = weno(S)%d
   rd = 0._R8P
#if defined DEV_OAC
   !$acc parallel loop copyin(a, p, d, v, vp, vm) copyout(rd) private(rmu)
#else
   !$omp target teams distribute parallel do map(to:a, p, d, v, vp, vm) map(from:rd) private(rmu)
#endif
   do ff=1, N
      rmu = seth_weno_descaling(S=s_S, sigma=s_sigma, vp=vp(1-s_S:s_S,ff), vm=vm(1-s_S:s_S,ff))
      call adam_dev_face(S=s_S, kind=s_k, rmu=rmu, zeps=s_eps, a=a, p=p, d=d, v=v(:,:,ff), r=rd(:,ff))
   enddo
   call compare(label='dev ADAM-FNL vs dev SETH '//trim(KINDS(k)), a=ra, b=rd, rtol=0._R8P, gate=.false.)
#endif
   endsubroutine device_checks

   subroutine hybrid_checks
   !< seth_keep_pe6 against compute_kep_flux; seth_ducros + the density-jump override against compute_ducros.
   real(R8P)    :: w(5,-2:3), fa(5), fs(5), g(3,3), u(3,3), pa, ps, jump, eps, dk, dd !< States, fluxes, sensors.
   integer(I4P) :: i, dir, nk, nd                                                  !< Counters.

   nk = 0 ; nd = 0 ; dk = 0._R8P ; dd = 0._R8P
   do i=1, N
      call random_number(w)
      w(1,:) = 0.1_R8P + 10._R8P*w(1,:) ; w(2:4,:) = 2._R8P*w(2:4,:) - 1._R8P ; w(5,:) = 0.1_R8P + 10._R8P*w(5,:)
      do dir=1, 3
         call compute_kep_flux(d=dir, gamma=1.4_R8P, w=w, f=fa)
         call seth_keep_pe6(d=dir, gamma=1.4_R8P, w=w, f=fs)
         if (any(fa /= fs)) nk = nk + 1
         dk = max(dk, maxval(abs(fa - fs)))
      enddo
      call random_number(u) ; call random_number(g) ; call random_number(jump)
      g = (2._R8P*g - 1._R8P)*10._R8P**(6._R8P*u - 3._R8P) ; jump = 0.2_R8P*jump ; eps = 1.e-8_R8P
      call compute_ducros(eps=eps, jump=jump, jump_max=0.1_R8P, g=g, psi=pa)
      call seth_ducros(eps=eps, g=g, phi=ps)
      if (jump > 0.1_R8P) ps = 1._R8P
      if (pa /= ps) nd = nd + 1
      dd = max(dd, abs(pa - ps))
   enddo
   print '(A,I0,A,I0,A,ES10.3)', 'host keep-pe6 vs compute_kep_flux gated  : ', nk, ' of ', 3*N, ' faces differ, largest ', dk
   print '(A,I0,A,I0,A,ES10.3)', 'host ducros vs compute_ducros gated      : ', nd, ' of ', N, ' cells differ, largest ', dd
   if (.not.dev_built .and. (nk > 0 .or. nd > 0)) then
      print '(A)', '   FAIL' ; nfail = nfail + 1
   endif
   endsubroutine hybrid_checks

   subroutine seam_sensor_check
   !< min(wrp, wrm) of seth_weno5_js_wr against seam_weno5_sensor; the troubled decision (>= 0.2) that flips.
   real(R8P)    :: c(-2:3), u(-2:3), sa, ss, rp, rm, wrp, wrm, sc, dsm !< Cells, sensors, scale, largest difference.
   integer(I4P) :: i, m, nb_, nflip                                   !< Counters.

   nb_ = 0 ; nflip = 0 ; dsm = 0._R8P
   do i=1, N
      call random_number(u) ; call random_number(sc)
      sc = 10._R8P**(6._R8P*sc - 3._R8P)
      do m=-2, 3
         if (mod(i, 2) == 0) then
            c(m) = sc*(1._R8P + 0.3_R8P*sin(0.4_R8P*m + 6.28_R8P*u(-2)) + 1.e-3_R8P*u(m))
         else
            c(m) = sc*merge(1._R8P, 0.125_R8P, m <= mod(i/2, 4) - 2) + sc*1.e-2_R8P*u(m)
         endif
      enddo
      sa = seam_weno5_sensor(v=c)
      call seth_weno5_js_wr(vp=c, vm=c, eps=SETH_WENO_JS_EPS, rmu=1._R8P, rp=rp, rm=rm, wrp=wrp, wrm=wrm)
      ss = min(wrp, wrm)
      if (sa /= ss) nb_ = nb_ + 1
      dsm = max(dsm, abs(sa - ss))
      if ((sa >= 0.2_R8P) .neqv. (ss >= 0.2_R8P)) nflip = nflip + 1
   enddo
   print '(A,I0,A,I0,A,ES10.3,A,ES9.2,A,I0)', 'host seam sensor vs seth_weno5_js_wr      : ', nb_, ' of ', N, &
      ' differ, largest ', dsm, ' (', dsm/epsilon(1._R8P), ' eps), troubled decisions flipped ', nflip
   ! gated on the decision only: the hand-expanded indicators and SETH's table sum the same products in another order,
   ! and on smooth data the sums cancel (IS << v^2), so the sensors differ far above one ulp (measured 3.3e-12 on this
   ! sample); what reflux_policy = troubled reads is the comparison with its threshold
   if (nflip > 0) then
      print '(A)', '   FAIL' ; nfail = nfail + 1
   endif
   endsubroutine seam_sensor_check

   subroutine fill_stencils(vp, vm)
   !< Random stencils, a fixed seed: half smooth (a sampled sine plus 1e-3 noise), half rough (a jump at a random cell),
   !< each face scaled by 10^k, k uniform in [-12, 12]; vp and vm independent.
   real(R8P), intent(out) :: vp(1-SMAX:,:), vm(1-SMAX:,:) !< Stencils.
   real(R8P)              :: u(2*SMAX), a, ph, sc, j       !< Random numbers.
   integer, allocatable   :: seed(:)                      !< Seed.
   integer                :: ns, ff, m, h                 !< Seed size, counters.

   call random_seed(size=ns) ; allocate(seed(ns)) ; seed = 139 ; call random_seed(put=seed)
   do ff=1, size(vp, dim=2)
      do h=1, 2
         call random_number(u) ; call random_number(a) ; call random_number(ph) ; call random_number(sc)
         call random_number(j)
         sc = 10._R8P**(24._R8P*sc - 12._R8P)
         do m=1-SMAX, SMAX
            if (mod(ff, 2) == 0) then
               u(m+SMAX) = sc*(1._R8P + a*sin(0.7_R8P*m + 6.28_R8P*ph) + 1.e-3_R8P*u(m+SMAX))
            else
               u(m+SMAX) = sc*merge(1._R8P, 0.1_R8P, m <= int(j*2*SMAX) - SMAX) + sc*1.e-2_R8P*u(m+SMAX)
            endif
         enddo
         if (h == 1) vp(:,ff) = u
         if (h == 2) vm(:,ff) = u
      enddo
   enddo
   endsubroutine fill_stencils
endprogram test_weno_seth_parity
