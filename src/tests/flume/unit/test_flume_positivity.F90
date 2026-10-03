!< Unit test PV-0 of the FLUME positivity limiter (issue #47, D-9, M3-P5a).
program test_flume_positivity
!< Unit test PV-0 of the FLUME positivity limiter (issue #47, D-9, M3-P5a).
!<
!< **Why this test exists.** The limiter's guarantee is a statement about one forward-Euler step: whatever the high-order
!< face fluxes, the blended update of a cell whose first-order backbone is admissible keeps a positive density and
!< internal energy. A sign error in an antidiffusive difference, a face taking the wrong cell's factor, a missing corner
!< or a source left out of the backbone all still run, and fail only on the hard problem that needs them; here they
!< fail on random data.
!<
!< **What it pins**, for the three models that accept the limiter (Euler, MHD without divergence control, MHD with
!< EGLM), on `N` random trials: one block of 4^3 cells with 3 ghost layers holding random admissible states (density and
!< pressure log-uniform over 4-5 decades, velocity within 2 sound speeds, field up to 10 sqrt(rho) a, so beta down to
!< 1e-2; EGLM psi of the field's size), `dt = 0.4 dx / (3 sigma_max)` (the largest Wu speed over the faces for MHD),
!< and high-order fluxes equal to the backbone ones
!< plus random perturbations twice the flux scale. Through the CPU kernels (`compute_positivity_factors`, the ghost
!< factors left at 1, `blend_positivity_fluxes`, the flux difference, and for EGLM the damping and
!< `add_eglm_sources_limited`):
!< 1. every interior cell with an admissible backbone has positive density and internal energy after the update (count of
!<    failures, must be 0);
!< 2. the limiter is needed and acts: without it some updates are inadmissible, and some cells get `Lambda < 1` (both
!<    must be > 0, else the test exercises nothing);
!< 3. inadmissible backbones (count, reported: the backbone's own CFL-type condition is not part of the guarantee);
!< 4. the relative floor (issue #47, M3-P5c): the same cells keep density and internal energy above
!<    `POSITIVITY_LIMITER_KAPPA` times those of their backbone update (count of failures, must be 0; the internal energy
!<    within the round-off of the total energy);
!< 5. non-finite high-order fluxes (issue #47, M3-P6): each trial sets one face flux component to NaN and one to
!<    +infinity; the kernel flags their cells (count, must be > 0) and after the blend no face flux holds a non-finite
!<    value (count, must be 0): those faces take the backbone flux.
!< 6. the ghost positivity blend (issue #50, D2, `blend_inadmissible_ghosts`): on each trial's block three face ghosts
!<    are made inadmissible (negative internal energy, a NaN momentum, negative density) and one edge ghost too; the
!<    blend must count exactly the three face ghosts, take each above `POSITIVITY_LIMITER_KAPPA` times the density and
!<    internal energy of its interior anchor (within the round-off of the total energy), and leave every other value,
!<    the edge ghost included, bitwise unchanged (count of failures, must be 0).
!<
!< **Not pinned: the source term of the corners** (`dt (s_hi - s_lo)`, EGLM). Measured (M3-P5a): removing it is not
!< detected, and a case with the high-order fluxes equal to the backbone ones never needs limiting, not even on cold,
!< strongly magnetised states. On random cell data the field jumps that make `div B` large also make the backbone's
!< Lax-Friedrichs dissipation large, which heats the cell far more than `dt (div B)(u.B)` cools it; on smooth data the
!< two source orders differ at O(h^2). The source path is exercised by the blast (EGLM, splitting scheme).

use :: adam_flume_cpu_euler_kernels,    only : blend_euler=>blend_positivity_fluxes,                &
                                               factors_euler=>compute_positivity_factors,           &
                                               ghosts_euler=>blend_inadmissible_ghosts
use :: adam_flume_cpu_mhd_kernels,      only : blend_mhd=>blend_positivity_fluxes,                  &
                                               factors_mhd=>compute_positivity_factors,             &
                                               ghosts_mhd=>blend_inadmissible_ghosts
use :: adam_flume_cpu_mhd_eglm_kernels, only : add_eglm_sources_limited, add_glm_damping,           &
                                               blend_eglm=>blend_positivity_fluxes,                 &
                                               ghosts_eglm=>blend_inadmissible_ghosts,              &
                                               factors_eglm=>compute_positivity_factors
use :: adam_flume_euler_library,        only : compute_riemann_llf, conservative_to_auxiliary
use :: adam_flume_mhd_library,          only : mhd_conservative_to_auxiliary, mhd_eglm_conservative_to_auxiliary, &
                                               mhd_fast_speed
use :: adam_flume_mhd_riemann_library,  only : mhd_backbone_flux, mhd_eglm_backbone_flux
use :: adam_flume_parameters,           only : IA_A, IA_BX, IA_BZ, IA_R, IA_U, IQ_BX, IQ_BY, IQ_BZ, IQ_PSI, IQ_R,      &
                                               IQ_RE, IQ_RU, IQ_RV, IQ_RW, NV_AUX, NV_AUX_MHD, NV_EULER, NV_MHD,       &
                                               NV_MHD_EGLM, POSITIVITY_LIMITER_KAPPA
use :: penf,                            only : I4P, I8P, R8P, str
use, intrinsic :: ieee_arithmetic,      only : ieee_positive_inf, ieee_quiet_nan, ieee_value

implicit none

integer(I4P), parameter :: N=200_I4P           !< Random trials per model.
integer(I4P), parameter :: NC=4_I4P            !< Cells per direction.
integer(I4P), parameter :: NGC=3_I4P           !< Ghost cells.
integer(I4P), parameter :: HS=3_I4P            !< Half stencil of the high-order EGLM sources.
real(R8P),    parameter :: GAMMA=5._R8P/3._R8P !< Specific heats ratio.
real(R8P),    parameter :: CH=4._R8P           !< GLM cleaning speed (EGLM).
real(R8P),    parameter :: DAMPING=0.5_R8P     !< GLM damping rate (EGLM).
real(R8P),    parameter :: DX=0.25_R8P         !< Cell size (all directions).
integer(I4P)            :: fails(3)            !< Inadmissible updates per model.
integer(I4P)            :: limited(3)          !< Cells with Lambda < 1 per model.
integer(I4P)            :: bad(3)              !< Inadmissible backbones per model.
integer(I4P)            :: unlim(3)            !< Inadmissible updates without the limiter per model.
integer(I4P)            :: below(3)            !< Updates below the relative floor per model.
integer(I4P)            :: nonf(3)             !< Cells flagged for a non-finite high-order flux per model.
integer(I4P)            :: leak(3)             !< Non-finite face flux values after the blend per model.
integer(I4P)            :: gfail(3)            !< Ghost blend failures per model.
integer(I4P)            :: seed(64)            !< Random generator seed.
integer(I4P)            :: ns, n_, m           !< Seed size, counters.
logical                 :: test_passed         !< Aggregate pass flag.
character(len=5)        :: names(3)            !< Model names.

names = ['Euler', 'MHD  ', 'EGLM ']
call random_seed(size=ns)
if (ns > size(seed)) error stop 'random seed larger than expected'
seed = [(20261002_I4P + n_, n_=1, size(seed))]
call random_seed(put=seed(1:ns))
fails = 0 ; limited = 0 ; bad = 0 ; unlim = 0 ; below = 0 ; nonf = 0 ; leak = 0 ; gfail = 0
do n_=1, N
   do m=1, 3
      call trial(model=m, fails=fails(m), limited=limited(m), bad=bad(m), unlim=unlim(m), below=below(m), &
                 nonf=nonf(m), leak=leak(m), gfail=gfail(m))
   enddo
enddo
test_passed = .true.
do m=1, 3
   print '(A)', names(m)//': inadmissible updates '//trim(str(fails(m)))//', limited cells '//trim(str(limited(m)))// &
                ', inadmissible backbones '//trim(str(bad(m)))//', unlimited inadmissible '//        &
                trim(str(unlim(m)))//', below the relative floor '//trim(str(below(m)))//      &
                ', non-finite flux cells '//trim(str(nonf(m)))//', non-finite fluxes after the blend '//  &
                trim(str(leak(m)))//', ghost blend failures '//trim(str(gfail(m)))//' ('//trim(str(N))//      &
                ' trials of '//trim(str(NC**3))//' cells)'
   test_passed = test_passed .and. fails(m) == 0_I4P .and. limited(m) > 0_I4P .and. unlim(m) > 0_I4P .and. &
                 below(m) == 0_I4P .and. nonf(m) > 0_I4P .and. leak(m) == 0_I4P .and. gfail(m) == 0_I4P
enddo
if (test_passed) then
   print '(A)', 'TEST PASSED: flume positivity limiter'
else
   print '(A)', 'TEST FAILED: flume positivity limiter'
   error stop 1
endif

contains
   function uniform() result(x)
   !< Return a pseudo-random number in [0, 1) from the intrinsic generator (seeded once, deterministic).
   real(R8P) :: x !< Random number.

   call random_number(x)
   endfunction uniform

   pure function count_nonfinite(f) result(n)
   !< Return the number of non-finite values of a flux array (bit pattern: all the exponent bits set).
   real(R8P), intent(in) :: f(:,:,:,:,:) !< Face fluxes.
   integer(I4P)          :: n            !< Non-finite values.

   n = count(iand(ishft(transfer(f, 0_I8P, size(f)), -52), 2047_I8P) == 2047_I8P)
   endfunction count_nonfinite

   pure function energy(q) result(e)
   !< Return the internal energy per unit volume of a state (the model from its size: 5, 8 or 9 variables).
   real(R8P), intent(in) :: q(:) !< Conservative variables.
   real(R8P)             :: e    !< Internal energy.

   e = q(IQ_RE) - 0.5_R8P * (q(IQ_RU)**2 + q(IQ_RV)**2 + q(IQ_RW)**2) / q(IQ_R)
   if (size(q) >= NV_MHD) e = e - 0.5_R8P * (q(IQ_BX)**2 + q(IQ_BY)**2 + q(IQ_BZ)**2)
   if (size(q) == NV_MHD_EGLM) e = e - 0.5_R8P * q(IQ_PSI)**2
   endfunction energy

   subroutine trial(model, fails, limited, bad, unlim, below, nonf, leak, gfail)
   !< One random trial of a model (1 Euler, 2 MHD, 3 EGLM).
   integer(I4P), intent(in)    :: model                     !< Model.
   integer(I4P), intent(inout) :: fails, limited, bad       !< Counters.
   integer(I4P), intent(inout) :: unlim                     !< Inadmissible updates without the limiter.
   integer(I4P), intent(inout) :: below                     !< Updates below the relative floor.
   integer(I4P), intent(inout) :: nonf                      !< Cells flagged for a non-finite high-order flux.
   integer(I4P), intent(inout) :: leak                      !< Non-finite face flux values after the blend.
   integer(I4P), intent(inout) :: gfail                     !< Ghost blend failures.
   integer(I4P)                :: nn                        !< Flagged cells of the trial.
   real(R8P), allocatable      :: fbx(:,:,:,:,:), fby(:,:,:,:,:), fbz(:,:,:,:,:) !< Backbone face fluxes.
   real(R8P), allocatable      :: q(:,:,:,:,:)              !< State.
   real(R8P), allocatable      :: qa(:,:,:,:,:)             !< Auxiliary variables.
   real(R8P), allocatable      :: dq(:,:,:,:,:)             !< Residuals.
   real(R8P), allocatable      :: lam(:,:,:,:,:)            !< Cell factors.
   real(R8P), allocatable      :: flx(:,:,:,:,:), fly(:,:,:,:,:), flz(:,:,:,:,:) !< Face fluxes.
   real(R8P)                   :: prim(9)                   !< Primitive state.
   real(R8P)                   :: smax, a, dt               !< Speed bound, sound speed, step.
   real(R8P)                   :: dxyz(3,1)                 !< Space steps.
   logical                     :: is_null(3)                !< Null directions.
   integer(I4P)                :: nv, na                    !< Variables numbers.
   integer(I4P)                :: i, j, k, c, d, nb, nl     !< Counters.

   select case(model)
   case(1)
      nv = NV_EULER
   case(2)
      nv = NV_MHD
   case default
      nv = NV_MHD_EGLM
   endselect
   na = merge(NV_AUX, NV_AUX_MHD, model == 1)
   dxyz = DX
   is_null = .false.
   allocate(q(nv,1-NGC:NC+NGC,1-NGC:NC+NGC,1-NGC:NC+NGC,1), qa(na,1-NGC:NC+NGC,1-NGC:NC+NGC,1-NGC:NC+NGC,1))
   allocate(dq(nv,1-NGC:NC+NGC,1-NGC:NC+NGC,1-NGC:NC+NGC,1), lam(nv,1-NGC:NC+NGC,1-NGC:NC+NGC,1-NGC:NC+NGC,1))
   allocate(flx(nv,0:NC,1:NC,1:NC,1), fly(nv,1:NC,0:NC,1:NC,1), flz(nv,1:NC,1:NC,0:NC,1))
   allocate(fbx(nv,0:NC,1:NC,1:NC,1), fby(nv,1:NC,0:NC,1:NC,1), fbz(nv,1:NC,1:NC,0:NC,1))
   smax = 0._R8P
   do k=1-NGC, NC+NGC
      do j=1-NGC, NC+NGC
         do i=1-NGC, NC+NGC
            prim = 0._R8P
            prim(1) = 10._R8P**(-3._R8P + 4._R8P * uniform())
            prim(5) = 10._R8P**(-4._R8P + 5._R8P * uniform())
            a = sqrt(GAMMA * prim(5) / prim(1))
            do c=2, 4
               prim(c) = 2._R8P * a * (2._R8P * uniform() - 1._R8P)
            enddo
            if (model > 1) then
               do c=6, 8
                  prim(c) = 10._R8P * uniform() * sqrt(prim(1)) * a * (2._R8P * uniform() - 1._R8P)
               enddo
            endif
            if (model == 3) prim(9) = 10._R8P * uniform() * sqrt(prim(1)) * a * (2._R8P * uniform() - 1._R8P)
            q(IQ_R,i,j,k,1)  = prim(1)
            q(IQ_RU,i,j,k,1) = prim(1) * prim(2)
            q(IQ_RV,i,j,k,1) = prim(1) * prim(3)
            q(IQ_RW,i,j,k,1) = prim(1) * prim(4)
            q(IQ_RE,i,j,k,1) = prim(5) / (GAMMA - 1._R8P) + 0.5_R8P * prim(1) * (prim(2)**2 + prim(3)**2 + prim(4)**2)
            select case(model)
            case(1)
               call conservative_to_auxiliary(gamma=GAMMA, R=1._R8P, q=q(:,i,j,k,1), qa=qa(:,i,j,k,1))
               smax = max(smax, sqrt(prim(2)**2 + prim(3)**2 + prim(4)**2) + qa(IA_A,i,j,k,1))
            case(2, 3)
               q(IQ_BX:IQ_BZ,i,j,k,1) = prim(6:8)
               q(IQ_RE,i,j,k,1) = q(IQ_RE,i,j,k,1) + 0.5_R8P * (prim(6)**2 + prim(7)**2 + prim(8)**2)
               if (model == 3) then
                  q(IQ_PSI,i,j,k,1) = prim(9)
                  q(IQ_RE,i,j,k,1) = q(IQ_RE,i,j,k,1) + 0.5_R8P * prim(9)**2
                  call mhd_eglm_conservative_to_auxiliary(gamma=GAMMA, R=1._R8P, q=q(:,i,j,k,1), qa=qa(:,i,j,k,1))
               else
                  call mhd_conservative_to_auxiliary(gamma=GAMMA, R=1._R8P, q=q(:,i,j,k,1), qa=qa(:,i,j,k,1))
               endif
               do d=1, 3
                  smax = max(smax, abs(qa(IA_U+d-1,i,j,k,1)) + mhd_fast_speed(d=d, qa=qa(:,i,j,k,1)))
               enddo
            endselect
         enddo
      enddo
   enddo
   ! the largest Wu speed over the faces of the block (MHD): the fastest wave plus the field-jump term
   if (model > 1) then
      do k=0, NC+1
         do j=0, NC+1
            do i=0, NC+1
               do d=1, 3
                  smax = max(smax, wu(d=d, qaL=qa(:,i,j,k,1),                                                &
                                      qaR=qa(:,i+merge(1,0,d==1),j+merge(1,0,d==2),k+merge(1,0,d==3),1)))
               enddo
            enddo
         enddo
      enddo
   endif
   if (model == 3) smax = max(smax, CH)
   dt = 0.4_R8P * DX / (3._R8P * smax)
   do k=1, NC ; do j=1, NC ; do i=0, NC
      call backbone(model=model, d=1, qL=q(:,i,j,k,1), qR=q(:,i+1,j,k,1), f=flx(:,i,j,k,1))
      fbx(:,i,j,k,1) = flx(:,i,j,k,1)
      call perturb(f=flx(:,i,j,k,1), q=q(:,i,j,k,1), sigma=smax)
   enddo ; enddo ; enddo
   do k=1, NC ; do j=0, NC ; do i=1, NC
      call backbone(model=model, d=2, qL=q(:,i,j,k,1), qR=q(:,i,j+1,k,1), f=fly(:,i,j,k,1))
      fby(:,i,j,k,1) = fly(:,i,j,k,1)
      call perturb(f=fly(:,i,j,k,1), q=q(:,i,j,k,1), sigma=smax)
   enddo ; enddo ; enddo
   do k=0, NC ; do j=1, NC ; do i=1, NC
      call backbone(model=model, d=3, qL=q(:,i,j,k,1), qR=q(:,i,j,k+1,1), f=flz(:,i,j,k,1))
      fbz(:,i,j,k,1) = flz(:,i,j,k,1)
      call perturb(f=flz(:,i,j,k,1), q=q(:,i,j,k,1), sigma=smax)
   enddo ; enddo ; enddo
   call ghost_check(model, q, is_null, gfail)
   lam = 1._R8P
   call update(model, q, qa, lam, flx, fly, flz, dt, dxyz, is_null, unlim)
   flx(IQ_R,2,2,2,1) = ieee_value(1._R8P, ieee_quiet_nan)
   fly(IQ_RE,3,1,3,1) = ieee_value(1._R8P, ieee_positive_inf)
   select case(model)
   case(1)
      call factors_euler(ni=NC, nj=NC, nk=NC, ngc=NGC, blocks_number=1, gamma=GAMMA, ch=CH, damping=0._R8P, hs=HS, &
                         dt=dt, dxyz=dxyz, is_null=is_null, q=q, q_aux=qa, flx=flx, fly=fly, flz=flz, lam=lam, bad=nb, &
                         nonfinite=nn)
      call blend_euler(d=1, di=1, dj=0, dk=0, ni=NC, nj=NC, nk=NC, ngc=NGC, blocks_number=1, gamma=GAMMA, ch=CH, q=q, &
                       lam=lam, fl=flx, limited=nl)
      call blend_euler(d=2, di=0, dj=1, dk=0, ni=NC, nj=NC, nk=NC, ngc=NGC, blocks_number=1, gamma=GAMMA, ch=CH, q=q, &
                       lam=lam, fl=fly, limited=nl)
      call blend_euler(d=3, di=0, dj=0, dk=1, ni=NC, nj=NC, nk=NC, ngc=NGC, blocks_number=1, gamma=GAMMA, ch=CH, q=q, &
                       lam=lam, fl=flz, limited=nl)
   case(2)
      call factors_mhd(ni=NC, nj=NC, nk=NC, ngc=NGC, blocks_number=1, gamma=GAMMA, ch=CH, damping=0._R8P, hs=HS, &
                       dt=dt, dxyz=dxyz, is_null=is_null, q=q, q_aux=qa, flx=flx, fly=fly, flz=flz, lam=lam, bad=nb, &
                         nonfinite=nn)
      call blend_mhd(d=1, di=1, dj=0, dk=0, ni=NC, nj=NC, nk=NC, ngc=NGC, blocks_number=1, gamma=GAMMA, ch=CH, q=q, &
                     lam=lam, fl=flx, limited=nl)
      call blend_mhd(d=2, di=0, dj=1, dk=0, ni=NC, nj=NC, nk=NC, ngc=NGC, blocks_number=1, gamma=GAMMA, ch=CH, q=q, &
                     lam=lam, fl=fly, limited=nl)
      call blend_mhd(d=3, di=0, dj=0, dk=1, ni=NC, nj=NC, nk=NC, ngc=NGC, blocks_number=1, gamma=GAMMA, ch=CH, q=q, &
                     lam=lam, fl=flz, limited=nl)
   case(3)
      call factors_eglm(ni=NC, nj=NC, nk=NC, ngc=NGC, blocks_number=1, gamma=GAMMA, ch=CH, damping=DAMPING, hs=HS, &
                        dt=dt, dxyz=dxyz, is_null=is_null, q=q, q_aux=qa, flx=flx, fly=fly, flz=flz, lam=lam, bad=nb, &
                         nonfinite=nn)
      call blend_eglm(d=1, di=1, dj=0, dk=0, ni=NC, nj=NC, nk=NC, ngc=NGC, blocks_number=1, gamma=GAMMA, ch=CH, q=q, &
                      lam=lam, fl=flx, limited=nl)
      call blend_eglm(d=2, di=0, dj=1, dk=0, ni=NC, nj=NC, nk=NC, ngc=NGC, blocks_number=1, gamma=GAMMA, ch=CH, q=q, &
                      lam=lam, fl=fly, limited=nl)
      call blend_eglm(d=3, di=0, dj=0, dk=1, ni=NC, nj=NC, nk=NC, ngc=NGC, blocks_number=1, gamma=GAMMA, ch=CH, q=q, &
                      lam=lam, fl=flz, limited=nl)
   endselect
   bad = bad + nb
   nonf = nonf + nn
   leak = leak + count_nonfinite(flx) + count_nonfinite(fly) + count_nonfinite(flz)
   call update(model, q, qa, lam, flx, fly, flz, dt, dxyz, is_null, fails)
   call floor_check(model, q, qa, lam, flx, fly, flz, fbx, fby, fbz, dt, dxyz, is_null, below)
   do k=1, NC
      do j=1, NC
         do i=1, NC
            if (lam(1,i,j,k,1) < 1._R8P) limited = limited + 1_I4P
         enddo
      enddo
   enddo
   endsubroutine trial

   subroutine ghost_check(model, q, is_null, gfail)
   !< Check 6: corrupt three face ghosts and one edge ghost of a copy of `q`, blend, and count the failures.
   integer(I4P), intent(in)    :: model                     !< Model.
   real(R8P),    intent(in)    :: q(:,1-NGC:,1-NGC:,1-NGC:,:) !< Random admissible state (ghosts included).
   logical,      intent(in)    :: is_null(3)                !< Null directions.
   integer(I4P), intent(inout) :: gfail                     !< Failures.
   real(R8P), allocatable      :: qg(:,:,:,:,:), q0(:,:,:,:,:) !< Blended copy, corrupted reference.
   integer(I4P)                :: g(3,3), a(3,3)            !< Corrupted face ghosts and their interior anchors.
   integer(I4P)                :: blended, n                !< Ghosts blended, counter.
   real(R8P)                   :: e                         !< Internal energy.

   g = reshape([0, 2, 3, NC+2, 1, 1, 2, 0, 4], [3, 3])
   a = reshape([1, 2, 3, NC, 1, 1, 2, 1, 4], [3, 3])
   q0 = q
   e = energy(q0(:,g(1,1),g(2,1),g(3,1),1))
   q0(IQ_RE,g(1,1),g(2,1),g(3,1),1) = q0(IQ_RE,g(1,1),g(2,1),g(3,1),1) - 2._R8P * e
   q0(IQ_RU,g(1,2),g(2,2),g(3,2),1) = ieee_value(1._R8P, ieee_quiet_nan)
   q0(IQ_R,g(1,3),g(2,3),g(3,3),1) = -q0(IQ_R,g(1,3),g(2,3),g(3,3),1)
   q0(IQ_R,0,0,2,1) = -q0(IQ_R,0,0,2,1)
   qg = q0
   select case(model)
   case(1)
      call ghosts_euler(ni=NC, nj=NC, nk=NC, ngc=NGC, blocks_number=1, is_null=is_null, q=qg, blended=blended)
   case(2)
      call ghosts_mhd(ni=NC, nj=NC, nk=NC, ngc=NGC, blocks_number=1, is_null=is_null, q=qg, blended=blended)
   case default
      call ghosts_eglm(ni=NC, nj=NC, nk=NC, ngc=NGC, blocks_number=1, is_null=is_null, q=qg, blended=blended)
   endselect
   if (blended /= 3_I4P) gfail = gfail + 1_I4P
   do n=1, 3
      associate(u=>qg(:,g(1,n),g(2,n),g(3,n),1), v=>q0(:,a(1,n),a(2,n),a(3,n),1))
      if (.not.(u(IQ_R) >= POSITIVITY_LIMITER_KAPPA * v(IQ_R) * (1._R8P - 1.e-12_R8P))) gfail = gfail + 1_I4P
      if (.not.(energy(u) >= POSITIVITY_LIMITER_KAPPA * energy(v) - 1.e-12_R8P * abs(v(IQ_RE)))) gfail = gfail + 1_I4P
      endassociate
      qg(:,g(1,n),g(2,n),g(3,n),1) = q0(:,g(1,n),g(2,n),g(3,n),1)
   enddo
   if (any(transfer(qg, 0_I8P, size(qg)) /= transfer(q0, 0_I8P, size(q0)))) gfail = gfail + 1_I4P
   endsubroutine ghost_check

   subroutine backbone(model, d, qL, qR, f)
   !< The backbone flux of a model (1 Euler, 2 MHD, 3 EGLM).
   integer(I4P), intent(in)  :: model        !< Model.
   integer(I4P), intent(in)  :: d            !< Direction.
   real(R8P),    intent(in)  :: qL(:), qR(:) !< States.
   real(R8P),    intent(out) :: f(:)         !< Flux.

   select case(model)
   case(1)
      call compute_riemann_llf(gamma=GAMMA, d=d, qL=qL, qR=qR, f=f)
   case(2)
      call mhd_backbone_flux(gamma=GAMMA, d=d, qL=qL, qR=qR, f=f)
   case(3)
      call mhd_eglm_backbone_flux(ch=CH, gamma=GAMMA, d=d, qL=qL, qR=qR, f=f)
   endselect
   endsubroutine backbone

   pure function wu(d, qaL, qaR) result(sigma)
   !< Return the Wu (2018) speed of two MHD states (as `adam_flume_mhd_riemann_library`, which keeps it private).
   integer(I4P), intent(in) :: d               !< Direction.
   real(R8P),    intent(in) :: qaL(:), qaR(:)  !< Auxiliary variables.
   real(R8P)                :: sigma           !< Speed.
   real(R8P)                :: srL, srR, un    !< sqrt(rho), weighted normal velocity.

   srL = sqrt(qaL(IA_R))
   srR = sqrt(qaR(IA_R))
   un  = (srL * qaL(IA_U+d-1) + srR * qaR(IA_U+d-1)) / (srL + srR)
   sigma = max(abs(qaL(IA_U+d-1)) + mhd_fast_speed(d=d, qa=qaL), abs(qaR(IA_U+d-1)) + mhd_fast_speed(d=d, qa=qaR), &
               abs(un) + max(mhd_fast_speed(d=d, qa=qaL), mhd_fast_speed(d=d, qa=qaR))) +                       &
           sqrt(sum((qaL(IA_BX:IA_BZ) - qaR(IA_BX:IA_BZ))**2)) / (srL + srR)
   endfunction wu

   subroutine perturb(f, q, sigma)
   !< Add to a backbone flux a random perturbation of twice the flux scale per component.
   real(R8P), intent(inout) :: f(:)  !< Flux.
   real(R8P), intent(in)    :: q(:)  !< A state of the face (scale).
   real(R8P), intent(in)    :: sigma !< Speed scale.
   integer(I4P)             :: v     !< Counter.

   do v=1, size(f)
      f(v) = f(v) + 2._R8P * (2._R8P * uniform() - 1._R8P) * (abs(f(v)) + sigma * abs(q(v)))
   enddo
   endsubroutine perturb

   subroutine update(model, q, qa, lam, flx, fly, flz, dt, dxyz, is_null, count)
   !< Count the interior cells (with an admissible backbone) whose forward-Euler update with the given fluxes and
   !< factors is inadmissible (the EGLM damping and the sources limited by `lam`).
   integer(I4P), intent(in)    :: model                          !< Model.
   real(R8P),    intent(in)    :: q(1:,1-NGC:,1-NGC:,1-NGC:,1:)   !< State.
   real(R8P),    intent(in)    :: qa(1:,1-NGC:,1-NGC:,1-NGC:,1:)  !< Auxiliary variables.
   real(R8P),    intent(in)    :: lam(1:,1-NGC:,1-NGC:,1-NGC:,1:) !< Cell factors.
   real(R8P),    intent(in)    :: flx(1:,0:,1:,1:,1:)             !< X-face fluxes.
   real(R8P),    intent(in)    :: fly(1:,1:,0:,1:,1:)             !< Y-face fluxes.
   real(R8P),    intent(in)    :: flz(1:,1:,1:,0:,1:)             !< Z-face fluxes.
   real(R8P),    intent(in)    :: dt                              !< Time step.
   real(R8P),    intent(in)    :: dxyz(3,1)                       !< Space steps.
   logical,      intent(in)    :: is_null(3)                      !< Null directions.
   integer(I4P), intent(inout) :: count                           !< Counter.
   real(R8P)                   :: dq(size(q,1),1-NGC:NC+NGC,1-NGC:NC+NGC,1-NGC:NC+NGC,1) !< Residuals.
   integer(I4P)                :: i, j, k                         !< Counters.

   call residual(model, q, qa, lam, flx, fly, flz, dxyz, is_null, dq)
   do k=1, NC
      do j=1, NC
         do i=1, NC
            if (lam(1,i,j,k,1) == 0._R8P) cycle ! an inadmissible backbone (counted in bad) carries no guarantee
            if (.not.(q(IQ_R,i,j,k,1) + dt * dq(IQ_R,i,j,k,1) > 0._R8P)) then
               count = count + 1_I4P
            elseif (.not.(energy(q(:,i,j,k,1) + dt * dq(:,i,j,k,1)) > 0._R8P)) then
               count = count + 1_I4P
            endif
         enddo
      enddo
   enddo
   endsubroutine update

   subroutine residual(model, q, qa, lam, flx, fly, flz, dxyz, is_null, dq)
   !< Return the residuals of the given fluxes and factors (the EGLM damping and the sources limited by `lam`).
   integer(I4P), intent(in)  :: model                           !< Model.
   real(R8P),    intent(in)  :: q(1:,1-NGC:,1-NGC:,1-NGC:,1:)    !< State.
   real(R8P),    intent(in)  :: qa(1:,1-NGC:,1-NGC:,1-NGC:,1:)   !< Auxiliary variables.
   real(R8P),    intent(in)  :: lam(1:,1-NGC:,1-NGC:,1-NGC:,1:)  !< Cell factors.
   real(R8P),    intent(in)  :: flx(1:,0:,1:,1:,1:)              !< X-face fluxes.
   real(R8P),    intent(in)  :: fly(1:,1:,0:,1:,1:)              !< Y-face fluxes.
   real(R8P),    intent(in)  :: flz(1:,1:,1:,0:,1:)              !< Z-face fluxes.
   real(R8P),    intent(in)  :: dxyz(3,1)                        !< Space steps.
   logical,      intent(in)  :: is_null(3)                       !< Null directions.
   real(R8P),    intent(out) :: dq(1:,1-NGC:,1-NGC:,1-NGC:,1:)   !< Residuals.
   integer(I4P)              :: i, j, k                          !< Counters.

   dq = 0._R8P
   do k=1, NC
      do j=1, NC
         do i=1, NC
            dq(:,i,j,k,1) = -(flx(:,i,j,k,1) - flx(:,i-1,j,k,1) + fly(:,i,j,k,1) - fly(:,i,j-1,k,1) + &
                              flz(:,i,j,k,1) - flz(:,i,j,k-1,1)) / DX
         enddo
      enddo
   enddo
   if (model == 3) then
      call add_glm_damping(ni=NC, nj=NC, nk=NC, ngc=NGC, blocks_number=1, damping=DAMPING, q=q, dq=dq)
      call add_eglm_sources_limited(ni=NC, nj=NC, nk=NC, ngc=NGC, blocks_number=1, hs=HS, dxyz=dxyz, &
                                    is_null=is_null, q=q, q_aux=qa, lam=lam, dq=dq)
   endif
   endsubroutine residual

   subroutine floor_check(model, q, qa, lam, flx, fly, flz, fbx, fby, fbz, dt, dxyz, is_null, count)
   !< Count the interior cells (with an admissible backbone) whose limited update has a density or an internal energy
   !< below `POSITIVITY_LIMITER_KAPPA` times that of the backbone update (backbone fluxes, zero factors: the EGLM sources
   !< at second order); the internal energy within the round-off of the total energy.
   integer(I4P), intent(in)    :: model                            !< Model.
   real(R8P),    intent(in)    :: q(1:,1-NGC:,1-NGC:,1-NGC:,1:)     !< State.
   real(R8P),    intent(in)    :: qa(1:,1-NGC:,1-NGC:,1-NGC:,1:)    !< Auxiliary variables.
   real(R8P),    intent(in)    :: lam(1:,1-NGC:,1-NGC:,1-NGC:,1:)   !< Cell factors.
   real(R8P),    intent(in)    :: flx(1:,0:,1:,1:,1:)               !< X-face limited fluxes.
   real(R8P),    intent(in)    :: fly(1:,1:,0:,1:,1:)               !< Y-face limited fluxes.
   real(R8P),    intent(in)    :: flz(1:,1:,1:,0:,1:)               !< Z-face limited fluxes.
   real(R8P),    intent(in)    :: fbx(1:,0:,1:,1:,1:)               !< X-face backbone fluxes.
   real(R8P),    intent(in)    :: fby(1:,1:,0:,1:,1:)               !< Y-face backbone fluxes.
   real(R8P),    intent(in)    :: fbz(1:,1:,1:,0:,1:)               !< Z-face backbone fluxes.
   real(R8P),    intent(in)    :: dt                                !< Time step.
   real(R8P),    intent(in)    :: dxyz(3,1)                         !< Space steps.
   logical,      intent(in)    :: is_null(3)                        !< Null directions.
   integer(I4P), intent(inout) :: count                             !< Counter.
   real(R8P)                   :: dq(size(q,1),1-NGC:NC+NGC,1-NGC:NC+NGC,1-NGC:NC+NGC,1)  !< Limited residuals.
   real(R8P)                   :: dqb(size(q,1),1-NGC:NC+NGC,1-NGC:NC+NGC,1-NGC:NC+NGC,1) !< Backbone residuals.
   real(R8P)                   :: lam0(1,1-NGC:NC+NGC,1-NGC:NC+NGC,1-NGC:NC+NGC,1)        !< Zero factors.
   real(R8P)                   :: qn(size(q,1)), qb(size(q,1))      !< Limited and backbone updates.
   integer(I4P)                :: i, j, k                           !< Counters.

   lam0 = 0._R8P
   call residual(model, q, qa, lam, flx, fly, flz, dxyz, is_null, dq)
   call residual(model, q, qa, lam0, fbx, fby, fbz, dxyz, is_null, dqb)
   do k=1, NC
      do j=1, NC
         do i=1, NC
            if (lam(1,i,j,k,1) == 0._R8P) cycle ! an inadmissible backbone carries no guarantee
            qn = q(:,i,j,k,1) + dt * dq(:,i,j,k,1)
            qb = q(:,i,j,k,1) + dt * dqb(:,i,j,k,1)
            if (qn(IQ_R) < POSITIVITY_LIMITER_KAPPA * qb(IQ_R) * (1._R8P - 1.e-12_R8P)) then
               count = count + 1_I4P
            elseif (energy(qn) < POSITIVITY_LIMITER_KAPPA * energy(qb) - 1.e-11_R8P * abs(qn(IQ_RE))) then
               count = count + 1_I4P
            endif
         enddo
      enddo
   enddo
   endsubroutine floor_check
endprogram test_flume_positivity
