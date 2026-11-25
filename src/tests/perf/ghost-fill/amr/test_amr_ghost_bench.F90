!< AMR ghost-fill benchmark (2:1 neighbourhoods): ADAM's local ghost fill on a real 2:1 grid against descriptor kernels.
program test_amr_ghost_bench
!< AMR ghost-fill benchmark (2:1 neighbourhoods, 2026-11-25): ADAM's local ghost fill on a real 2:1 grid, built by
!< ADAM's own tree, grid and maps objects, against kernels that run over ADAM's descriptor table instead of its
!< per-cell rows. The uniform benchmark (`../ghost_fill_bench.F90`) showed the per-cell map costs a factor 3-5.5 on
!< same-level copies; this one asks whether that holds where a row is an interpolation or a restriction stencil.
!<
!< The grid: the unit cube, octree, blocks of b^3 cells (5 variables, ngc 3), refined uniformly to level `base`, then one
!< more level on every block inside the box [0,s]^3 (s = `box`): 2:1 faces, edges and corners around the box. Domain
!< boundaries are not exchanged (their ghosts are boundary-condition rows, not local ones). The density guard of the
!< seam fill (#101, `maps%seam_density_var` = 1) is on, as in FLUME. The field is smooth and positive.
!<
!< ADAM builds a table `maps%local_map_ghost(f, 1:13)`: one descriptor per (receiver block, donor block, region), its
!< ranges, its index shifts and its refinement relation (`portion`: 0 same level, > 0 finer donor, < 0 coarser donor),
!< then expands it into `maps%local_map_ghost_cell`, one row per ghost cell with its seam metadata, deduplicated so the
!< last descriptor covering a ghost wins (#22 F1). Only the rows go to the GPU. Variants:
!<   V0   update_ghost_local_gpu_dev on the rows: today, the real kernel (copy, interpolation and restriction rows mixed
!<        in one launch, a branch per row);
!<   V0s  the same kernel on the rows split by kind, one launch per kind (copy / prolongation / restriction): how much
!<        of today's cost is the mixing (warp divergence between a copy and a 216-point stencil);
!<   D1   one kernel per kind over the descriptors: the work items of each descriptor decoded from its ranges, the seam
!<        metadata computed on the device by ADAM's own pure routines, the value by ADAM's own operators
!<        (interp_seam_ghost_gpu_dev, restrict_seam_high_order_gpu_dev, the mean as in V0);
!<   D2   D1 with the descriptors batched into shape classes (same kind, ranges and shifts: only the blocks differ), the
!<        member innermost so that consecutive threads touch consecutive blocks (the uniform benchmark's V3/V5 layout);
!<   V0ng V0 with the density guard off (timing of the interpolation alone, and the reference of D3);
!<   D3   prolongation per coarse anchor (tricubic only, guard off): the union of its children's 4^3 footprints gathered
!<        once, then every child interpolated by ADAM's seam_interpolate_tricubic on its slice: the reads the 8 children
!<        of one anchor share are made once; bitwise equal to V0ng's prolongation.
!< Overlapping descriptors: an owner mask (int32 per stored cell of every block, built on the host from the descriptor
!< order) gives each ghost the descriptor that wins; D1/D2 write a ghost only from its owner. This is ADAM's
!< last-writer rule, so every variant must equal V0 bit for bit; each is checked after its timing.
!<
!< Usage: test_amr_ghost_bench b base box fill restriction nrep mem_gb [fields]
!<   fields: the block-sized fields realm_object%initialize budgets per block (default 170, its own default): the
!<   block capacity is 0.85 mem_gb / (fields x block bytes); a larger grid needs fewer fields or more memory
!<   fill:        tricubic | quintic | weno6 | quintic-mp        restriction: mean | lagrange6 | weno6
!< Output: the grid and its descriptor and row counts by kind, then per variant the median time of a fill (and per kind
!< for V0s, D1, D2), and the bitwise check.

#include "fundal.H"
use :: adam_common_library,             only : realm_object
use :: adam_parameters,                 only : TO_BE_REFINED, TO_NOT_TOUCH
use :: adam_mpih_global,                only : mpih
use :: adam_fnl_field_object,           only : field_fnl_object
use :: adam_fnl_field_kernels,          only : update_ghost_local_gpu_dev, interp_seam_ghost_gpu_dev, &
                                               restrict_seam_high_order_gpu_dev
use :: adam_fnl_mpih_global,            only : mpih_fnl
use :: adam_seam_interpolation_library, only : SEAM_FILL_TRICUBIC, seam_interpolate_tricubic
use :: adam_seam_interpolation_library, only : SEAM_RESTRICT_MEAN, SEAM_SUB_COINCIDENT, seam_shift_anchor_pos,      &
                                               seam_tricubic_centered_pos, seam_compatible_centered_pos,            &
                                               seam_quintic_centered_pos, seam_meta_pack, seam_restriction_code,   &
                                               seam_restriction_meta_pack
use :: fundal
use :: mpi
use :: penf,                            only : I4P, I8P, R8P, str

implicit none

integer(I4P), parameter :: NV = 5_I4P              !< Variables (FLUME Euler).
integer(I4P), parameter :: KCOPY = 0_I4P           !< Kind: same-level copy.
integer(I4P), parameter :: KREST = 1_I4P           !< Kind: fine->coarse restriction (finer donor).
integer(I4P), parameter :: KPROL = 2_I4P           !< Kind: coarse->fine prolongation (coarser donor).
character(len=13), parameter :: KNAME(0:2) = [character(len=13) :: 'copy', 'restriction', 'prolongation']

type(realm_object)          :: realm                !< The realm.
type(field_fnl_object)      :: field_fnl            !< FNL field helper.
real(R8P), pointer          :: q_gpu(:,:,:,:,:)     !< Device field (b, i, j, k, v).
real(R8P),    allocatable   :: q(:,:,:,:,:)         !< Host field (v, i, j, k, b).
real(R8P),    allocatable   :: q_t(:,:,:,:,:)       !< Host transposed field (b, i, j, k, v), the initial state.
real(R8P),    allocatable   :: q_ref(:,:,:,:,:)     !< V0's result (b, i, j, k, v).
real(R8P),    allocatable   :: q_chk(:,:,:,:,:)     !< A variant's result.
integer(I4P)                :: b, nb, ngc, n, base, nrep, ri, rj, rk, fill, restr, ierr, nfields
real(R8P)                   :: box, mem_gb
character(len=32)           :: arg, fill_name, restr_name
character(:), allocatable   :: filename
! descriptors (from maps%local_map_ghost), host and device
integer(I4P)                :: nf                   !< Descriptors.
integer(I4P), allocatable   :: dbr(:), dbs(:), dkind(:), dlo(:,:), dhi(:,:), ddl(:,:)
integer(I4P), allocatable   :: owner(:,:,:,:)       !< Winning descriptor per stored cell of every block (0: none).
! D1: per kind, the descriptor list and the prefix of its work items
integer(I4P)                :: d1n(0:2)             !< Descriptors per kind.
integer(I4P), allocatable   :: d1l(:), d1p(:)       !< Concatenated lists and prefixes (kind segments).
integer(I4P)                :: d1o(0:3)             !< Segment offsets into d1l.
integer(I4P)                :: d1w(0:2)             !< Work items per kind.
integer(I4P), allocatable   :: d3p(:)               !< D3: prefix of the anchors of the prolongation descriptors.
integer(I4P)                :: d3w                  !< D3: anchors.
! D2: shape classes
integer(I4P)                :: ncls                 !< Classes.
integer(I4P), allocatable   :: crep(:), cnm(:), cmo(:), cmem(:), cpre(:) !< Representative, members, member offset, members, prefix.
integer(I4P)                :: c2o(0:3)             !< Class segment offsets per kind.
integer(I4P)                :: c2w(0:2)             !< Work items per kind.
! V0s: rows split by kind
integer(I8P), pointer       :: rows_k0(:,:) => null(), rows_k1(:,:) => null(), rows_k2(:,:) => null()
integer(I4P)                :: nrows_k(0:2)
real(R8P),    allocatable   :: t(:)

call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)
call mpih_fnl%initialize(do_mpi_init=.false., do_device_init=.true.)
call get_args()
call build_realm()
call describe()
call build_device_state()
allocate(t(nrep))

call time_v0()
call time_v0s()
call time_d(2_I4P)
call time_d(1_I4P)
call time_v0ng_d3()
call mpih%finalize

contains
   subroutine get_args()
   !< b base box fill restriction nrep mem_gb.
   character(len=32) :: a(8)
   integer(I4P)      :: k

   a = [character(len=32) :: '24', '3', '0.5', 'tricubic', 'mean', '30', '8.0', '170']
   do k=1, min(8, command_argument_count())
      call get_command_argument(k, a(k))
   enddo
   read(a(1), *) n ; read(a(2), *) base ; read(a(3), *) box ; fill_name = a(4) ; restr_name = a(5)
   read(a(6), *) nrep ; read(a(7), *) mem_gb ; read(a(8), *) nfields
   endsubroutine get_args

   subroutine build_realm()
   !< Write the input, build the realm, refine uniformly to `base`, refine the box once more, fill the interior.
   integer(I4P) :: i, j, k
   real(R8P)    :: x, y, z, pi

   filename = 'amr-ghost-bench.ini'
   call write_input()
   call realm%initialize(filename=filename, memory_avail=mem_gb, nv=NV, fields_number=nfields)
   ngc = realm%adam%grid%ngc
   allocate(q(NV,1-ngc:n+ngc,1-ngc:n+ngc,1-ngc:n+ngc,1:realm%adam%field%nb))
   q = 1._R8P
   call realm%adam%refine_uniform(refinement_levels=realm%adam%tree%iu_ref_levels, q=q, do_mpi_redistribute=.true., &
                                  do_blocks_reorder=.false.)
   realm%adam%field%refinements_needed = [(TO_NOT_TOUCH, b=1, realm%adam%field%blocks_number)]
   do b=1, realm%adam%field%blocks_number
      if (all(realm%adam%field%emax(:,b) <= box + 1.e-12_R8P)) realm%adam%field%refinements_needed(b) = TO_BE_REFINED
   enddo
   call realm%adam%amr_update(q=q, is_marked_by_field=.true., do_mpi_redistribute=.true., do_blocks_reorder=.false.)
   nb = realm%adam%field%blocks_number
   pi = acos(-1._R8P)
   q = -1._R8P
   do b=1, nb ; do k=1, n ; do j=1, n ; do i=1, n
      x = realm%adam%field%x_cell(i,b) ; y = realm%adam%field%y_cell(j,b) ; z = realm%adam%field%z_cell(k,b)
      q(1,i,j,k,b) = 1._R8P + 0.2_R8P * sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z)  ! density, positive
      q(2,i,j,k,b) = 0.3_R8P * cos(2*pi*x) * sin(2*pi*z)
      q(3,i,j,k,b) = 0.2_R8P * sin(2*pi*y)
      q(4,i,j,k,b) = -0.1_R8P * cos(2*pi*z)
      q(5,i,j,k,b) = 2.5_R8P + 0.1_R8P * cos(2*pi*(x + y + z))
   enddo ; enddo ; enddo ; enddo
   ri = realm%adam%maps%refine_ratio(1) ; rj = realm%adam%maps%refine_ratio(2) ; rk = realm%adam%maps%refine_ratio(3)
   fill = realm%adam%maps%seam_ghost_fill ; restr = realm%adam%maps%seam_restriction
   endsubroutine build_realm

   subroutine write_input()
   !< The sections realm_object%initialize reads; [amr] carries the seam regimes under test.
   integer(I4P) :: u

   open(newunit=u, file=filename, action='write', status='replace')
   write(u, '(A)') '[IO]'
   write(u, '(A)') 'output_basename = amr-ghost-bench'
   write(u, '(A)') 'it_save = 999999'
   write(u, '(A)') 'restart = .false.'
   write(u, '(A)') 'restart_basename = amr-ghost-bench-restart'
   write(u, '(A)') 'restart_save = 999999'
   write(u, '(A)') 'residuals_save = 999999'
   write(u, '(A)') 'save_memory_status = .false.'
   write(u, '(A)') 'save_residual_fields = .false.'
   write(u, '(A)') 'save_curl_fields = .false.'
   write(u, '(A)') 'save_divergence_fields = .false.'
   write(u, '(A)') 'save_gradient_fields = .false.'
   write(u, '(A)') 'save_laplacian_fields = .false.'
   write(u, '(A)') '[grid]'
   write(u, '(A)') 'ni = '//trim(str(n, .true.))
   write(u, '(A)') 'nj = '//trim(str(n, .true.))
   write(u, '(A)') 'nk = '//trim(str(n, .true.))
   write(u, '(A)') 'ngc = 3'
   write(u, '(A)') 'emin_x = 0.0'
   write(u, '(A)') 'emin_y = 0.0'
   write(u, '(A)') 'emin_z = 0.0'
   write(u, '(A)') 'emax_x = 1.0'
   write(u, '(A)') 'emax_y = 1.0'
   write(u, '(A)') 'emax_z = 1.0'
   write(u, '(A)') 'null_x = .false.'
   write(u, '(A)') 'null_y = .false.'
   write(u, '(A)') 'null_z = .false.'
   write(u, '(A)') '[amr]'
   write(u, '(A)') 'max_level = '//trim(str(base + 1, .true.))
   write(u, '(A)') 'ratio = 8'
   write(u, '(A)') 'iu_ref_levels = '//trim(str(base, .true.))
   write(u, '(A)') 'i_prune = 0'
   write(u, '(A)') 'j_prune = 0'
   write(u, '(A)') 'k_prune = 0'
   write(u, '(A)') 'l_prune = -1'
   write(u, '(A)') 'frequency = 999999'
   write(u, '(A)') 'iters = 1'
   write(u, '(A)') 'markers_number = 0'
   write(u, '(A)') 'seam_ghost_fill = '//trim(fill_name)
   write(u, '(A)') 'seam_restriction = '//trim(restr_name)
   write(u, '(A)') '[field]'
   write(u, '(A)') 'nv = '//trim(str(NV, .true.))
   write(u, '(A)') '[runge_kutta]'
   write(u, '(A)') 'scheme = runge-kutta-ssp-33'
   write(u, '(A)') '[weno]'
   write(u, '(A)') 'scheme = weno-u-5'
   write(u, '(A)') 'ror_number = 0'
   write(u, '(A)') 'ror_threshold = 0.9'
   write(u, '(A)') 'ror_vars_number = 0'
   write(u, '(A)') 'enable_ror_stats = .false.'
   write(u, '(A)') 'ib_reduction_extent = 0'
   write(u, '(A)') 'ib_reduced_order = 2'
   write(u, '(A)') '[fdv]'
   write(u, '(A)') 'fdv_scheme = fd'
   write(u, '(A)') 'fdv_order = 2'
   write(u, '(A)') '[solids]'
   write(u, '(A)') 'number = 0'
   write(u, '(A)') 'n_eikonal = 0'
   write(u, '(A)') '[slices]'
   write(u, '(A)') 'slices_number = 0'
   close(u)
   endsubroutine write_input

   subroutine describe()
   !< The grid, the descriptors and rows by kind, the duplicates the dedupe removed; build owner, D1 and D2 tables.
   integer(I8P) :: covered(0:2), cells
   integer(I4P) :: f, k, i, j, kk, ic, jc, kc, rr(3), lo(3), hi(3), dl(3), m, c, w
   integer(I4P), allocatable :: key(:,:), order(:)
   logical      :: same

   associate(lmg => realm%adam%maps%local_map_ghost, cellmap => realm%adam%maps%local_map_ghost_cell)
   nf = size(lmg, dim=1)
   allocate(dbr(nf), dbs(nf), dkind(nf), dlo(3,nf), dhi(3,nf), ddl(3,nf))
   rr = [ri, rj, rk]
   do f=1, nf
      dbr(f) = int(lmg(f,1), I4P) ; dbs(f) = int(lmg(f,2), I4P)
      dlo(:,f) = int(lmg(f,5:7), I4P) ; dhi(:,f) = int(lmg(f,8:10), I4P) ; ddl(:,f) = int(lmg(f,11:13), I4P)
      if (lmg(f,4) == 0) then
         dkind(f) = KCOPY
      elseif (lmg(f,4) > 0) then
         dkind(f) = KREST
      else
         dkind(f) = KPROL
      endif
   enddo
   ! owner mask: the last descriptor covering a receiver cell wins (the dedupe's last-writer rule)
   allocate(owner(nb,1-ngc:n+ngc,1-ngc:n+ngc,1-ngc:n+ngc)) ; owner = 0 ; covered = 0
   do f=1, nf
      do kk=dlo(3,f), dhi(3,f) ; do j=dlo(2,f), dhi(2,f) ; do i=dlo(1,f), dhi(1,f)
         if (dkind(f) /= KPROL) then
            owner(dbr(f),i,j,kk) = f ; covered(dkind(f)) = covered(dkind(f)) + 1
         else
            do kc=0,rk-1 ; do jc=0,rj-1 ; do ic=0,ri-1
               lo = [ri*i + ddl(1,f) + ic, rj*j + ddl(2,f) + jc, rk*kk + ddl(3,f) + kc]
               if (any(lo < 1-ngc) .or. any(lo > n+ngc)) cycle
               owner(dbr(f),lo(1),lo(2),lo(3)) = f ; covered(KPROL) = covered(KPROL) + 1
            enddo ; enddo ; enddo
         endif
      enddo ; enddo ; enddo
   enddo
   ! rows by kind (column 9: 1 copy, 4 prolongation, 8 restriction)
   nrows_k(0) = count(cellmap(:,9) == 1) ; nrows_k(2) = count(cellmap(:,9) == 4) ; nrows_k(1) = count(cellmap(:,9) == 8)
   cells = int(nb, I8P) * int(n, I8P)**3
   print '(A,I3,A,I2,A,F5.3,A,A,A,A)', '# b ', n, ' base level ', base, ' box ', box, ' fill ', trim(fill_name), &
         ' restriction ', trim(restr_name)
   print '(A,I7,A,F9.3,A,I3,A,I3)', '# blocks ', nb, ' cells ', real(cells, R8P) / 1e6_R8P, ' M, levels ', base, '-', base+1
   print '(A,I8,A,3I8,A)', '# descriptors ', nf, ' (copy, restriction, prolongation: ', count(dkind == KCOPY), &
         count(dkind == KREST), count(dkind == KPROL), ')'
   print '(A,I10,A,3I10,A)', '# rows ', size(cellmap, dim=1), ' (copy, restriction, prolongation: ', nrows_k(0), &
         nrows_k(1), nrows_k(2), ')'
   print '(A,3I10,A)', '# cells the descriptors cover (copy, restriction, prolongation) ', covered, &
         ' (minus the rows: the duplicates the dedupe removed)'
   ! D1 tables: per kind the descriptor list and the prefix of its work items (cells; anchors x r^3 for prolongation)
   d1n = [count(dkind == KCOPY), count(dkind == KREST), count(dkind == KPROL)]
   d1o(0) = 0 ; d1o(1) = d1n(0) ; d1o(2) = d1n(0) + d1n(1) ; d1o(3) = nf
   allocate(d1l(nf), d1p(0:nf+2))
   m = 0
   do k=0, 2
      do f=1, nf
         if (dkind(f) == k) then ; m = m + 1 ; d1l(m) = f ; endif
      enddo
   enddo
   do k=0, 2
      w = 0
      d1p(d1o(k) + k) = 0
      do m=d1o(k)+1, d1o(k+1)
         f = d1l(m)
         w = w + product(dhi(:,f) - dlo(:,f) + 1) * merge(ri*rj*rk, 1, k == KPROL)
         d1p(m + k) = w
      enddo
      d1w(k) = w
   enddo
   ! D3 table: the anchors of each prolongation descriptor (no r^3 factor)
   allocate(d3p(0:d1n(2)))
   d3p(0) = 0
   do m=1, d1n(2)
      f = d1l(d1o(2) + m)
      d3p(m) = d3p(m-1) + product(dhi(:,f) - dlo(:,f) + 1)
   enddo
   d3w = d3p(d1n(2))
   ! D2 tables: shape classes (kind, lo, hi, dl), members sorted by receiver block
   allocate(key(10,nf), order(nf))
   do f=1, nf
      key(:,f) = [dkind(f), dlo(:,f), dhi(:,f), ddl(:,f)]
   enddo
   order = [(f, f=1, nf)]
   call sort_by_key(order)
   allocate(crep(nf), cnm(nf), cmo(nf+1), cmem(nf), cpre(0:nf+2))
   ncls = 0 ; c2o = 0
   do m=1, nf
      f = order(m)
      same = .false.
      if (ncls > 0) same = all(key(:,f) == key(:,crep(ncls)))
      if (.not.same) then
         ncls = ncls + 1 ; crep(ncls) = f ; cnm(ncls) = 0 ; cmo(ncls) = m - 1
      endif
      cnm(ncls) = cnm(ncls) + 1 ; cmem(m) = f
   enddo
   cmo(ncls+1) = nf
   do k=0, 2
      c2o(k+1) = c2o(k)
      do c=1, ncls
         if (dkind(crep(c)) == k) c2o(k+1) = c2o(k+1) + 1
      enddo
   enddo
   do k=0, 2
      w = 0
      cpre(c2o(k) + k) = 0
      do c=c2o(k)+1, c2o(k+1)
         f = crep(c)
         w = w + product(dhi(:,f) - dlo(:,f) + 1) * merge(ri*rj*rk, 1, k == KPROL) * cnm(c)
         cpre(c + k) = w
      enddo
      c2w(k) = w
   enddo
   print '(A,I6,A,3I6,A)', '# shape classes ', ncls, ' (copy, restriction, prolongation: ', c2o(1)-c2o(0), &
         c2o(2)-c2o(1), c2o(3)-c2o(2), ')'
   endassociate
   endsubroutine describe

   subroutine sort_by_key(order)
   !< Sort the descriptors by (kind, lo, hi, dl), then receiver block: classes become contiguous runs (merge sort).
   integer(I4P), intent(inout) :: order(:)
   integer(I4P), allocatable   :: tmp(:)
   integer(I4P)                :: width, lo, mid, hi, a, z, k

   allocate(tmp(size(order)))
   width = 1
   do while (width < size(order))
      lo = 1
      do while (lo <= size(order))
         mid = min(lo + width - 1, size(order)) ; hi = min(lo + 2 * width - 1, size(order))
         a = lo ; z = mid + 1 ; k = lo
         do while (k <= hi)
            if (z > hi) then
               tmp(k) = order(a) ; a = a + 1
            elseif (a > mid) then
               tmp(k) = order(z) ; z = z + 1
            elseif (less(order(a), order(z)) .or. .not.less(order(z), order(a))) then
               tmp(k) = order(a) ; a = a + 1
            else
               tmp(k) = order(z) ; z = z + 1
            endif
            k = k + 1
         enddo
         lo = hi + 1
      enddo
      order = tmp
      width = 2 * width
   enddo
   endsubroutine sort_by_key

   logical function less(f, g)
   !< Key order: kind, lo, hi, dl, then the receiver block.
   integer(I4P), intent(in) :: f, g
   integer(I4P)             :: a(11), z(11), k

   a = [dkind(f), dlo(:,f), dhi(:,f), ddl(:,f), dbr(f)]
   z = [dkind(g), dlo(:,g), dhi(:,g), ddl(:,g), dbr(g)]
   less = .false.
   do k=1, 11
      if (a(k) /= z(k)) then
         less = a(k) < z(k) ; return
      endif
   enddo
   endfunction less

   subroutine build_device_state()
   !< The device field and its initial state; the per-kind row maps (V0s); the descriptor tables (D1, D2).
   integer(I8P), allocatable :: tmp(:,:)
   integer(I4P)              :: k, flag(0:2)

   realm%adam%maps%seam_density_var = 1_I4P ! the density guard of the seam fill, on as in FLUME (#101)
   nullify(q_gpu)
   call field_fnl%initialize(grid=realm%adam%grid, field=realm%adam%field, maps=realm%adam%maps, q_gpu=q_gpu)
   allocate(q_t(1:realm%adam%field%nb,1-ngc:n+ngc,1-ngc:n+ngc,1-ngc:n+ngc,1:NV))
   do b=1, realm%adam%field%nb
      q_t(b,:,:,:,:) = reshape(q(:,:,:,:,b), [n+2*ngc, n+2*ngc, n+2*ngc, NV], order=[4, 1, 2, 3])
   enddo
   flag = [1, 8, 4]
   associate(cellmap => realm%adam%maps%local_map_ghost_cell)
   do k=0, 2
      tmp = cellmap(pack([(b, b=1, size(cellmap, dim=1))], cellmap(:,9) == flag(k)), :)
      select case(k)
      case(0) ; call dev_assign_to_device(dst=rows_k0, src=tmp)
      case(1) ; call dev_assign_to_device(dst=rows_k1, src=tmp)
      case(2) ; call dev_assign_to_device(dst=rows_k2, src=tmp)
      endselect
   enddo
   endassociate
   !$acc enter data copyin(dbr, dbs, dkind, dlo, dhi, ddl, owner, d1l, d1p, d1o, crep, cnm, cmo, cmem, cpre, c2o, d3p)
   endsubroutine build_device_state

   subroutine reset()
   !< The initial state on the device: every ghost back to -1.
   call dev_memcpy_to_device(dst=q_gpu, src=q_t)
   endsubroutine reset

   subroutine fetch(dst)
   !< The device field into dst.
   real(R8P), allocatable, intent(inout) :: dst(:,:,:,:,:)

   if (.not.allocated(dst)) allocate(dst, mold=q_t)
   call dev_memcpy_from_device(dst=dst, src=q_gpu)
   endsubroutine fetch

   real(R8P) function wall()
   !< Wall clock, s.
   integer(I8P) :: c, rate

   call system_clock(c, rate)
   wall = real(c, R8P) / real(rate, R8P)
   endfunction wall

   real(R8P) function median(a)
   !< Median of the timings.
   real(R8P), intent(in) :: a(:)
   real(R8P)             :: s(size(a)), x
   integer(I4P)          :: k, l

   s = a
   do k=2, size(s)
      x = s(k) ; l = k - 1
      do while (l >= 1)
         if (s(l) <= x) exit
         s(l+1) = s(l) ; l = l - 1
      enddo
      s(l+1) = x
   enddo
   median = s((size(s) + 1) / 2)
   endfunction median

   subroutine time_v0()
   !< V0: ADAM's kernel on its rows.
   integer(I4P) :: r
   real(R8P)    :: t0

   call reset()
   do r=1, 3 ; call field_fnl%update_ghost_local_gpu(q_gpu=q_gpu) ; enddo
   do r=1, nrep
      !$acc wait
      t0 = wall() ; call field_fnl%update_ghost_local_gpu(q_gpu=q_gpu) ; !$acc wait
      t(r) = wall() - t0
   enddo
   call fetch(q_ref)
   print '(A,F10.3,A)', 'V0  ADAM rows, one launch       median ', median(t)*1e3_R8P, ' ms'
   endsubroutine time_v0

   subroutine run_rows(k)
   !< ADAM's kernel on the rows of one kind.
   integer(I4P), intent(in) :: k

   select case(k)
   case(0) ; if (nrows_k(0) > 0) call update_ghost_local_gpu_dev(ngc=ngc, seam_ghost_fill=fill, seam_restriction=restr,  &
                                           ri=ri, rj=rj, rk=rk, l_map_ghost_cell_gpu=rows_k0, q_gpu=q_gpu, seam_density_var=1)
   case(1) ; if (nrows_k(1) > 0) call update_ghost_local_gpu_dev(ngc=ngc, seam_ghost_fill=fill, seam_restriction=restr,  &
                                           ri=ri, rj=rj, rk=rk, l_map_ghost_cell_gpu=rows_k1, q_gpu=q_gpu, seam_density_var=1)
   case(2) ; if (nrows_k(2) > 0) call update_ghost_local_gpu_dev(ngc=ngc, seam_ghost_fill=fill, seam_restriction=restr,  &
                                           ri=ri, rj=rj, rk=rk, l_map_ghost_cell_gpu=rows_k2, q_gpu=q_gpu, seam_density_var=1)
   endselect
   endsubroutine run_rows

   subroutine time_v0s()
   !< V0s: ADAM's kernel on the rows split by kind, one launch per kind; the per-kind times, then the whole.
   integer(I4P) :: r, k
   real(R8P)    :: t0, tk(0:2)

   call reset()
   do k=0, 2
      do r=1, 2 ; call run_rows(k) ; enddo
      do r=1, nrep
         !$acc wait
         t0 = wall() ; call run_rows(k) ; !$acc wait
         t(r) = wall() - t0
      enddo
      tk(k) = median(t)
   enddo
   do r=1, nrep
      !$acc wait
      t0 = wall() ; call run_rows(0) ; call run_rows(1) ; call run_rows(2) ; !$acc wait
      t(r) = wall() - t0
   enddo
   call fetch(q_chk)
   print '(A,F10.3,A,3F10.3,A,L1)', 'V0s ADAM rows, one launch/kind   median ', median(t)*1e3_R8P, &
         ' ms  (copy, restr, prol ', tk*1e3_R8P, ')  bitwise ', all(q_chk == q_ref)
   endsubroutine time_v0s

   subroutine time_d(variant)
   !< D1 (variant 1) or D2 (variant 2): the descriptor kernels per kind; the per-kind times, then the whole.
   integer(I4P), intent(in) :: variant
   integer(I4P)             :: r, k
   real(R8P)                :: t0, tk(0:2)
   character(len=40)        :: label

   call reset()
   do k=0, 2
      do r=1, 2 ; call run_desc(variant, k, q_gpu) ; enddo
      do r=1, nrep
         !$acc wait
         t0 = wall() ; call run_desc(variant, k, q_gpu) ; !$acc wait
         t(r) = wall() - t0
      enddo
      tk(k) = median(t)
   enddo
   do r=1, nrep
      !$acc wait
      t0 = wall()
      do k=0, 2 ; call run_desc(variant, k, q_gpu) ; enddo
      !$acc wait
      t(r) = wall() - t0
   enddo
   call fetch(q_chk)
   if (variant == 1) then
      label = 'D1  descriptors, one launch/kind'
   else
      label = 'D2  shape classes, one launch/kind'
   endif
   print '(A,A,F10.3,A,3F10.3,A,L1)', label(1:34), ' median ', median(t)*1e3_R8P, ' ms  (copy, restr, prol ', &
         tk*1e3_R8P, ')  bitwise ', all(q_chk == q_ref)
   if (.not.all(q_chk == q_ref)) print '(A,I12)', '   cells that differ: ', count(q_chk /= q_ref)
   endsubroutine time_d

   subroutine time_v0ng_d3()
   !< V0ng (guard off): the whole fill and its prolongation rows; D3 (tricubic only): the prolongation per anchor,
   !< checked bit for bit against V0ng's prolongation.
   integer(I4P) :: r
   real(R8P)    :: t0, tall, tprol
   real(R8P), allocatable :: q_ng(:,:,:,:,:)

   call reset()
   do r=1, 3 ; call field_fnl%update_ghost_local_gpu(q_gpu=q_gpu, seam_guard=.false.) ; enddo
   do r=1, nrep
      !$acc wait
      t0 = wall() ; call field_fnl%update_ghost_local_gpu(q_gpu=q_gpu, seam_guard=.false.) ; !$acc wait
      t(r) = wall() - t0
   enddo
   tall = median(t)
   ! V0ng's prolongation rows alone, from the initial state (the reference of D3)
   call reset()
   do r=1, nrep
      !$acc wait
      t0 = wall()
      if (nrows_k(2) > 0) call update_ghost_local_gpu_dev(ngc=ngc, seam_ghost_fill=fill, seam_restriction=restr, ri=ri, &
                                  rj=rj, rk=rk, l_map_ghost_cell_gpu=rows_k2, q_gpu=q_gpu, seam_density_var=0)
      !$acc wait
      t(r) = wall() - t0
   enddo
   tprol = median(t)
   call fetch(q_ng)
   print '(A,F10.3,A,F10.3,A)', 'V0ng ADAM rows, guard off         median ', tall*1e3_R8P, ' ms  (prolongation rows ', &
         tprol*1e3_R8P, ')'
   if (fill /= SEAM_FILL_TRICUBIC .or. d3w == 0) return
   call reset()
   do r=1, 2 ; call run_d3(q_gpu) ; enddo
   do r=1, nrep
      !$acc wait
      t0 = wall() ; call run_d3(q_gpu) ; !$acc wait
      t(r) = wall() - t0
   enddo
   call fetch(q_chk)
   print '(A,F10.3,A,L1)', 'D3  prolongation per anchor (tricubic, guard off) median ', median(t)*1e3_R8P, &
         ' ms  bitwise vs V0ng prolongation ', all(q_chk == q_ng)
   endsubroutine time_v0ng_d3

   subroutine run_d3(qd)
   !< Tricubic prolongation per coarse anchor: gather the union of the children's footprints, interpolate each child.
   real(R8P), intent(inout) :: qd(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P) :: v, w, lo_, hi_, mid, m, f, loc, cnt(3), e(3), i, j, kk, ic, jc, kc, d, sub(3), p4(3), rcv(3)
   integer(I4P) :: olo(3), ohi(3), ii, jj, kq, nijk(3), rr(3), s0
   logical      :: alive(0:1,0:1,0:1)
   real(R8P)    :: fpu(-3:3,-3:3,-3:3), fp4(4,4,4)

   nijk = [n, n, n] ; rr = [ri, rj, rk] ; s0 = d1o(2)
   !$acc parallel loop independent gang vector collapse(2) DEVICEVAR(qd) present(dbr, dbs, dlo, dhi, ddl, owner, d1l, d3p) &
   !$acc& private(lo_, hi_, mid, m, f, loc, cnt, e, i, j, kk, ic, jc, kc, d, sub, p4, rcv, olo, ohi, ii, jj, kq, alive, &
   !$acc&         fpu, fp4)
   do v=1, NV
      do w=0, d3w-1
         lo_ = 1 ; hi_ = d1n(2)
         do while (lo_ < hi_)
            mid = (lo_ + hi_) / 2
            if (d3p(mid) <= w) then
               lo_ = mid + 1
            else
               hi_ = mid
            endif
         enddo
         m = lo_ ; f = d1l(s0 + m) ; loc = w - d3p(m - 1)
         cnt = dhi(:,f) - dlo(:,f) + 1
         e(1) = mod(loc, cnt(1)) ; e(2) = mod(loc / cnt(1), cnt(2)) ; e(3) = loc / (cnt(1) * cnt(2))
         i = dlo(1,f) + e(1) ; j = dlo(2,f) + e(2) ; kk = dlo(3,f) + e(3)
         ! the children this anchor writes, and the union of their footprint offsets
         olo = 99 ; ohi = -99
         do kc=0, rk-1 ; do jc=0, rj-1 ; do ic=0, ri-1
            rcv = [rr(1)*i + ddl(1,f) + ic, rr(2)*j + ddl(2,f) + jc, rr(3)*kk + ddl(3,f) + kc]
            alive(ic,jc,kc) = all(rcv >= 1-ngc) .and. all(rcv <= n+ngc)
            if (alive(ic,jc,kc)) alive(ic,jc,kc) = owner(dbr(f),rcv(1),rcv(2),rcv(3)) == f
            if (.not.alive(ic,jc,kc)) cycle
            sub = [ic, jc, kc] + 1
            do d=1, 3
               p4(d) = seam_shift_anchor_pos(anchor=merge(i, merge(j, kk, d == 2), d == 1), n_cells=nijk(d),           &
                                             p_centered=seam_tricubic_centered_pos(sub(d)), footprint_n=4)
               olo(d) = min(olo(d), 1 - p4(d)) ; ohi(d) = max(ohi(d), 4 - p4(d))
            enddo
         enddo ; enddo ; enddo
         if (olo(1) > ohi(1)) cycle
         do kq=olo(3), ohi(3) ; do jj=olo(2), ohi(2) ; do ii=olo(1), ohi(1)
            fpu(ii,jj,kq) = qd(dbs(f), i+ii, j+jj, kk+kq, v)
         enddo ; enddo ; enddo
         do kc=0, rk-1 ; do jc=0, rj-1 ; do ic=0, ri-1
            if (.not.alive(ic,jc,kc)) cycle
            sub = [ic, jc, kc] + 1
            do d=1, 3
               p4(d) = seam_shift_anchor_pos(anchor=merge(i, merge(j, kk, d == 2), d == 1), n_cells=nijk(d),           &
                                             p_centered=seam_tricubic_centered_pos(sub(d)), footprint_n=4)
            enddo
            do kq=1, 4 ; do jj=1, 4 ; do ii=1, 4
               fp4(ii,jj,kq) = fpu(ii-p4(1), jj-p4(2), kq-p4(3))
            enddo ; enddo ; enddo
            rcv = [rr(1)*i + ddl(1,f) + ic, rr(2)*j + ddl(2,f) + jc, rr(3)*kk + ddl(3,f) + kc]
            qd(dbr(f),rcv(1),rcv(2),rcv(3),v) = seam_interpolate_tricubic(footprint=fp4, sub=sub, anchor_pos=p4)
         enddo ; enddo ; enddo
      enddo
   enddo
   endsubroutine run_d3

   subroutine run_desc(variant, k, qd)
   !< One kind, over the descriptors: D1 flat (one work item per cell of a descriptor), D2 by shape class (the
   !< members innermost). The value of each ghost: the copy, ADAM's restriction (mean or high order with the meta of
   !< restriction_meta) or ADAM's interpolation (the meta of the row expansion), written only by its owner descriptor.
   integer(I4P), intent(in)    :: variant, k
   real(R8P),    intent(inout) :: qd(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P) :: v, w, lo_, hi_, mid, s, f, loc, e(3), cnt(3), i, j, kk, ic, jc, kc, anchor(3), rcv(3), sub(3), d
   integer(I4P) :: p4(3), p3(3), p6(3), code(3), nijk(3), rr(3), seg0, seg1, pre0, c, nmem, mm, meta, nw
   real(R8P)    :: val

   nw = merge(d1w(k), c2w(k), variant == 1)
   if (nw == 0) return
   nijk = [n, n, n] ; rr = [ri, rj, rk]
   if (variant == 1) then
      seg0 = d1o(k) + 1 ; seg1 = d1o(k+1) ; pre0 = k  ! d1p(m + k) is the prefix after descriptor m of kind k
   else
      seg0 = c2o(k) + 1 ; seg1 = c2o(k+1) ; pre0 = k
   endif
   !$acc parallel loop independent gang vector collapse(2) DEVICEVAR(qd)                                              &
   !$acc& present(dbr, dbs, dlo, dhi, ddl, owner, d1l, d1p, crep, cnm, cmo, cmem, cpre)                                 &
   !$acc& private(lo_, hi_, mid, s, f, loc, e, cnt, i, j, kk, ic, jc, kc, anchor, rcv, sub, d, p4, p3, p6, code, c,    &
   !$acc&         nmem, mm, meta, val)
   do v=1, NV
      do w=0, nw-1
         ! the descriptor (D1) or class (D2) holding work item w: binary search of the kind's prefix segment
         lo_ = seg0 ; hi_ = seg1
         do while (lo_ < hi_)
            mid = (lo_ + hi_) / 2
            if (variant == 1) then
               s = d1p(mid + pre0)
            else
               s = cpre(mid + pre0)
            endif
            if (s <= w) then
               lo_ = mid + 1
            else
               hi_ = mid
            endif
         enddo
         if (variant == 1) then
            f = d1l(lo_) ; loc = w - d1p(lo_ - 1 + pre0)
         else
            c = lo_ ; nmem = cnm(c) ; loc = w - cpre(c - 1 + pre0)
            mm = mod(loc, nmem) ; loc = loc / nmem ; f = cmem(cmo(c) + 1 + mm)
         endif
         cnt = dhi(:,f) - dlo(:,f) + 1
         if (k == KPROL) then
            sub = 0
            sub(1) = mod(loc, ri) ; loc = loc / ri ; sub(2) = mod(loc, rj) ; loc = loc / rj ; sub(3) = mod(loc, rk)
            loc = loc / rk
         endif
         e(1) = mod(loc, cnt(1)) ; e(2) = mod(loc / cnt(1), cnt(2)) ; e(3) = loc / (cnt(1) * cnt(2))
         i = dlo(1,f) + e(1) ; j = dlo(2,f) + e(2) ; kk = dlo(3,f) + e(3)
         select case(k)
         case(KCOPY)
            if (owner(dbr(f),i,j,kk) == f) qd(dbr(f),i,j,kk,v) = qd(dbs(f),i+ddl(1,f),j+ddl(2,f),kk+ddl(3,f),v)
         case(KREST)
            if (owner(dbr(f),i,j,kk) == f) then
               anchor = [rr(1)*i + ddl(1,f), rr(2)*j + ddl(2,f), rr(3)*kk + ddl(3,f)] ! first fine donor cell
               meta = 0
               if (restr /= SEAM_RESTRICT_MEAN) then
                  do d=1, 3
                     code(d) = seam_restriction_code(first=anchor(d), n_cells=nijk(d), ratio=rr(d))
                  enddo
                  if (all(code >= 0)) meta = seam_restriction_meta_pack(code=code, regime=restr)
               endif
               if (meta /= 0) then
                  qd(dbr(f),i,j,kk,v) = restrict_seam_high_order_gpu_dev(meta=meta, ngc=ngc, q_gpu=qd, v=v,         &
                                           b_send=dbs(f), i_send=anchor(1), j_send=anchor(2), k_send=anchor(3))
               else
                  val = 0._R8P
                  do kc=0,rk-1 ; do jc=0,rj-1 ; do ic=0,ri-1
                     val = val + qd(dbs(f),anchor(1)+ic,anchor(2)+jc,anchor(3)+kc,v)
                  enddo ; enddo ; enddo
                  qd(dbr(f),i,j,kk,v) = val / real(ri*rj*rk, R8P)
               endif
            endif
         case(KPROL)
            rcv = [rr(1)*i + ddl(1,f) + sub(1), rr(2)*j + ddl(2,f) + sub(2), rr(3)*kk + ddl(3,f) + sub(3)]
            if (all(rcv >= 1-ngc) .and. all(rcv <= n+ngc)) then
               if (owner(dbr(f),rcv(1),rcv(2),rcv(3)) == f) then
                  anchor = [i, j, kk]
                  sub = sub + 1
                  do d=1, 3
                     if (rr(d) == 1) then
                        sub(d) = SEAM_SUB_COINCIDENT ; p4(d) = 2 ; p3(d) = 2
                        p6(d) = seam_quintic_centered_pos(sub(d))
                     else
                        p4(d) = seam_shift_anchor_pos(anchor=anchor(d), n_cells=nijk(d),                          &
                                                      p_centered=seam_tricubic_centered_pos(sub(d)), footprint_n=4)
                        p3(d) = seam_shift_anchor_pos(anchor=anchor(d), n_cells=nijk(d),                          &
                                                      p_centered=seam_compatible_centered_pos(sub(d)), footprint_n=3)
                        p6(d) = seam_shift_anchor_pos(anchor=anchor(d), n_cells=nijk(d),                          &
                                                      p_centered=seam_quintic_centered_pos(sub(d)), footprint_n=6)
                     endif
                  enddo
                  meta = seam_meta_pack(sub=sub, p4=p4, p3=p3, p6=p6)
                  qd(dbr(f),rcv(1),rcv(2),rcv(3),v) = interp_seam_ghost_gpu_dev(regime=fill, meta=meta, ngc=ngc,    &
                                                         q_gpu=qd, v=v, b_send=dbs(f), i_send=i, j_send=j, k_send=kk, &
                                                         density_var=1)
               endif
            endif
         endselect
      enddo
   enddo
   endsubroutine run_desc
endprogram test_amr_ghost_bench
