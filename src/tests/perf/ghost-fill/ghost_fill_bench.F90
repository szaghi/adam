!< Ghost-fill microbenchmark (AMR redesign, tracks 1 and 2, 2026-11-25): how far ADAM's local ghost fill is from the
!< memory floor of a GPU, and which design choice explains the gap.
!<
!< A uniform periodic brick of nbx x nby x nbz blocks of b^3 cells, ngc ghost layers, nv variables. Every ghost cell of
!< every block (faces, edges, corners) is filled from its same-level neighbour, the rows ADAM's local map holds. Variants,
!< each changing one thing:
!<   V0  ADAM's kernel (update_ghost_local_gpu_dev, adam_fnl_field_kernels.F90): one row per ghost cell, 10 integer(I8P)
!<       columns of which 9 are read, collapse(2) with the variable loop outermost, the branch on column 9, layout
!<       q(b,i,j,k,v) with the block index fastest; rows ordered as adam_maps_object builds them (region by region, i
!<       fastest);
!<   V1  V0 with integer(int32) map columns (the bytes of the indices);
!<   V2  one precomputed linear offset per row for the donor and the receiver (16 B per row instead of 72), the variable
!<       loop inside each row (the map read once per cell, not once per variable);
!<   V3  structured slabs: no per-cell map; one kernel per neighbour direction (26) over (v, k, j, i, block) of that
!<       direction's ghost region, the donor block from a per-block neighbour table, the donor cell by a constant shift;
!<       layout q(b,i,j,k,v) as V0;
!<   V4  V3 on the layout q(i,j,k,v,b), the cell index fastest;
!<   V5  V3 fused into one launch: every ghost cell of a block numbered direction by direction, its direction and cell
!<       decoded from the number with a 26-entry table (no per-cell map, one launch), layout q(b,i,j,k,v);
!<   V6  V5 on the layout q(i,j,k,v,b);
!<   ref a contiguous copy of as many values (dst(:) = src(:)): the floor this GPU reaches, measured.
!< Every variant's array is compared bit for bit with a host fill after its timing (copies are exact).
!<
!< Usage: ghost_fill_bench nbx nby nbz b ngc nv nrep
!< Output: one line per variant: median and minimum time of one fill (ms), effective bandwidth (16 B per ghost value:
!< one read, one write), fraction of the ref copy's bandwidth, bitwise check.
program ghost_fill_bench
use, intrinsic :: iso_fortran_env, only : int32, int64, real64, output_unit
use openacc
implicit none

integer(int32)              :: nbx, nby, nbz, b, g, nv, nrep, nb, n, nrows, d, nd
integer(int32)              :: dirs(3,26)                    !< The 26 neighbour directions.
integer(int32)              :: rlo(3,26), rsh(3,26), rn(3,26) !< Per direction: ghost range start, donor shift, extent.
integer(int32)              :: pre(0:26)                     !< Ghost cells of a block before direction d (V5, V6).
integer(int32)              :: ng                            !< Ghost cells per block.
integer(int32), allocatable :: nbr(:,:)                      !< nbr(dir, block): the neighbour block.
real(real64),   allocatable :: qa(:,:,:,:,:)                 !< q(b,i,j,k,v), ADAM's layout.
real(real64),   allocatable :: qb(:,:,:,:,:)                 !< q(i,j,k,v,b), cell index fastest.
real(real64),   allocatable :: q1(:)                         !< qa flattened (V2).
real(real64),   allocatable :: ref(:,:,:,:,:)                !< The host fill (qa layout), the bitwise reference.
real(real64),   allocatable :: src(:), dst(:)                !< The ref copy.
integer(int64), allocatable :: map8(:,:)                     !< V0 map, (rows, 10).
integer(int32), allocatable :: map4(:,:)                     !< V1 map, (rows, 9).
integer(int64), allocatable :: offs(:), offr(:)              !< V2 offsets (donor, receiver) of variable 1.
integer(int64)              :: sv                            !< V2 stride of one variable.
integer(int64)              :: nvalues                       !< Ghost values per fill.
real(real64),   allocatable :: t(:)                          !< The timings of one variant.
character(len=32)           :: arg
integer(int32)              :: i

call get_args()
nb = nbx * nby * nbz
n  = b + 2 * g
nvalues = int(nb, int64) * (int(n, int64)**3 - int(b, int64)**3) * nv
write(output_unit, '(A,4I6,A,I3,A,I3,A,I5,A,F10.3,A,F8.4,A)') '# blocks ', nbx, nby, nbz, nb, ' b ', b, ' ngc ', g,  &
   ' nv ', nv, ' cells ', real(nb, real64) * real(b, real64)**3 / 1e6_real64, ' M, ghost values ', &
   real(nvalues, real64) / 1e9_real64, ' G'
call describe_device()
call build_neighbours()
call build_maps()
allocate(t(nrep))

allocate(src(nvalues), dst(nvalues))
src = 1._real64 ; dst = 0._real64
!$acc enter data copyin(src, dst)
call time_variant('ref contiguous copy', 9)
!$acc exit data delete(src, dst)
deallocate(src, dst)

allocate(qa(nb, 1-g:b+g, 1-g:b+g, 1-g:b+g, nv), ref(nb, 1-g:b+g, 1-g:b+g, 1-g:b+g, nv))
call init_a(qa)
ref = qa
call host_fill(ref)
!$acc enter data copyin(qa, map8, map4, nbr, rlo, rsh, rn, pre)

call time_variant('V0 ADAM map I8P', 0)
call time_variant('V1 map I4P', 1)
call time_variant('V3 slabs, b fastest', 3)
call time_variant('V5 slabs 1 launch, b fast', 5)
!$acc exit data delete(qa, map8, map4)
deallocate(map8, map4)

allocate(q1(int(nb, int64) * int(n, int64)**3 * nv))
call init_a(qa)
q1 = reshape(qa, [size(q1, kind=int64)])
!$acc enter data copyin(q1, offs, offr)
call time_variant('V2 offsets, v inner', 2)
!$acc exit data delete(q1, offs, offr)
deallocate(q1)

allocate(qb(1-g:b+g, 1-g:b+g, 1-g:b+g, nv, nb))
call init_b(qb)
!$acc enter data copyin(qb)
call time_variant('V4 slabs, cell fastest', 4)
call time_variant('V6 slabs 1 launch, cell f', 6)
!$acc exit data delete(qb, nbr, rlo, rsh, rn, pre)
deallocate(qb, qa)

contains
   subroutine get_args()
   !< Read the command line: nbx nby nbz b ngc nv nrep.
   integer(int32) :: k, v(7)

   v = [4, 4, 8, 38, 3, 5, 50]
   do k=1, min(7, command_argument_count())
      call get_command_argument(k, arg)
      read(arg, *) v(k)
   enddo
   nbx = v(1) ; nby = v(2) ; nbz = v(3) ; b = v(4) ; g = v(5) ; nv = v(6) ; nrep = v(7)
   if (g > b) error stop 'ghost_fill_bench: ngc > b'
   endsubroutine get_args

   subroutine describe_device()
   !< Print the device the run uses.
   integer :: ndev

   ndev = acc_get_num_devices(acc_device_nvidia)
   write(output_unit, '(A,I2,A)') '# NVIDIA devices visible: ', ndev, ' (device 0 used)'
   endsubroutine describe_device

   subroutine build_neighbours()
   !< The 26 directions (dx fastest) and the periodic neighbour of every block in each.
   integer(int32) :: dx, dy, dz, ix, iy, iz, blk

   nd = 0
   do dz=-1, 1 ; do dy=-1, 1 ; do dx=-1, 1
      if (dx == 0 .and. dy == 0 .and. dz == 0) cycle
      nd = nd + 1 ; dirs(:,nd) = [dx, dy, dz]
   enddo ; enddo ; enddo
   allocate(nbr(26, nb))
   do iz=0, nbz-1 ; do iy=0, nby-1 ; do ix=0, nbx-1
      blk = 1 + ix + nbx * (iy + nby * iz)
      do d=1, 26
         nbr(d, blk) = 1 + modulo(ix + dirs(1,d), nbx) + nbx * (modulo(iy + dirs(2,d), nby) + nby * modulo(iz + dirs(3,d), nbz))
      enddo
   enddo ; enddo ; enddo
   pre(0) = 0
   do d=1, 26
      call region(dirs(1,d), rlo(1,d), rn(1,d), rsh(1,d))
      call region(dirs(2,d), rlo(2,d), rn(2,d), rsh(2,d))
      call region(dirs(3,d), rlo(3,d), rn(3,d), rsh(3,d))
      rn(:,d) = rn(:,d) - rlo(:,d) + 1 ! region() returns the upper bound: keep the extent
      pre(d) = pre(d-1) + rn(1,d) * rn(2,d) * rn(3,d)
   enddo
   ng = pre(26)
   endsubroutine build_neighbours

   pure subroutine region(dd, lo, hi, sh)
   !< The ghost range of a block along one axis for direction component dd, and the shift to the donor index.
   integer(int32), intent(in)  :: dd
   integer(int32), intent(out) :: lo, hi, sh

   select case(dd)
   case(-1) ; lo = 1 - g ; hi = 0     ; sh = b
   case( 0) ; lo = 1     ; hi = b     ; sh = 0
   case( 1) ; lo = b + 1 ; hi = b + g ; sh = -b
   endselect
   endsubroutine region

   subroutine build_maps()
   !< V0/V1 rows (region by region, k, j, i with i fastest, as adam_maps_object) and V2 offsets.
   integer(int32) :: blk, ilo, ihi, ish, jlo, jhi, jsh, klo, khi, ksh, ii, jj, kk
   integer(int64) :: r

   nrows = 0
   do d=1, 26
      call region(dirs(1,d), ilo, ihi, ish) ; call region(dirs(2,d), jlo, jhi, jsh) ; call region(dirs(3,d), klo, khi, ksh)
      nrows = nrows + (ihi - ilo + 1) * (jhi - jlo + 1) * (khi - klo + 1)
   enddo
   nrows = nrows * nb
   allocate(map8(nrows, 10), map4(nrows, 9), offs(nrows), offr(nrows))
   sv = int(nb, int64) * int(n, int64)**3
   r = 0
   do blk=1, nb
      do d=1, 26
         call region(dirs(1,d), ilo, ihi, ish) ; call region(dirs(2,d), jlo, jhi, jsh) ; call region(dirs(3,d), klo, khi, ksh)
         do kk=klo, khi ; do jj=jlo, jhi ; do ii=ilo, ihi
            r = r + 1
            map8(r,:) = [int(nbr(d,blk), int64), int(blk, int64), int(ii+ish, int64), int(jj+jsh, int64),            &
                         int(kk+ksh, int64), int(ii, int64), int(jj, int64), int(kk, int64), 1_int64, 0_int64]
            map4(r,:) = int(map8(r,1:9), int32)
            offs(r) = lin(nbr(d,blk), ii+ish, jj+jsh, kk+ksh)
            offr(r) = lin(blk, ii, jj, kk)
         enddo ; enddo ; enddo
      enddo
   enddo
   write(output_unit, '(A,I12,A,F8.3,A)') '# map rows ', nrows, ' (V0 map ', real(nrows, real64) * 80 / 1e9_real64, ' GB)'
   endsubroutine build_maps

   pure integer(int64) function lin(blk, ii, jj, kk)
   !< The 1-based linear index of (blk, ii, jj, kk, v=1) in qa(nb, 1-g:b+g, 1-g:b+g, 1-g:b+g, nv).
   integer(int32), intent(in) :: blk, ii, jj, kk

   lin = int(blk, int64) + int(nb, int64) * (int(ii+g-1, int64) + int(n, int64) * (int(jj+g-1, int64) +             &
         int(n, int64) * int(kk+g-1, int64)))
   endfunction lin

   pure real(real64) function val(blk, ii, jj, kk, v)
   !< A distinct value per interior cell and variable (exact in binary64).
   integer(int32), intent(in) :: blk, ii, jj, kk, v

   val = real(v, real64) + real(blk, real64) * 1e-3_real64 + real(ii + 64 * (jj + 64 * kk), real64) * 1e-9_real64
   endfunction val

   subroutine init_a(q)
   !< Interior values, ghosts set to -1 (qa layout).
   real(real64), intent(out) :: q(nb, 1-g:b+g, 1-g:b+g, 1-g:b+g, nv)
   integer(int32)            :: blk, ii, jj, kk, v

   q = -1._real64
   do v=1, nv ; do kk=1, b ; do jj=1, b ; do ii=1, b ; do blk=1, nb
      q(blk,ii,jj,kk,v) = val(blk, ii, jj, kk, v)
   enddo ; enddo ; enddo ; enddo ; enddo
   endsubroutine init_a

   subroutine init_b(q)
   !< Interior values, ghosts set to -1 (qb layout).
   real(real64), intent(out) :: q(1-g:b+g, 1-g:b+g, 1-g:b+g, nv, nb)
   integer(int32)            :: blk, ii, jj, kk, v

   q = -1._real64
   do blk=1, nb ; do v=1, nv ; do kk=1, b ; do jj=1, b ; do ii=1, b
      q(ii,jj,kk,v,blk) = val(blk, ii, jj, kk, v)
   enddo ; enddo ; enddo ; enddo ; enddo
   endsubroutine init_b

   subroutine host_fill(q)
   !< The reference fill on the host (qa layout).
   real(real64), intent(inout) :: q(nb, 1-g:b+g, 1-g:b+g, 1-g:b+g, nv)
   integer(int32)              :: blk, ilo, ihi, ish, jlo, jhi, jsh, klo, khi, ksh, ii, jj, kk, v

   do v=1, nv ; do blk=1, nb ; do d=1, 26
      call region(dirs(1,d), ilo, ihi, ish) ; call region(dirs(2,d), jlo, jhi, jsh) ; call region(dirs(3,d), klo, khi, ksh)
      do kk=klo, khi ; do jj=jlo, jhi ; do ii=ilo, ihi
         q(blk,ii,jj,kk,v) = q(nbr(d,blk), ii+ish, jj+jsh, kk+ksh, v)
      enddo ; enddo ; enddo
   enddo ; enddo ; enddo
   endsubroutine host_fill

   subroutine time_variant(name, kind)
   !< Warm up, time nrep fills one by one, check bitwise against the host fill, print one line.
   character(*),   intent(in) :: name
   integer(int32), intent(in) :: kind
   integer(int64)             :: c0, c1, rate
   integer(int32)             :: k
   logical                    :: ok
   real(real64)               :: tmed, tmin, bw, bwref
   real(real64), save         :: bw_ref = -1._real64

   call reset(kind) ! every ghost back to -1: the check below sees this variant's fills only
   do k=1, 3
      call fill(kind)
   enddo
   !$acc wait
   call system_clock(count_rate=rate)
   do k=1, nrep
      call system_clock(c0)
      call fill(kind)
      !$acc wait
      call system_clock(c1)
      t(k) = real(c1 - c0, real64) / real(rate, real64)
   enddo
   call sort(t)
   tmed = t((nrep + 1) / 2) ; tmin = t(1)
   bw = 16._real64 * real(nvalues, real64) / tmed / 1e9_real64
   ok = check(kind)
   if (kind == 9) bw_ref = bw
   bwref = merge(bw / bw_ref, -1._real64, bw_ref > 0._real64)
   write(output_unit, '(A28,A,F10.3,A,F10.3,A,F9.1,A,F7.3,A,L1)') name, ' median ', tmed*1e3_real64, ' ms  min ',  &
      tmin*1e3_real64, ' ms  ', bw, ' GB/s  of ref ', bwref, '  bitwise ', ok
   flush(output_unit)
   endsubroutine time_variant

   subroutine reset(kind)
   !< Ghosts to -1 on the device before a variant runs.
   integer(int32), intent(in) :: kind

   select case(kind)
   case(0, 1, 3, 5)
      call init_a(qa)
      !$acc update device(qa)
   case(2)
      call init_a(qa)
      q1 = reshape(qa, [size(q1, kind=int64)])
      !$acc update device(q1)
   case(4, 6)
      call init_b(qb)
      !$acc update device(qb)
   case(9)
      dst = 0._real64
      !$acc update device(dst)
   endselect
   endsubroutine reset

   subroutine fill(kind)
   !< One ghost fill with variant kind.
   integer(int32), intent(in) :: kind

   select case(kind)
   case(0) ; call fill_v0()
   case(1) ; call fill_v1()
   case(2) ; call fill_v2()
   case(3) ; call fill_v3()
   case(4) ; call fill_v4()
   case(5) ; call fill_v5()
   case(6) ; call fill_v6()
   case(9) ; call copy_ref()
   endselect
   endsubroutine fill

   subroutine fill_v0()
   !< ADAM's update_ghost_local_gpu_dev, same-level rows (flag 1).
   integer(int32) :: v, mf, b_send, b_recv, i_send, j_send, k_send, i_recv, j_recv, k_recv, one_or_eight

   !$acc parallel loop independent gang vector collapse(2) present(map8, qa)                                        &
   !$acc& private(b_send, b_recv, i_send, j_send, k_send, i_recv, j_recv, k_recv, one_or_eight)
   do v=1, nv
      do mf=1, nrows
         b_send       = int(map8(mf,1), int32)
         b_recv       = int(map8(mf,2), int32)
         i_send       = int(map8(mf,3), int32)
         j_send       = int(map8(mf,4), int32)
         k_send       = int(map8(mf,5), int32)
         i_recv       = int(map8(mf,6), int32)
         j_recv       = int(map8(mf,7), int32)
         k_recv       = int(map8(mf,8), int32)
         one_or_eight = int(map8(mf,9), int32)
         if (one_or_eight == 1) then
            qa(b_recv,i_recv,j_recv,k_recv,v) = qa(b_send,i_send,j_send,k_send,v)
         else
            qa(b_recv,i_recv,j_recv,k_recv,v) = 0._real64
         endif
      enddo
   enddo
   endsubroutine fill_v0

   subroutine fill_v1()
   !< V0 with int32 map columns.
   integer(int32) :: v, mf

   !$acc parallel loop independent gang vector collapse(2) present(map4, qa)
   do v=1, nv
      do mf=1, nrows
         if (map4(mf,9) == 1) then
            qa(map4(mf,2),map4(mf,6),map4(mf,7),map4(mf,8),v) = qa(map4(mf,1),map4(mf,3),map4(mf,4),map4(mf,5),v)
         else
            qa(map4(mf,2),map4(mf,6),map4(mf,7),map4(mf,8),v) = 0._real64
         endif
      enddo
   enddo
   endsubroutine fill_v1

   subroutine fill_v2()
   !< One donor and one receiver offset per row, the variables inside the row.
   integer(int32) :: v, mf
   integer(int64) :: os, orc

   !$acc parallel loop independent gang vector present(offs, offr, q1) private(os, orc)
   do mf=1, nrows
      os = offs(mf) ; orc = offr(mf)
      !$acc loop seq
      do v=0, nv-1
         q1(orc + v * sv) = q1(os + v * sv)
      enddo
   enddo
   endsubroutine fill_v2

   subroutine fill_v3()
   !< Structured slabs, layout q(b,i,j,k,v): one kernel per direction, the block index innermost.
   integer(int32) :: ilo, ihi, ish, jlo, jhi, jsh, klo, khi, ksh, v, ii, jj, kk, blk, dd

   do dd=1, 26
      call region(dirs(1,dd), ilo, ihi, ish) ; call region(dirs(2,dd), jlo, jhi, jsh) ; call region(dirs(3,dd), klo, khi, ksh)
      !$acc parallel loop independent gang vector collapse(5) present(qa, nbr) async(1)
      do v=1, nv
         do kk=klo, khi
            do jj=jlo, jhi
               do ii=ilo, ihi
                  do blk=1, nb
                     qa(blk,ii,jj,kk,v) = qa(nbr(dd,blk),ii+ish,jj+jsh,kk+ksh,v)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   !$acc wait(1)
   endsubroutine fill_v3

   subroutine fill_v4()
   !< Structured slabs, layout q(i,j,k,v,b): one kernel per direction, the cell index i innermost.
   integer(int32) :: ilo, ihi, ish, jlo, jhi, jsh, klo, khi, ksh, v, ii, jj, kk, blk, dd, bs

   do dd=1, 26
      call region(dirs(1,dd), ilo, ihi, ish) ; call region(dirs(2,dd), jlo, jhi, jsh) ; call region(dirs(3,dd), klo, khi, ksh)
      !$acc parallel loop independent gang vector collapse(5) present(qb, nbr) private(bs) async(1)
      do blk=1, nb
         do v=1, nv
            do kk=klo, khi
               do jj=jlo, jhi
                  do ii=ilo, ihi
                     bs = nbr(dd,blk)
                     qb(ii,jj,kk,v,blk) = qb(ii+ish,jj+jsh,kk+ksh,v,bs)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   !$acc wait(1)
   endsubroutine fill_v4

   subroutine fill_v5()
   !< Structured slabs in one launch, layout q(b,i,j,k,v): (v, ghost number, block), the block index innermost.
   integer(int32) :: v, e, blk, dd, loc, ii, jj, kk

   !$acc parallel loop independent gang vector collapse(3) present(qa, nbr, rlo, rsh, rn, pre) private(dd, loc, ii, jj, kk)
   do v=1, nv
      do e=0, ng-1
         do blk=1, nb
            dd = 1
            do while (pre(dd) <= e)
               dd = dd + 1
            enddo
            loc = e - pre(dd-1)
            ii = rlo(1,dd) + mod(loc, rn(1,dd))
            jj = rlo(2,dd) + mod(loc / rn(1,dd), rn(2,dd))
            kk = rlo(3,dd) + loc / (rn(1,dd) * rn(2,dd))
            qa(blk,ii,jj,kk,v) = qa(nbr(dd,blk),ii+rsh(1,dd),jj+rsh(2,dd),kk+rsh(3,dd),v)
         enddo
      enddo
   enddo
   endsubroutine fill_v5

   subroutine fill_v6()
   !< Structured slabs in one launch, layout q(i,j,k,v,b): (block, v, ghost number), the ghost number innermost.
   integer(int32) :: v, e, blk, dd, loc, ii, jj, kk

   !$acc parallel loop independent gang vector collapse(3) present(qb, nbr, rlo, rsh, rn, pre) private(dd, loc, ii, jj, kk)
   do blk=1, nb
      do v=1, nv
         do e=0, ng-1
            dd = 1
            do while (pre(dd) <= e)
               dd = dd + 1
            enddo
            loc = e - pre(dd-1)
            ii = rlo(1,dd) + mod(loc, rn(1,dd))
            jj = rlo(2,dd) + mod(loc / rn(1,dd), rn(2,dd))
            kk = rlo(3,dd) + loc / (rn(1,dd) * rn(2,dd))
            qb(ii,jj,kk,v,blk) = qb(ii+rsh(1,dd),jj+rsh(2,dd),kk+rsh(3,dd),v,nbr(dd,blk))
         enddo
      enddo
   enddo
   endsubroutine fill_v6

   subroutine copy_ref()
   !< The contiguous copy of nvalues doubles: the floor.
   integer(int64) :: k

   !$acc parallel loop independent gang vector present(src, dst)
   do k=1, nvalues
      dst(k) = src(k)
   enddo
   endsubroutine copy_ref

   logical function check(kind)
   !< Bitwise comparison of the variant's array with the host fill.
   integer(int32), intent(in) :: kind
   integer(int32)             :: blk, ii, jj, kk, v

   check = .true.
   select case(kind)
   case(0, 1, 3, 5)
      !$acc update self(qa)
      check = all(qa == ref)
   case(2)
      !$acc update self(q1)
      check = all(reshape(q1, shape(ref)) == ref)
   case(4, 6)
      !$acc update self(qb)
      do blk=1, nb ; do v=1, nv ; do kk=1-g, b+g ; do jj=1-g, b+g ; do ii=1-g, b+g
         if (qb(ii,jj,kk,v,blk) /= ref(blk,ii,jj,kk,v)) check = .false.
      enddo ; enddo ; enddo ; enddo ; enddo
   case(9)
      !$acc update self(dst)
      check = all(dst == 1._real64)
   endselect
   endfunction check

   subroutine sort(a)
   !< Insertion sort (nrep is small).
   real(real64), intent(inout) :: a(:)
   integer(int32)              :: k, l
   real(real64)                :: x

   do k=2, size(a)
      x = a(k) ; l = k - 1
      do while (l >= 1)
         if (a(l) <= x) exit
         a(l+1) = a(l) ; l = l - 1
      enddo
      a(l+1) = x
   enddo
   endsubroutine sort
endprogram ghost_fill_bench
