!< ADAM, FNL multigrid-preconditioned conjugate-gradient Poisson solver on AMR grids (OpenACC).
module adam_fnl_poisson_mg_object
!< ADAM, FNL multigrid-preconditioned conjugate-gradient Poisson solver on AMR grids (OpenACC).
!<
!< Extends `poisson_mg_object` (issue #84, E6a): the leaves, the hierarchy and every operator are assembled on the host
!< as there; `work_place` then puts the operators (the composite operator, the ghost fill, the level operators,
!< ghost terms and transfers, CSR), the level diagonals and colours, the cell volumes, the maps of the leaves into the
!< level spaces and the work vectors on the device, where they stay until `destroy`. The FCG iteration and the V-cycle
!< run there through the hooks:
!<   - `apply`: the halo values packed on the device, exchanged by `MPI_ALLTOALLV` on the device buffers (GPU-aware
!<     MPI, as the ghost exchange; none on one rank), then one thread per row, its entries summed in the host's order;
!<   - `smooth`: the red-black Gauss-Seidel sweeps, one kernel per colour after its halo exchange. A row couples only
!<     to cells of the other colour (the colour is the global cell parity, the level operators are 7-point: measured,
!<     no same-colour entry), so the parallel update of a colour is the host's sequential one;
!<   - the vector hooks (fill, copy, combinations, dots and the max-norm by device reductions, the permutations of the
!<     leaves into the levels).
!< The root goes through the host to the FFT solver (cuFFT by default: `initialize` passes `poisson_fft_fnl_object`
!< when the caller gives no `coarse`), and so do the right-hand side, the James step and the ghost fill's placement.
!< The results equal the host's to round-off: the dots sum in another order and the device contracts products into
!< FMAs.
use :: adam_field_object, only : field_object
use :: adam_fnl_poisson_fft_object, only : poisson_fft_fnl_object
use :: adam_grid_object, only : grid_object
use :: adam_mpih_global, only : mpih
use :: adam_poisson_fft_object, only : poisson_fft_object
use :: adam_poisson_mg_object, only : poisson_mg_object, dcsr
use :: mpi
use :: penf, only : I4P, R8P

implicit none
private
public :: poisson_mg_fnl_object

type, extends(poisson_mg_object) :: poisson_mg_fnl_object
   !< FNL multigrid-preconditioned CG Poisson solver on AMR grids (OpenACC).
   logical                :: on_device=.false. !< The operators and the work vectors are on the device.
   real(R8P), allocatable :: dsb(:)            !< Device: the halo send buffer (the largest operator's).
   real(R8P), allocatable :: dhb(:)            !< Device: the halo receive buffer (the largest operator's).
   contains
      ! public methods
      procedure, pass(self) :: initialize => initialize_fnl !< As the host's, the coarsest level on cuFFT by default.
      ! backend hooks
      procedure, pass(self) :: apply       => apply_fnl       !< y = M x on the device.
      procedure, pass(self) :: smooth      => smooth_fnl      !< Red-black Gauss-Seidel on the device.
      procedure, nopass     :: vec_fill    => vec_fill_fnl    !< x(1:n) = a on the device.
      procedure, nopass     :: vec_copy    => vec_copy_fnl    !< y(1:n) = x(1:n) on the device.
      procedure, nopass     :: vec_lincomb => vec_lincomb_fnl !< w = x + a y on the device.
      procedure, nopass     :: vec_axpy    => vec_axpy_fnl    !< y = y + a x on the device.
      procedure, nopass     :: vec_dot     => vec_dot_fnl     !< sum x y on the device and the ranks.
      procedure, nopass     :: vec_relmax  => vec_relmax_fnl  !< max |r / vol| on the device and the ranks.
      procedure, nopass     :: vec_permute => vec_permute_fnl !< y(iy) = x(ix) on the device.
      procedure, nopass     :: vec_update  => vec_update_fnl  !< A work vector to or from the device.
      procedure, pass(self) :: work_place  => work_place_fnl  !< The operators and work vectors to or from the device.
endtype poisson_mg_fnl_object

contains
   ! public methods
   subroutine initialize_fnl(self, grid, field, tolerance, max_iterations, coarse, isolated, bc, root_cells)
   !< As `poisson_mg_object%initialize`; without `coarse`, the coarsest level (and the James step's Green's function) on
   !< cuFFT.
   class(poisson_mg_fnl_object), intent(inout)        :: self           !< The solver.
   type(grid_object),            intent(in)           :: grid           !< The grid.
   type(field_object),           intent(in)           :: field          !< The field (block coordinates).
   real(R8P),                    intent(in), optional :: tolerance      !< Relative residual to stop.
   integer(I4P),                 intent(in), optional :: max_iterations !< FCG iterations at most.
   class(poisson_fft_object),    intent(in), optional :: coarse         !< Mold of the coarsest-level solver.
   logical,                      intent(in), optional :: isolated       !< Free-space boundaries (James).
   integer(I4P),                 intent(in), optional :: bc(2,3)        !< Kinds per face [lo/hi, axis].
   integer(I4P),                 intent(in), optional :: root_cells     !< Most cells of the root.
   type(poisson_fft_fnl_object)                       :: cufft_mold     !< The cuFFT root, by default.

   ! through the parent component the passed object is a poisson_mg_object: its destroy and work_place are the host's
   ! (no-ops), so the device copies of a previous setup are released here first and the new ones placed after
   call self%destroy
   if (present(coarse)) then
      call self%poisson_mg_object%initialize(grid=grid, field=field, tolerance=tolerance,                        &
                                             max_iterations=max_iterations, coarse=coarse, isolated=isolated, &
                                             bc=bc, root_cells=root_cells)
   else
      call self%poisson_mg_object%initialize(grid=grid, field=field, tolerance=tolerance,                            &
                                             max_iterations=max_iterations, coarse=cufft_mold, isolated=isolated, &
                                             bc=bc, root_cells=root_cells)
   endif
   call self%work_place(attach=.true.)
   endsubroutine initialize_fnl

   ! backend hooks
   subroutine apply_fnl(self, m, x, y)
   !< y = M x on the device: the halo of x packed, exchanged (GPU-aware MPI_ALLTOALLV), then the local rows.
   class(poisson_mg_fnl_object), intent(inout) :: self  !< The solver.
   type(dcsr),                   intent(in)    :: m     !< Operator (on the device).
   real(R8P),                    intent(in)    :: x(1:) !< Local column-space vector (on the device).
   real(R8P),                    intent(inout) :: y(1:) !< Result (rows, on the device).

   call halo_dev(self, m, x)
   if (m%nrow > 0) call spmv_dev(nrow=m%nrow, nloc=m%nloc, nnz=size(m%col), nx=size(x), ny=size(y), nh=size(self%dhb), &
                                 rowptr=m%rowptr, col=m%col, val=m%val, x=x, h=self%dhb, y=y)
   endsubroutine apply_fnl

   subroutine smooth_fnl(self, lv, e, rhs, order)
   !< Red-black Gauss-Seidel sweeps on level lv on the device: per colour, the halo of e, then every row of the colour
   !< (they read only the other colour's cells).
   class(poisson_mg_fnl_object), intent(inout) :: self    !< The solver.
   integer(I4P),                 intent(in)    :: lv      !< Level.
   real(R8P),                    intent(inout) :: e(1:)   !< Correction (on the device).
   real(R8P),                    intent(in)    :: rhs(1:) !< Right-hand side (on the device).
   integer(I4P),                 intent(in)    :: order   !< 0: colour 0 first, 1: colour 1 first.
   integer(I4P)                                :: sw, cc_ !< Counters.

   do sw=1, self%sweeps
      do cc_=0, 1
         call halo_dev(self, self%al(lv), e)
         if (self%al(lv)%nrow > 0) &
            call sweep_dev(nrow=self%al(lv)%nrow, nloc=self%al(lv)%nloc, nnz=size(self%al(lv)%col), ne=size(e),  &
                           nr=size(rhs), nh=size(self%dhb), nd=size(self%diag(lv)%v), colr=mod(cc_ + order, 2), &
                           rowptr=self%al(lv)%rowptr, col=self%al(lv)%col, val=self%al(lv)%val,                &
                           diag=self%diag(lv)%v, colour=self%colour(lv)%v, rhs=rhs, h=self%dhb, e=e)
      enddo
   enddo
   endsubroutine smooth_fnl

   subroutine vec_fill_fnl(x, a, n)
   !< x(1:n) = a on the device.
   real(R8P),    intent(inout) :: x(1:) !< Vector (on the device).
   real(R8P),    intent(in)    :: a     !< Value.
   integer(I4P), intent(in)    :: n     !< Entries.

   if (n > 0) call fill_dev(size(x), x, a, n)
   endsubroutine vec_fill_fnl

   subroutine vec_copy_fnl(x, y, n)
   !< y(1:n) = x(1:n) on the device.
   real(R8P),    intent(in)    :: x(1:) !< Source (on the device).
   real(R8P),    intent(inout) :: y(1:) !< Destination (on the device).
   integer(I4P), intent(in)    :: n     !< Entries.

   if (n > 0) call lincomb_dev(size(x), x, 0._R8P, size(x), x, size(y), y, n, .false.)
   endsubroutine vec_copy_fnl

   subroutine vec_lincomb_fnl(x, a, y, w, n)
   !< w(1:n) = x(1:n) + a y(1:n) on the device.
   real(R8P),    intent(in)    :: x(1:) !< First vector (on the device).
   real(R8P),    intent(in)    :: a     !< Coefficient.
   real(R8P),    intent(in)    :: y(1:) !< Second vector (on the device).
   real(R8P),    intent(inout) :: w(1:) !< Result (on the device).
   integer(I4P), intent(in)    :: n     !< Entries.

   if (n > 0) call lincomb_dev(size(x), x, a, size(y), y, size(w), w, n, .true.)
   endsubroutine vec_lincomb_fnl

   subroutine vec_axpy_fnl(a, x, y, n)
   !< y(1:n) = y(1:n) + a x(1:n) on the device.
   real(R8P),    intent(in)    :: a     !< Coefficient.
   real(R8P),    intent(in)    :: x(1:) !< Vector (on the device).
   real(R8P),    intent(inout) :: y(1:) !< Updated vector (on the device).
   integer(I4P), intent(in)    :: n     !< Entries.

   if (n > 0) call axpy_dev(a, size(x), x, size(y), y, n)
   endsubroutine vec_axpy_fnl

   function vec_dot_fnl(x, y, n) result(d)
   !< sum_i x_i y_i over the local n entries (a device reduction) and the ranks.
   real(R8P),    intent(in) :: x(1:) !< First vector (on the device).
   real(R8P),    intent(in) :: y(1:) !< Second vector (on the device).
   integer(I4P), intent(in) :: n     !< Local entries.
   real(R8P)                :: d     !< Result.
   integer(I4P)             :: ierr  !< MPI error status.

   d = 0._R8P
   if (n > 0) d = dot_dev(size(x), x, size(y), y, n)
   call MPI_ALLREDUCE(MPI_IN_PLACE, d, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
   endfunction vec_dot_fnl

   function vec_relmax_fnl(r, vol, n, fmax) result(rm)
   !< max |r / vol| over the local n entries (a device reduction) and the ranks, over fmax (when positive).
   real(R8P),    intent(in) :: r(1:)   !< Residual (on the device).
   real(R8P),    intent(in) :: vol(1:) !< Cell volumes (on the device).
   integer(I4P), intent(in) :: n       !< Local entries.
   real(R8P),    intent(in) :: fmax    !< Normalisation.
   real(R8P)                :: rm      !< Result.
   integer(I4P)             :: ierr    !< MPI error status.

   rm = 0._R8P
   if (n > 0) rm = relmax_dev(size(r), r, size(vol), vol, n)
   call MPI_ALLREDUCE(MPI_IN_PLACE, rm, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   if (fmax > 0._R8P) rm = rm / fmax
   endfunction vec_relmax_fnl

   subroutine vec_permute_fnl(x, ix, y, iy, n)
   !< y(iy(i)) = x(ix(i)), i = 1..n, on the device.
   real(R8P),    intent(in)    :: x(1:)  !< Source (on the device).
   integer(I4P), intent(in)    :: ix(1:) !< Source entries (on the device).
   real(R8P),    intent(inout) :: y(1:)  !< Destination (on the device).
   integer(I4P), intent(in)    :: iy(1:) !< Destination entries (on the device).
   integer(I4P), intent(in)    :: n      !< Entries.

   if (n > 0) call permute_dev(size(x), x, size(ix), ix, size(y), y, size(iy), iy, n)
   endsubroutine vec_permute_fnl

   subroutine vec_update_fnl(x, to_device)
   !< A work vector to (`to_device`) or from the device.
   real(R8P), intent(inout) :: x(1:)     !< Work vector (on the device).
   logical,   intent(in)    :: to_device !< Direction.

   call update_dev(size(x), x, to_device)
   endsubroutine vec_update_fnl

   subroutine work_place_fnl(self, attach)
   !< Put (`attach`) the operators, the level diagonals, colours and maps, the volumes and the work vectors on the
   !< device, with the halo buffers sized for the largest operator; or release them.
   class(poisson_mg_fnl_object), intent(inout) :: self   !< The solver.
   logical,                      intent(in)    :: attach !< Put, or release.
   integer(I4P)                                :: lv, ns, nh !< Level, buffer sizes.

   if (.not.self%is_initialized) call mpih%error_stop(msg='poisson_mg_fnl_object%work_place: not initialized')
   if (attach .eqv. self%on_device) return
   if (attach) then
      ns = max(size(self%a%sidx), size(self%g%sidx)) ; nh = max(self%a%nhalo, self%g%nhalo)
      do lv=self%lbot+1, self%lmax
         ns = max(ns, size(self%al(lv)%sidx), size(self%el(lv)%sidx), size(self%rl(lv)%sidx), size(self%pl(lv)%sidx))
         nh = max(nh, self%al(lv)%nhalo, self%el(lv)%nhalo, self%rl(lv)%nhalo, self%pl(lv)%nhalo)
      enddo
      allocate(self%dsb(max(ns,1)), self%dhb(max(nh,1)))
      self%dsb = 0._R8P ; self%dhb = 0._R8P
      call place_r(self%dsb, .true.) ; call place_r(self%dhb, .true.)
   endif
   call place_csr(self%a, attach) ; call place_csr(self%g, attach)
   do lv=self%lbot+1, self%lmax
      call place_csr(self%al(lv), attach) ; call place_csr(self%el(lv), attach)
      call place_csr(self%rl(lv), attach) ; call place_csr(self%pl(lv), attach)
      call place_r(self%diag(lv)%v, attach) ; call place_r(self%colour(lv)%v, attach)
   enddo
   do lv=self%lbot, self%lmax
      call place_r(self%lr(lv)%v, attach) ; call place_r(self%le(lv)%v, attach)
      call place_r(self%lt1(lv)%v, attach) ; call place_r(self%lt2(lv)%v, attach)
      call place_i(self%mlev(lv)%v, attach) ; call place_i(self%mleaf(lv)%v, attach)
   enddo
   call place_r(self%vol, attach) ; call place_r(self%u, attach) ; call place_r(self%u0, attach)
   call place_r(self%wf, attach) ; call place_r(self%wr, attach) ; call place_r(self%wz, attach)
   call place_r(self%wp, attach) ; call place_r(self%wq, attach) ; call place_r(self%wpp, attach)
   call place_r(self%wqp, attach) ; call place_r(self%w1, attach) ; call place_r(self%wg, attach)
   if (.not.attach) then
      call place_r(self%dsb, .false.) ; call place_r(self%dhb, .false.)
      deallocate(self%dsb, self%dhb)
   endif
   self%on_device = attach
   endsubroutine work_place_fnl

   ! private procedures
   subroutine halo_dev(self, m, x)
   !< The halo of x for operator m into self%dhb: x(sidx) packed into self%dsb, exchanged on the device buffers.
   class(poisson_mg_fnl_object), intent(inout) :: self  !< The solver.
   type(dcsr),                   intent(in)    :: m     !< Operator.
   real(R8P),                    intent(in)    :: x(1:) !< Local column-space vector (on the device).
   integer(I4P)                                :: ierr  !< MPI error status.

   if (size(m%sidx) > 0) call pack_dev(size(m%sidx), m%sidx, size(x), x, size(self%dsb), self%dsb)
   if (self%nr > 1) call alltoallv_dev(size(self%dsb), self%dsb, m%scount, m%sdispl, size(self%dhb), self%dhb, &
                                       m%rcount, m%rdispl, ierr)
   endsubroutine halo_dev

   subroutine place_csr(m, attach)
   !< A CSR operator's arrays to (`attach`) or from the device (the exchange counts stay on the host).
   type(dcsr), intent(inout) :: m      !< Operator.
   logical,    intent(in)    :: attach !< Put, or release.

   call place_i(m%rowptr, attach) ; call place_i(m%col, attach) ; call place_r(m%val, attach)
   call place_i(m%sidx, attach)
   endsubroutine place_csr

   subroutine place_r(x, attach)
   !< A real array to (`attach`, copied in) or from the device (deleted); zero-size arrays stay on the host.
   real(R8P), intent(inout) :: x(:)   !< Array.
   logical,   intent(in)    :: attach !< Put, or release.

   if (size(x) == 0) return
   if (attach) then
      !$acc enter data copyin(x)
   else
      !$acc exit data delete(x)
   endif
   endsubroutine place_r

   subroutine place_i(x, attach)
   !< An integer array to (`attach`, copied in) or from the device (deleted); zero-size arrays stay on the host.
   integer(I4P), intent(inout) :: x(:)   !< Array.
   logical,      intent(in)    :: attach !< Put, or release.

   if (size(x) == 0) return
   if (attach) then
      !$acc enter data copyin(x)
   else
      !$acc exit data delete(x)
   endif
   endsubroutine place_i

   subroutine pack_dev(ns, sidx, nx, x, nb, sb)
   !< sb(i) = x(sidx(i)), i = 1..ns, on the device.
   integer(I4P), intent(in)    :: ns, nx, nb !< Sizes.
   integer(I4P), intent(in)    :: sidx(ns)   !< Local entries to send (present).
   real(R8P),    intent(in)    :: x(nx)      !< Vector (present).
   real(R8P),    intent(inout) :: sb(nb)     !< Send buffer (present).
   integer(I4P)                :: i          !< Counter.

   !$acc parallel loop present(sidx, x, sb)
   do i=1, ns
      sb(i) = x(sidx(i))
   enddo
   endsubroutine pack_dev

   subroutine alltoallv_dev(ns, sb, sc, sd, nh, hb, rc, rd, ierr)
   !< MPI_ALLTOALLV of the device buffers (GPU-aware MPI).
   integer(I4P), intent(in)    :: ns, nh                         !< Buffer sizes.
   real(R8P),    intent(in)    :: sb(ns)                         !< Send buffer (present).
   integer(I4P), intent(in)    :: sc(0:), sd(0:), rc(0:), rd(0:) !< Counts and displacements (values).
   real(R8P),    intent(inout) :: hb(nh)                         !< Receive buffer (present).
   integer(I4P), intent(out)   :: ierr                           !< MPI error status.

   !$acc host_data use_device(sb, hb)
   call MPI_ALLTOALLV(sb, sc, sd, MPI_REAL8, hb, rc, rd, MPI_REAL8, MPI_COMM_WORLD, ierr)
   !$acc end host_data
   endsubroutine alltoallv_dev

   subroutine spmv_dev(nrow, nloc, nnz, nx, ny, nh, rowptr, col, val, x, h, y)
   !< y = M x on the device, one row per thread, the entries in the host's order.
   integer(I4P), intent(in)    :: nrow, nloc, nnz, nx, ny, nh !< Sizes.
   integer(I4P), intent(in)    :: rowptr(nrow+1)              !< Row pointers (present).
   integer(I4P), intent(in)    :: col(nnz)                    !< Column slots (present).
   real(R8P),    intent(in)    :: val(nnz)                    !< Values (present).
   real(R8P),    intent(in)    :: x(nx)                       !< Local columns (present).
   real(R8P),    intent(in)    :: h(nh)                       !< Halo columns (present).
   real(R8P),    intent(inout) :: y(ny)                       !< Result (present).
   real(R8P)                   :: acc                         !< Accumulator.
   integer(I4P)                :: i, k                        !< Counters.

   !$acc parallel loop gang vector present(rowptr, col, val, x, h, y) private(acc)
   do i=1, nrow
      acc = 0._R8P
      !$acc loop seq
      do k=rowptr(i), rowptr(i+1)-1
         if (col(k) <= nloc) then
            acc = acc + val(k) * x(col(k))
         else
            acc = acc + val(k) * h(col(k) - nloc)
         endif
      enddo
      y(i) = acc
   enddo
   endsubroutine spmv_dev

   subroutine sweep_dev(nrow, nloc, nnz, ne, nr, nh, nd, colr, rowptr, col, val, diag, colour, rhs, h, e)
   !< One colour of a red-black Gauss-Seidel sweep on the device: e(i) = (rhs(i) - sum_{j /= i} a_ij e_j) / a_ii for the
   !< rows of colour colr, the off-diagonal entries in the host's order.
   integer(I4P), intent(in)    :: nrow, nloc, nnz, ne, nr, nh, nd !< Sizes.
   integer(I4P), intent(in)    :: colr                            !< Colour of this pass.
   integer(I4P), intent(in)    :: rowptr(nrow+1)                  !< Row pointers (present).
   integer(I4P), intent(in)    :: col(nnz)                        !< Column slots (present).
   real(R8P),    intent(in)    :: val(nnz)                        !< Values (present).
   real(R8P),    intent(in)    :: diag(nd)                        !< Diagonal (present).
   real(R8P),    intent(in)    :: colour(nd)                      !< Colour of each row, 0 or 1 (present).
   real(R8P),    intent(in)    :: rhs(nr)                         !< Right-hand side (present).
   real(R8P),    intent(in)    :: h(nh)                           !< Halo of e (present).
   real(R8P),    intent(inout) :: e(ne)                           !< Correction (present).
   real(R8P)                   :: acc                             !< Off-diagonal sum.
   integer(I4P)                :: i, k                            !< Counters.

   !$acc parallel loop gang vector present(rowptr, col, val, diag, colour, rhs, h, e) private(acc)
   do i=1, nrow
      if (nint(colour(i)) == colr) then
         acc = 0._R8P
         !$acc loop seq
         do k=rowptr(i), rowptr(i+1)-1
            if (col(k) == i) cycle
            if (col(k) <= nloc) then
               acc = acc + val(k) * e(col(k))
            else
               acc = acc + val(k) * h(col(k) - nloc)
            endif
         enddo
         e(i) = (rhs(i) - acc) / diag(i)
      endif
   enddo
   endsubroutine sweep_dev

   subroutine fill_dev(nx, x, a, n)
   !< x(1:n) = a on the device.
   integer(I4P), intent(in)    :: nx, n !< Sizes.
   real(R8P),    intent(inout) :: x(nx) !< Vector (present).
   real(R8P),    intent(in)    :: a     !< Value.
   integer(I4P)                :: i     !< Counter.

   !$acc parallel loop present(x)
   do i=1, n
      x(i) = a
   enddo
   endsubroutine fill_dev

   subroutine lincomb_dev(nx, x, a, ny, y, nw, w, n, comb)
   !< w(1:n) = x(1:n) + a y(1:n) (comb) or w(1:n) = x(1:n) on the device.
   integer(I4P), intent(in)    :: nx, ny, nw, n !< Sizes.
   real(R8P),    intent(in)    :: x(nx)         !< First vector (present).
   real(R8P),    intent(in)    :: a             !< Coefficient.
   real(R8P),    intent(in)    :: y(ny)         !< Second vector (present).
   real(R8P),    intent(inout) :: w(nw)         !< Result (present).
   logical,      intent(in)    :: comb          !< Combine, or copy x.
   integer(I4P)                :: i             !< Counter.

   if (comb) then
      !$acc parallel loop present(x, y, w)
      do i=1, n
         w(i) = x(i) + a * y(i)
      enddo
   else
      !$acc parallel loop present(x, w)
      do i=1, n
         w(i) = x(i)
      enddo
   endif
   endsubroutine lincomb_dev

   subroutine axpy_dev(a, nx, x, ny, y, n)
   !< y(1:n) = y(1:n) + a x(1:n) on the device.
   real(R8P),    intent(in)    :: a         !< Coefficient.
   integer(I4P), intent(in)    :: nx, ny, n !< Sizes.
   real(R8P),    intent(in)    :: x(nx)     !< Vector (present).
   real(R8P),    intent(inout) :: y(ny)     !< Updated vector (present).
   integer(I4P)                :: i         !< Counter.

   !$acc parallel loop present(x, y)
   do i=1, n
      y(i) = y(i) + a * x(i)
   enddo
   endsubroutine axpy_dev

   function dot_dev(nx, x, ny, y, n) result(d)
   !< sum_{i<=n} x_i y_i on the device.
   integer(I4P), intent(in) :: nx, ny, n !< Sizes.
   real(R8P),    intent(in) :: x(nx)     !< First vector (present).
   real(R8P),    intent(in) :: y(ny)     !< Second vector (present).
   real(R8P)                :: d         !< Result.
   integer(I4P)             :: i         !< Counter.

   d = 0._R8P
   !$acc parallel loop present(x, y) reduction(+:d)
   do i=1, n
      d = d + x(i) * y(i)
   enddo
   endfunction dot_dev

   function relmax_dev(nr, r, nv, vol, n) result(rm)
   !< max_{i<=n} |r_i / vol_i| on the device.
   integer(I4P), intent(in) :: nr, nv, n !< Sizes.
   real(R8P),    intent(in) :: r(nr)     !< Residual (present).
   real(R8P),    intent(in) :: vol(nv)   !< Cell volumes (present).
   real(R8P)                :: rm        !< Result.
   integer(I4P)             :: i         !< Counter.

   rm = 0._R8P
   !$acc parallel loop present(r, vol) reduction(max:rm)
   do i=1, n
      rm = max(rm, abs(r(i) / vol(i)))
   enddo
   endfunction relmax_dev

   subroutine permute_dev(nx, x, nix, ix, ny, y, niy, iy, n)
   !< y(iy(i)) = x(ix(i)), i = 1..n, on the device.
   integer(I4P), intent(in)    :: nx, nix, ny, niy, n !< Sizes.
   real(R8P),    intent(in)    :: x(nx)               !< Source (present).
   integer(I4P), intent(in)    :: ix(nix)             !< Source entries (present).
   real(R8P),    intent(inout) :: y(ny)               !< Destination (present).
   integer(I4P), intent(in)    :: iy(niy)             !< Destination entries (present).
   integer(I4P)                :: i                   !< Counter.

   !$acc parallel loop present(x, ix, y, iy)
   do i=1, n
      y(iy(i)) = x(ix(i))
   enddo
   endsubroutine permute_dev

   subroutine update_dev(nx, x, to_device)
   !< x to (`to_device`) or from the device.
   integer(I4P), intent(in)    :: nx        !< Size.
   real(R8P),    intent(inout) :: x(nx)     !< Vector (present).
   logical,      intent(in)    :: to_device !< Direction.

   if (to_device) then
      !$acc update device(x)
   else
      !$acc update self(x)
   endif
   endsubroutine update_dev
endmodule adam_fnl_poisson_mg_object
