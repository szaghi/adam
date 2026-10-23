!< ADAM, FLAIL: the elliptic solver library (issue #84).
module adam_flail_object
!< ADAM, FLAIL: the elliptic solver library (issue #84).
!<
!< One entry point for the elliptic solves of every ADAM application. It is built on the two solvers verified in issues
!< #81 and #87: the FFT solver `poisson_fft_object` (uniform grids, periodic or isolated, FFTW on the host; its FNL
!< extension runs the transforms on cuFFT) and the multigrid-preconditioned conjugate-gradient solver
!< `poisson_mg_object` (2:1 AMR grids, the symmetric composite operator, periodic or isolated through the James step,
!< its coarsest level solved by the FFT solver). The legacy smoothers of `adam_flail_legacy_object` retire once PRISM
!< moves here (issue #84, E3-E4).
!<
!< Phase E0 is a facade over the two solvers: `setup` chooses one for the current grid, as FLUME's self-gravity did,
!< and `solve`/`solve_packed` forward to it with the call sequence unchanged, so a client gets bitwise the results it
!< got from the solvers directly. The operator is the compact second-order Laplacian `L_h phi = factor (f - mean f)`
!< (periodic) or `L_h phi = factor f` (isolated); further operators and boundary conditions come in E1-E2.
!<
!< A backend's FFT solver is chosen by allocating the `fft` component before `setup` (the FNL backend allocates the
!< cuFFT extension); otherwise `setup` allocates the host one.
use :: adam_field_object,       only : field_object
use :: adam_grid_object,        only : grid_object
use :: adam_mpih_global,        only : mpih
use :: adam_poisson_fft_object, only : poisson_fft_object
use :: adam_poisson_mg_object,  only : poisson_mg_object
use :: mpi
use :: penf,                    only : I4P, R8P

implicit none
private
public :: flail_object
public :: FLAIL_BACKEND_AUTO
public :: FLAIL_BACKEND_FFT
public :: FLAIL_BACKEND_MULTIGRID

integer(I4P), parameter :: FLAIL_BACKEND_AUTO      = 0_I4P !< FFT on a uniform grid, multigrid otherwise.
integer(I4P), parameter :: FLAIL_BACKEND_FFT       = 1_I4P !< The FFT solver (uniform grids only).
integer(I4P), parameter :: FLAIL_BACKEND_MULTIGRID = 2_I4P !< The multigrid solver, on any 2:1 grid.

type :: flail_object
   !< FLAIL: the elliptic solver of a realm.
   class(poisson_fft_object), allocatable :: fft !< FFT solver; on a multigrid setup, the mold of its coarsest level.
   type(poisson_mg_object),   allocatable :: mg  !< Multigrid solver, allocated only when the current setup uses it.
   contains
      ! public methods
      procedure, pass(self) :: destroy      !< Destroy the solvers.
      procedure, pass(self) :: is_multigrid !< True when the current setup uses the multigrid solver.
      procedure, pass(self) :: iterations   !< FCG iterations of the last multigrid solve (0 for the FFT solver).
      procedure, pass(self) :: setup        !< Choose and set up the solver for the current grid.
      procedure, pass(self) :: solve        !< Solve for phi from one variable of q.
      procedure, pass(self) :: solve_packed !< Solve for phi from the packed interior source.
endtype flail_object

contains
   subroutine destroy(self)
   !< Destroy the solvers.
   class(flail_object), intent(inout) :: self !< FLAIL.

   if (allocated(self%fft)) then
      call self%fft%destroy
      deallocate(self%fft)
   endif
   if (allocated(self%mg)) then
      call self%mg%destroy
      deallocate(self%mg)
   endif
   endsubroutine destroy

   pure function is_multigrid(self)
   !< True when the current setup uses the multigrid solver.
   class(flail_object), intent(in) :: self         !< FLAIL.
   logical                         :: is_multigrid !< The multigrid solver is in use.

   is_multigrid = allocated(self%mg)
   endfunction is_multigrid

   pure function iterations(self)
   !< FCG iterations of the last multigrid solve (0 for the FFT solver).
   class(flail_object), intent(in) :: self       !< FLAIL.
   integer(I4P)                    :: iterations !< Iterations.

   iterations = 0_I4P
   if (allocated(self%mg)) iterations = self%mg%iterations
   endfunction iterations

   subroutine setup(self, grid, field, isolated, backend, option_name)
   !< Choose and set up the solver for the current grid: call it on the initial grid and after every regrid.
   !<
   !< `FLAIL_BACKEND_AUTO` uses the FFT solver when every leaf is at one level and the multigrid solver otherwise;
   !< `FLAIL_BACKEND_FFT` on a grid with more than one level is fatal, the message naming `option_name` (the client's
   !< input key, default `[flail].(backend)`). The multigrid solver is rebuilt at every call, its coarsest level solved
   !< by a solver of the `fft` component's type.
   class(flail_object), intent(inout)        :: self        !< FLAIL.
   type(grid_object),   intent(in)           :: grid        !< The grid.
   type(field_object),  intent(in)           :: field       !< The field (block coordinates).
   logical,             intent(in)           :: isolated    !< Free-space boundaries, else periodic.
   integer(I4P),        intent(in)           :: backend     !< FLAIL_BACKEND_AUTO, _FFT or _MULTIGRID.
   character(*),        intent(in), optional :: option_name !< The client's input key, for the messages.
   integer(I4P)                              :: lev(2)      !< Lowest and (negated) highest leaf level over the ranks.
   integer(I4P)                              :: b           !< Counter.
   integer(I4P)                              :: ierr        !< MPI error status.

   if (.not.allocated(self%fft)) allocate(poisson_fft_object :: self%fft)
   lev = [huge(1_I4P), huge(1_I4P)]
   do b=1, field%blocks_number
      lev = min(lev, [field%coordinates(4,b), -field%coordinates(4,b)])
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, lev, 2, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
   if (allocated(self%mg)) then
      call self%mg%destroy
      deallocate(self%mg)
   endif
   if (lev(1) /= -lev(2) .and. backend == FLAIL_BACKEND_FFT) then
      if (present(option_name)) then
         call mpih%error_stop(msg=': '//option_name//' = fft needs a uniform grid (the leaves are at more than one '// &
                                  'level): use auto or multigrid')
      else
         call mpih%error_stop(msg=': [flail].(backend) = fft needs a uniform grid (the leaves are at more than one '// &
                                  'level): use auto or multigrid')
      endif
   endif
   if (backend == FLAIL_BACKEND_FFT .or. (backend == FLAIL_BACKEND_AUTO .and. lev(1) == -lev(2))) then
      call self%fft%initialize(grid=grid, field=field, isolated=isolated)
   else
      allocate(self%mg)
      call self%mg%initialize(grid=grid, field=field, coarse=self%fft, isolated=isolated)
   endif
   endsubroutine setup

   subroutine solve(self, field, q, iv, factor, phi, energy)
   !< Solve for phi from the variable `iv` of q, every ghost layer of phi filled; optionally the energy
   !< `sum_cells f phi dV` (every rank).
   class(flail_object), intent(inout)         :: self       !< FLAIL.
   type(field_object),  intent(in)            :: field      !< The field (coordinates).
   real(R8P),           intent(in)            :: q(1:,1:,1:,1:,1:) !< Conservative [nv,i,j,k,nb], ghosts included.
   integer(I4P),        intent(in)            :: iv         !< Variable of q, the source.
   real(R8P),           intent(in)            :: factor     !< Source factor (4 pi G for gravity).
   real(R8P),           intent(inout)         :: phi(1:,1:,1:,1:) !< Potential [i,j,k,nb], ghosts included.
   real(R8P),           intent(out), optional :: energy     !< sum_cells f phi dV (all ranks).

   if (allocated(self%mg)) then
      call self%mg%solve(field=field, q=q, iv=iv, factor=factor, phi=phi, energy=energy)
   else
      call self%fft%solve(field=field, q=q, iv=iv, factor=factor, phi=phi, energy=energy)
      if (present(energy)) energy = energy * product(self%fft%dxyz) ! uniform grid: one cell volume
   endif
   endsubroutine solve

   subroutine solve_packed(self, field, buf, factor, phi, energy)
   !< Solve for phi from the packed interior source of the active blocks (block-major, `i` fastest), every ghost layer
   !< of phi filled; optionally the energy `sum_cells f phi dV` (every rank).
   class(flail_object), intent(inout)         :: self    !< FLAIL.
   type(field_object),  intent(in)            :: field   !< The field (coordinates).
   real(R8P),           intent(in)            :: buf(1:) !< Packed interior source.
   real(R8P),           intent(in)            :: factor  !< Source factor (4 pi G for gravity).
   real(R8P),           intent(inout)         :: phi(1:,1:,1:,1:) !< Potential [i,j,k,nb], ghosts included.
   real(R8P),           intent(out), optional :: energy  !< sum_cells f phi dV (all ranks).

   if (allocated(self%mg)) then
      call self%mg%solve_packed(field=field, buf=buf, factor=factor, phi=phi, energy=energy)
   else
      call self%fft%solve_packed(field=field, buf=buf, factor=factor, phi=phi, energy=energy)
      if (present(energy)) energy = energy * product(self%fft%dxyz) ! uniform grid: one cell volume
   endif
   endsubroutine solve_packed
endmodule adam_flail_object
