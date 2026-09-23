!< Unit test of library periodicity: ghost cells across a periodic domain boundary (FLUME #35, P0-B).
program test_periodic_ghost
!< Unit test of library periodicity: ghost cells across a periodic domain boundary (FLUME #35, P0-B).
!<
!< **Why this test exists.** PRISM's `periodic` BC is an app constant (`BC_PERIOD = 5`) that never reaches the
!< library: `grid%set_bc_type` flags a direction periodic only for `BC_PERIODIC = -1`, and the tree wraps a
!< neighbor query across the domain only for a flagged direction (`adam_tree_object.F90`). PRISM therefore fakes
!< periodicity by a same-block wrap, correct only when one block spans the periodic direction. FLUME needs true
!< library periodicity; this test establishes whether it works through the real initialization path, with many
!< blocks and many ranks, and with `set_bc_type` called AFTER `realm_object%initialize` (the PRISM order).
!<
!< **What it pins.** A realm is initialized from an INI file, all six faces are set to `BC_PERIODIC`, the grid is
!< uniformly refined into many blocks and the maps are rebuilt. The interior is filled with the smooth periodic
!< function
!<```
!< f(x,y,z) = sin(2 pi x/Lx) + 2 sin(2 pi y/Ly) + 3 sin(2 pi z/Lz)
!<```
!< and the ghosts with a sentinel; after `update_ghost_local` + `update_ghost_mpi` every face-ghost cell must hold
!< `f` evaluated at its own (out-of-domain) center, to round-off. Edge and corner ghosts are reported, not
!< asserted: directional stencils never read them.
!<
!< **Negative control.** With `--no-periodic` the faces are left non-periodic: face ghosts on the domain boundary
!< must then stay unfilled (sentinel), proving the check can fail.
!<
!< Usage: `mpirun -np N exe/test_periodic_ghost [input.ini] [--no-periodic]`.

use :: adam_common_library, only : realm_object, BC_PERIODIC
use :: adam_mpih_global,    only : mpih
use :: mpi
use :: penf,                only : I4P, R8P, str

implicit none

real(R8P),    parameter   :: PI=acos(-1._R8P)     !< Pi greek.
real(R8P),    parameter   :: SENTINEL=-1.e30_R8P  !< Value marking a ghost cell not filled by the exchange.
real(R8P),    parameter   :: TOLERANCE=1.e-12_R8P !< Round-off tolerance on face ghosts.
type(realm_object)        :: realm                !< Realm under test.
real(R8P),    allocatable :: q(:,:,:,:,:)         !< Scalar field (1, i, j, k, b).
character(999)            :: arg                  !< Command line argument.
character(:), allocatable :: filename             !< Input file name.
logical                   :: is_periodic          !< Periodic leg (default) or negative control.
real(R8P)                 :: err_max(0:3)         !< Max |q - f| per class: 0 interior, 1 face, 2 edge, 3 corner.
integer(I4P)              :: unfilled(0:3)        !< Sentinel cells left per class.
integer(I4P)              :: checked(0:3)         !< Cells checked per class.
integer(I4P)              :: ngc, ni, nj, nk      !< Grid dimensions.
integer(I4P)              :: b, i, j, k, a        !< Counters.
integer(I4P)              :: cls                  !< Cell class (number of directions outside the interior).
integer(I4P)              :: ierr                 !< MPI error status.
real(R8P)                 :: f                    !< Exact value.
logical                   :: test_passed          !< Aggregate pass flag.

filename = 'test_periodic_ghost.ini'
is_periodic = .true.
do a=1, command_argument_count()
   call get_command_argument(a, arg)
   if (trim(arg) == '--no-periodic') then
      is_periodic = .false.
   else
      filename = trim(arg)
   endif
enddo

call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)
call realm%initialize(filename=filename, memory_avail=1._R8P, nv=1_I4P)
if (is_periodic) call realm%adam%grid%set_bc_type(bc_type=[(BC_PERIODIC, a=1, 6)])

ngc = realm%adam%grid%ngc ; ni = realm%adam%grid%ni ; nj = realm%adam%grid%nj ; nk = realm%adam%grid%nk
allocate(q(1:1,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:realm%adam%field%nb))
q = 0._R8P
call realm%adam%refine_uniform(refinement_levels=realm%adam%tree%iu_ref_levels, q=q, &
                               do_mpi_redistribute=.true., do_blocks_reorder=.false.)
call realm%adam%make_comm_local_maps_ghost_bc()

q = SENTINEL
do b=1, realm%adam%field%blocks_number
   do k=1, nk
      do j=1, nj
         do i=1, ni
            q(1,i,j,k,b) = exact(b=b, i=i, j=j, k=k)
         enddo
      enddo
   enddo
enddo
call realm%adam%field%update_ghost_local(grid=realm%adam%grid, maps=realm%adam%maps, q=q)
call realm%adam%field%update_ghost_mpi(grid=realm%adam%grid, maps=realm%adam%maps, q=q)

err_max = 0._R8P ; unfilled = 0_I4P ; checked = 0_I4P
do b=1, realm%adam%field%blocks_number
   do k=1-ngc, nk+ngc
      do j=1-ngc, nj+ngc
         do i=1-ngc, ni+ngc
            cls = count([i < 1 .or. i > ni, j < 1 .or. j > nj, k < 1 .or. k > nk])
            checked(cls) = checked(cls) + 1_I4P
            if (q(1,i,j,k,b) == SENTINEL) then
               unfilled(cls) = unfilled(cls) + 1_I4P
            else
               f = exact(b=b, i=i, j=j, k=k)
               err_max(cls) = max(err_max(cls), abs(q(1,i,j,k,b) - f))
            endif
         enddo
      enddo
   enddo
enddo
call MPI_ALLREDUCE(MPI_IN_PLACE, err_max,  4, MPI_REAL8,   MPI_MAX, MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE(MPI_IN_PLACE, unfilled, 4, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE(MPI_IN_PLACE, checked,  4, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

if (mpih%myrank == 0) then
   print '(A)', 'periodic leg: '//trim(str(is_periodic))//', ranks: '//trim(str(mpih%procs_number))// &
                ', blocks: '//trim(str(realm%adam%tree%nodes_number))
   print '(A)', '  class     checked    unfilled    max|q-f|'
   print '(A,3(1X,A))', '  interior', str(checked(0)), str(unfilled(0)), str(err_max(0))
   print '(A,3(1X,A))', '  face    ', str(checked(1)), str(unfilled(1)), str(err_max(1))
   print '(A,3(1X,A))', '  edge    ', str(checked(2)), str(unfilled(2)), str(err_max(2))
   print '(A,3(1X,A))', '  corner  ', str(checked(3)), str(unfilled(3)), str(err_max(3))
endif

if (is_periodic) then
   test_passed = unfilled(1) == 0_I4P .and. err_max(1) < TOLERANCE .and. err_max(0) < TOLERANCE
else
   test_passed = unfilled(1) > 0_I4P
endif
if (mpih%myrank == 0) then
   if (test_passed) then
      print '(A)', 'TEST PASSED: periodic ghost (periodic leg: '//trim(str(is_periodic))//')'
   else
      print '(A)', 'TEST FAILED: periodic ghost (periodic leg: '//trim(str(is_periodic))//')'
   endif
endif
call mpih%finalize
if (.not.test_passed) error stop 1

contains
   function exact(b, i, j, k) result(fx)
   !< Return the periodic test function at the center of cell (i, j, k) of block b.
   integer(I4P), intent(in) :: b, i, j, k !< Block and cell indexes.
   real(R8P)                :: fx         !< Function value.
   real(R8P)                :: L(3)       !< Domain extents.

   L = realm%adam%grid%domain_emax - realm%adam%grid%domain_emin
   fx =        sin(2._R8P*PI*realm%adam%field%x_cell(i,b)/L(1)) + &
        2._R8P*sin(2._R8P*PI*realm%adam%field%y_cell(j,b)/L(2)) + &
        3._R8P*sin(2._R8P*PI*realm%adam%field%z_cell(k,b)/L(3))
   endfunction exact
endprogram test_periodic_ghost
