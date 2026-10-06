!< Unit test of the 2:1 seam ghost exchange on octrees and quadtrees (issue #46).
program test_quadtree_seam_ghost
!< Unit test of the 2:1 seam ghost exchange on octrees and quadtrees (issue #46).
!<
!< **Why this test exists.** A quadtree (`ratio = 4`) refines x and y only: its blocks span the whole domain along z at
!< every level, so a 2:1 seam is 2:1 in x and y and 1:1 in z. The seam machinery was written for the octree (fine index
!< `2 i + delta` on every axis, 8-cell restriction, tricubic coarse->fine fill); on a quadtree it mixes cells of
!< different z, and the coarse cells beside a seam picked up a spurious z dependence (#46).
!<
!< **What it pins.** A realm refined uniformly to level 1, then once more on the block at the domain origin (a 2:1 seam),
!< holds the linear field `f = 1 + 0.3 x + 0.7 y + 1.1 z`; after `update_ghost_local` + `update_ghost_mpi` every ghost cell
!< inside the domain (the others belong to the boundary conditions, not exchanged) must hold `f` at its own centre:
!< the same-level copy, the fine->coarse restriction (mean of the fine cells) and the coarse->fine tricubic fill are all
!< exact on a linear field. Configurations: octree (`ratio = 8`, 8^3 cells per block), quadtree with `nk = 4` and with
!< `nk = 1` (true 2-D). Ghosts are reported by class (face, edge, corner); every class is asserted.
!<
!< Usage: `mpirun -np N exe/test_quadtree_seam_ghost`.

use :: adam_common_library, only : realm_object
use :: adam_parameters,     only : TO_BE_REFINED, TO_NOT_TOUCH
use :: adam_mpih_global,    only : mpih
use :: mpi
use :: penf,                only : I4P, R8P, str

implicit none

real(R8P), parameter :: SENTINEL=-1.e30_R8P  !< Value marking a ghost cell not filled by the exchange.
real(R8P), parameter :: TOLERANCE=1.e-12_R8P !< Round-off tolerance, relative to the field scale.
real(R8P), parameter :: SCALE=6._R8P         !< Field scale, max |f| over the domain [0,1] x [0,2] x [0,3].
logical              :: test_passed          !< Aggregate pass flag.

call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)
test_passed = .true.
call check_configuration(label='octree,   nk 8', ratio=8_I4P, nk=8_I4P)
call check_configuration(label='quadtree, nk 4', ratio=4_I4P, nk=4_I4P)
call check_configuration(label='quadtree, nk 1', ratio=4_I4P, nk=1_I4P)
if (mpih%myrank == 0) then
   if (test_passed) then
      print '(A)', 'test_quadtree_seam_ghost: PASSED'
   else
      print '(A)', 'test_quadtree_seam_ghost: FAILED'
   endif
endif
call mpih%finalize
if (.not.test_passed) error stop 1

contains
   subroutine check_configuration(label, ratio, nk)
   !< Build the 2:1 seam, exchange the ghosts of the linear field, check every ghost inside the domain.
   character(*), intent(in)  :: label            !< Configuration label.
   integer(I4P), intent(in)  :: ratio            !< Tree ratio.
   integer(I4P), intent(in)  :: nk               !< Cells per block along z.
   type(realm_object)        :: realm            !< Realm under test.
   character(:), allocatable :: filename         !< Input file name.
   real(R8P),    allocatable :: q(:,:,:,:,:)     !< Field (1, i, j, k, b).
   real(R8P)                 :: err_max(0:3)     !< Max |q - f| / scale per class: 0 interior, 1 face, 2 edge, 3 corner.
   integer(I4P)              :: unfilled(0:3)    !< Sentinel ghosts left per class.
   integer(I4P)              :: checked(0:3)     !< Ghosts checked per class.
   real(R8P)                 :: xyz(3)           !< Cell centre.
   real(R8P)                 :: dmin(3), dmax(3) !< Domain extents.
   integer(I4P)              :: ngc, ni, nj      !< Grid dimensions.
   integer(I4P)              :: b, i, j, k, cls  !< Counters, cell class.
   integer(I4P)              :: ierr             !< MPI error status.
   logical                   :: inside(3)        !< Cell inside the block interior per direction.
   logical                   :: ok               !< Configuration pass flag.

   filename = 'test_quadtree_seam_ghost-r'//trim(str(ratio, .true.))//'-nk'//trim(str(nk, .true.))//'.ini'
   if (mpih%myrank == 0) call write_input(filename=filename, ratio=ratio, nk=nk)
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   call realm%initialize(filename=filename, memory_avail=1._R8P, nv=1_I4P)
   ngc = realm%adam%grid%ngc ; ni = realm%adam%grid%ni ; nj = realm%adam%grid%nj
   allocate(q(1:1,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:realm%adam%field%nb))
   q = 0._R8P
   call realm%adam%refine_uniform(refinement_levels=realm%adam%tree%iu_ref_levels, q=q, do_mpi_redistribute=.true., &
                                  do_blocks_reorder=.false.)
   ! one more level on the block at the domain origin: a 2:1 seam
   dmin = realm%adam%grid%domain_emin
   dmax = realm%adam%grid%domain_emax
   realm%adam%field%refinements_needed = [(TO_NOT_TOUCH, b=1, realm%adam%field%blocks_number)]
   do b=1, realm%adam%field%blocks_number
      if (all(abs(realm%adam%field%emin(:,b) - dmin) < 1.e-12_R8P)) realm%adam%field%refinements_needed(b) = TO_BE_REFINED
   enddo
   call realm%adam%amr_update(q=q, is_marked_by_field=.true., do_mpi_redistribute=.true., do_blocks_reorder=.false.)

   q = SENTINEL
   do b=1, realm%adam%field%blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               q(1,i,j,k,b) = exact([realm%adam%field%x_cell(i,b), realm%adam%field%y_cell(j,b), &
                                     realm%adam%field%z_cell(k,b)])
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
               inside = [i >= 1 .and. i <= ni, j >= 1 .and. j <= nj, k >= 1 .and. k <= nk]
               cls = count(.not.inside)
               if (cls == 0) cycle
               xyz = [realm%adam%field%x_cell(i,b), realm%adam%field%y_cell(j,b), realm%adam%field%z_cell(k,b)]
               if (any(xyz < dmin) .or. any(xyz > dmax)) cycle ! a boundary-condition ghost
               checked(cls) = checked(cls) + 1_I4P
               if (q(1,i,j,k,b) == SENTINEL) then
                  unfilled(cls) = unfilled(cls) + 1_I4P
               else
                  err_max(cls) = max(err_max(cls), abs(q(1,i,j,k,b) - exact(xyz)) / SCALE)
               endif
            enddo
         enddo
      enddo
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, err_max, 4, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE, unfilled, 4, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE, checked, 4, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   ok = all(err_max(1:3) <= TOLERANCE) .and. all(unfilled(1:3) == 0_I4P)
   test_passed = test_passed .and. ok
   if (mpih%myrank == 0) then
      print '(A)', '   '//label//': '//trim(str(realm%adam%tree%nodes_number, .true.))//' nodes; ghosts checked '// &
                   '(face, edge, corner) '//trim(str(checked(1:3)))//', unfilled '//trim(str(unfilled(1:3)))//  &
                   ', max |q - f| / scale '//trim(str(err_max(1:3)))
      if (ok) then
         print '(A)', '   '//label//': PASS'
      else
         print '(A)', '   '//label//': FAIL'
      endif
   endif
   endsubroutine check_configuration

   pure function exact(xyz) result(f)
   !< Exact value of the linear field at a point.
   real(R8P), intent(in) :: xyz(3) !< Point.
   real(R8P)             :: f      !< Value.

   f = 1._R8P + 0.3_R8P * xyz(1) + 0.7_R8P * xyz(2) + 1.1_R8P * xyz(3)
   endfunction exact

   subroutine write_input(filename, ratio, nk)
   !< Write the input of one configuration: the sections realm_object%initialize reads, domain [0,1] x [0,2] x [0,3].
   character(*), intent(in) :: filename !< Input file name.
   integer(I4P), intent(in) :: ratio    !< Tree ratio.
   integer(I4P), intent(in) :: nk       !< Cells per block along z.
   integer(I4P)             :: u        !< File unit.

   open(newunit=u, file=filename, action='write', status='replace')
   write(u, '(A)') '[IO]'
   write(u, '(A)') 'output_basename        = test_quadtree_seam_ghost'
   write(u, '(A)') 'it_save                = 999999'
   write(u, '(A)') 'restart                = .false.'
   write(u, '(A)') 'restart_basename       = test_quadtree_seam_ghost-restart'
   write(u, '(A)') 'restart_save           = 999999'
   write(u, '(A)') 'residuals_save         = 999999'
   write(u, '(A)') 'save_memory_status     = .false.'
   write(u, '(A)') 'save_residual_fields   = .false.'
   write(u, '(A)') 'save_curl_fields       = .false.'
   write(u, '(A)') 'save_divergence_fields = .false.'
   write(u, '(A)') 'save_gradient_fields   = .false.'
   write(u, '(A)') 'save_laplacian_fields  = .false.'
   write(u, '(A)') '[grid]'
   write(u, '(A)') 'ni     = 8'
   write(u, '(A)') 'nj     = 8'
   write(u, '(A)') 'nk     = '//trim(str(nk, .true.))
   write(u, '(A)') 'ngc    = 3'
   write(u, '(A)') 'emin_x = 0.0'
   write(u, '(A)') 'emin_y = 0.0'
   write(u, '(A)') 'emin_z = 0.0'
   write(u, '(A)') 'emax_x = 1.0'
   write(u, '(A)') 'emax_y = 2.0'
   write(u, '(A)') 'emax_z = 3.0'
   write(u, '(A)') 'null_x = .false.'
   write(u, '(A)') 'null_y = .false.'
   if (nk == 1_I4P) then
      write(u, '(A)') 'null_z = .true.'
   else
      write(u, '(A)') 'null_z = .false.'
   endif
   write(u, '(A)') '[amr]'
   write(u, '(A)') 'max_level      = 2'
   write(u, '(A)') 'ratio          = '//trim(str(ratio, .true.))
   write(u, '(A)') 'iu_ref_levels  = 1'
   write(u, '(A)') 'i_prune        = 0'
   write(u, '(A)') 'j_prune        = 0'
   write(u, '(A)') 'k_prune        = 0'
   write(u, '(A)') 'l_prune        = -1'
   write(u, '(A)') 'frequency      = 999999'
   write(u, '(A)') 'iters          = 1'
   write(u, '(A)') 'markers_number = 0'
   write(u, '(A)') '[field]'
   write(u, '(A)') 'nv = 1'
   write(u, '(A)') '[runge_kutta]'
   write(u, '(A)') 'scheme = runge-kutta-ssp-33'
   write(u, '(A)') '[weno]'
   write(u, '(A)') 'scheme              = weno-u-5'
   write(u, '(A)') 'ror_number          = 0'
   write(u, '(A)') 'ror_threshold       = 0.9'
   write(u, '(A)') 'ror_vars_number     = 0'
   write(u, '(A)') 'enable_ror_stats    = .false.'
   write(u, '(A)') 'ib_reduction_extent = 0'
   write(u, '(A)') 'ib_reduced_order    = 2'
   write(u, '(A)') '[linear-algebra]'
   write(u, '(A)') 'smoothing         = gauss-seidel'
   write(u, '(A)') 'iterations_init   = 3'
   write(u, '(A)') 'iterations_coarse = 10'
   write(u, '(A)') 'iterations_fine   = 3'
   write(u, '(A)') 'iterations        = 10'
   write(u, '(A)') 'tolerance         = 1.e-20'
   write(u, '(A)') '[fdv]'
   write(u, '(A)') 'fdv_scheme = fd'
   write(u, '(A)') 'fdv_order  = 2'
   write(u, '(A)') '[solids]'
   write(u, '(A)') 'number    = 0'
   write(u, '(A)') 'n_eikonal = 0'
   write(u, '(A)') '[slices]'
   write(u, '(A)') 'slices_number = 0'
   close(u)
   endsubroutine write_input
endprogram test_quadtree_seam_ghost
