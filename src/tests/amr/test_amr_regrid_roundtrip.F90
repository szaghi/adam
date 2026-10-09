!< Unit test of a runtime regrid round trip: refine, then derefine back (issue #74, M5-P0).
program test_amr_regrid_roundtrip
!< Unit test of a runtime regrid round trip: refine, then derefine back (issue #74, M5-P0).
!<
!< **Why this test exists.** Runtime AMR (M5) calls `adam_object%amr_update` during the time loop: the tree adapts, the
!< interior data of the new blocks is prolonged from their parent (`field%refine*`) or restricted to it
!< (`field%derefine*`, the mean of the children), and the blocks are redistributed over the ranks. Until M5 nothing ran
!< it after initialisation, and no test ever derefined.
!<
!< **What it pins.** A realm refined uniformly to level 1 holds a field `f` in every cell, ghosts included (a perfect
!< ghost fill, so the transfers are tested alone). The block at the domain origin is refined once, then its children
!< are derefined; after each step every interior cell is compared with `f` at its centre:
!<
!< - linear `f`: the prolongation and the restriction are both exact, so the refined state and the round trip are exact
!<   (asserted to round-off);
!< - quadratic `f`: the library prolongation (tensor linear, 1/4 and 3/4) is second order and not conservative, so the
!<   children differ from `f` by O(h^2) and the round trip does not return the parent (reported; P1 adds a conservative
!<   prolongation for which it does);
!< - bookkeeping: after the refine the tree has 2^d - 1 more leaves (8 octree, 4 quadtree), after the derefine the
!<   original count, and the blocks summed over the ranks equal the leaves (asserted), on any number of ranks (the
!<   redistribution moves blocks between them).
!<
!< Configurations: octree (`ratio = 8`, 8^3 cells per block) and quadtree (`ratio = 4`, `nk = 1`).
!<
!< Usage: `mpirun -np N exe/test_amr_regrid_roundtrip`.

use :: adam_common_library, only : realm_object
use :: adam_parameters,     only : TO_BE_DEREFINED, TO_BE_REFINED, TO_NOT_TOUCH
use :: adam_mpih_global,    only : mpih
use :: mpi
use :: penf,                only : I4P, R8P, str

implicit none

real(R8P), parameter :: TOLERANCE=1.e-13_R8P !< Round-off tolerance, relative to the field scale.
logical              :: test_passed          !< Aggregate pass flag.
type(realm_object), allocatable :: realm     !< Realm under test (one per configuration).
real(R8P), allocatable :: q(:,:,:,:,:)       !< Field (1, i, j, k, b).
integer(I4P)           :: ngc, ni, nj, nk    !< Grid dimensions of the configuration.
logical                :: quadratic          !< Quadratic field (else linear).
real(R8P)              :: scale              !< Field scale.
integer(I4P)           :: ierr               !< MPI error status.

call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)
test_passed = .true.
call check_configuration(label='octree,   linear   ', ratio=8_I4P, nk_=8_I4P, quadratic_=.false.)
call check_configuration(label='octree,   quadratic', ratio=8_I4P, nk_=8_I4P, quadratic_=.true.)
call check_configuration(label='quadtree, linear   ', ratio=4_I4P, nk_=1_I4P, quadratic_=.false.)
call check_configuration(label='quadtree, quadratic', ratio=4_I4P, nk_=1_I4P, quadratic_=.true.)
if (mpih%myrank == 0) then
   if (test_passed) then
      print '(A)', 'test_amr_regrid_roundtrip: PASSED'
   else
      print '(A)', 'test_amr_regrid_roundtrip: FAILED'
   endif
endif
call mpih%finalize
if (.not.test_passed) error stop 1

contains
   subroutine check_configuration(label, ratio, nk_, quadratic_)
   !< Refine the origin block, derefine its children, check the field and the bookkeeping after each step.
   character(*), intent(in)  :: label          !< Configuration label.
   integer(I4P), intent(in)  :: ratio          !< Tree ratio.
   integer(I4P), intent(in)  :: nk_            !< Cells per block along z.
   logical,      intent(in)  :: quadratic_     !< Quadratic field (else linear).
   character(:), allocatable :: filename       !< Input file name.
   real(R8P)                 :: dmin(3)        !< Domain origin.
   real(R8P)                 :: err_refine     !< Max |q - f| / scale after the refine.
   real(R8P)                 :: err_trip       !< Max |q - f| / scale after the round trip.
   integer(I4P)              :: leaves(0:2)    !< Blocks summed over the ranks: initial, refined, round trip.
   integer(I4P)              :: b              !< Block counter.
   logical                   :: ok             !< Configuration pass flag.

   nk = nk_ ; quadratic = quadratic_
   if (allocated(realm)) deallocate(realm)
   if (allocated(q)) deallocate(q)
   allocate(realm)
   filename = 'test_amr_regrid_roundtrip-r'//trim(str(ratio, .true.))//'-nk'//trim(str(nk, .true.))//'.ini'
   if (mpih%myrank == 0) call write_input(filename=filename, ratio=ratio, nk=nk)
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   call realm%initialize(filename=filename, memory_avail=1._R8P, nv=1_I4P)
   ngc = realm%adam%grid%ngc ; ni = realm%adam%grid%ni ; nj = realm%adam%grid%nj
   allocate(q(1:1,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:realm%adam%field%nb))
   q = 0._R8P
   call realm%adam%refine_uniform(refinement_levels=realm%adam%tree%iu_ref_levels, q=q, do_mpi_redistribute=.true., &
                                  do_blocks_reorder=.false.)
   scale = merge(10._R8P, 4._R8P, quadratic)
   dmin = realm%adam%grid%domain_emin
   call fill(all_cells=.true.)
   leaves(0) = total_blocks()

   ! refine the block at the domain origin
   realm%adam%field%refinements_needed = [(TO_NOT_TOUCH, b=1, realm%adam%field%blocks_number)]
   do b=1, realm%adam%field%blocks_number
      if (all(abs(realm%adam%field%emin(:,b) - dmin) < 1.e-12_R8P)) realm%adam%field%refinements_needed(b) = TO_BE_REFINED
   enddo
   call realm%adam%amr_update(q=q, is_marked_by_field=.true., do_mpi_redistribute=.true., do_blocks_reorder=.false.)
   leaves(1) = total_blocks()
   err_refine = max_error()

   ! derefine the children (the level-2 blocks), back to the uniform tree; ghosts refreshed exactly first
   call fill(all_cells=.false.)
   realm%adam%field%refinements_needed = [(TO_NOT_TOUCH, b=1, realm%adam%field%blocks_number)]
   do b=1, realm%adam%field%blocks_number
      if (realm%adam%field%coordinates(4,b) == 2_I4P) realm%adam%field%refinements_needed(b) = TO_BE_DEREFINED
   enddo
   call realm%adam%amr_update(q=q, is_marked_by_field=.true., do_mpi_redistribute=.true., do_blocks_reorder=.false.)
   leaves(2) = total_blocks()
   err_trip = max_error()

   ok = leaves(1) == leaves(0) + (ratio - 1_I4P) .and. leaves(2) == leaves(0)
   if (.not.quadratic) ok = ok .and. err_refine <= TOLERANCE .and. err_trip <= TOLERANCE
   ok = ok .and. err_refine == err_refine .and. err_trip == err_trip ! finite
   test_passed = test_passed .and. ok
   if (mpih%myrank == 0) then
      print '(A)', '   '//label//': leaves '//trim(str(leaves(0), .true.))//' -> '//trim(str(leaves(1), .true.))//   &
                   ' -> '//trim(str(leaves(2), .true.))//', max |q - f| / scale after refine '//trim(str(err_refine))// &
                   ', after the round trip '//trim(str(err_trip))//': '//merge('PASS', 'FAIL', ok)
   endif
   endsubroutine check_configuration

   subroutine fill(all_cells)
   !< Set the field to `f` in every cell of every block (ghosts included), or in the ghosts only.
   logical, intent(in) :: all_cells !< Interior too.
   integer(I4P)        :: b, i, j, k !< Counters.

   do b=1, realm%adam%field%blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (.not.all_cells .and. i >= 1 .and. i <= ni .and. j >= 1 .and. j <= nj .and. k >= 1 .and. k <= nk) cycle
               q(1,i,j,k,b) = exact([realm%adam%field%x_cell(i,b), realm%adam%field%y_cell(j,b), &
                                     realm%adam%field%z_cell(k,b)], quadratic)
            enddo
         enddo
      enddo
   enddo
   endsubroutine fill

   function max_error() result(err)
   !< Max over the ranks of |q - f| / scale on the interior cells.
   real(R8P)    :: err        !< Error.
   integer(I4P) :: b, i, j, k !< Counters.

   err = 0._R8P
   do b=1, realm%adam%field%blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               err = max(err, abs(q(1,i,j,k,b) - exact([realm%adam%field%x_cell(i,b), realm%adam%field%y_cell(j,b), &
                                                        realm%adam%field%z_cell(k,b)], quadratic)) / scale)
            enddo
         enddo
      enddo
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, err, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   endfunction max_error

   function total_blocks() result(n)
   !< Blocks summed over the ranks; must equal the tree leaves.
   integer(I4P) :: n !< Count.

   n = realm%adam%field%blocks_number
   call MPI_ALLREDUCE(MPI_IN_PLACE, n, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   endfunction total_blocks

   pure function exact(xyz, quad) result(f)
   !< Exact field at a point: linear, or with a quadratic part.
   real(R8P), intent(in) :: xyz(3)    !< Point.
   logical,   intent(in) :: quad   !< Add the quadratic part.
   real(R8P)             :: f      !< Value.

   f = 1._R8P + 0.3_R8P * xyz(1) + 0.7_R8P * xyz(2) + 0.5_R8P * xyz(3)
   if (quad) f = f + 2._R8P * xyz(1)**2 + 1.5_R8P * xyz(2)**2 + xyz(1) * xyz(2)
   endfunction exact

   subroutine write_input(filename, ratio, nk)
   !< Write the input of one configuration: the sections realm_object%initialize reads, domain [0,1]^3.
   character(*), intent(in) :: filename !< Input file name.
   integer(I4P), intent(in) :: ratio    !< Tree ratio.
   integer(I4P), intent(in) :: nk       !< Cells per block along z.
   integer(I4P)             :: u        !< File unit.

   open(newunit=u, file=filename, action='write', status='replace')
   write(u, '(A)') '[IO]'
   write(u, '(A)') 'output_basename        = test_amr_regrid_roundtrip'
   write(u, '(A)') 'it_save                = 999999'
   write(u, '(A)') 'restart                = .false.'
   write(u, '(A)') 'restart_basename       = test_amr_regrid_roundtrip-restart'
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
   write(u, '(A)') 'emax_y = 1.0'
   write(u, '(A)') 'emax_z = 1.0'
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
   write(u, '(A)') 'frequency      = 0'
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
endprogram test_amr_regrid_roundtrip
