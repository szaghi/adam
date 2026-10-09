!< Unit test of a runtime regrid round trip: refine, then derefine back (issue #74, M5-P0 and P1).
program test_amr_regrid_roundtrip
!< Unit test of a runtime regrid round trip: refine, then derefine back (issue #74, M5-P0 and P1).
!<
!< **Why this test exists.** Runtime AMR (M5) calls `adam_object%amr_update` during the time loop: the tree adapts, the
!< interior data of the new blocks is prolonged from their parent (`field%refine*`) or restricted to it
!< (`field%derefine*`, the mean of the children), and the blocks are redistributed over the ranks. Until M5 nothing ran
!< it after initialisation, and no test ever derefined.
!<
!< **What it pins.** A realm refined uniformly to level 1 holds a field `f` in every cell, ghosts included (a perfect
!< ghost fill, so the transfers are tested alone). The block at the domain origin is refined once, then its children
!< are derefined; after each step every interior cell is compared with `f` at its centre, and the discrete integral
!< sum(q dV) over the ranks with its initial value. Fields: linear, quadratic, and a steep front from 1 to 0.01 (a tanh
!< of width 0.03, half a coarse cell, across the diagonal plane x + y + z = 0.45 on the octree, the line x + y = 0.3 on
!< the quadtree, through the refined block). On the octree front the monotonized-central slopes alone would put corner
!< children up to 0.10 outside the data, so below zero: the bound scaling of the conservative prolongation is what
!< keeps it positive. For both prolongations (`[amr] regrid_prolongation`):
!<
!< - `linear` (tensor linear, 1/4 and 3/4; the library default): a linear `f` exact after the refine and the round
!<   trip (asserted); not conservative, so on the quadratic `f` the integral drifts and the round trip does not
!<   return the parent (both asserted, the negative control of the conservative one);
!< - `conservative` (limited linear slopes, children's mean = parent, #74 D-M5-2): a linear `f` exact after the refine
!<   (asserted); for every field the round trip returns the initial state and the integral is unchanged after the
!<   refine and after the round trip, to round-off (asserted); on the front every child stays in (0.01, 1), the range of
!<   the data, so a positive variable stays positive (asserted);
!< - bookkeeping: after the refine the tree has 2^d - 1 more leaves (8 octree, 4 quadtree), after the derefine the
!<   original count, and the blocks summed over the ranks equal the leaves (asserted), on any number of ranks (the
!<   redistribution moves blocks between them, so the integral checks cover it).
!<
!< Configurations: octree (`ratio = 8`, 8^3 cells per block) and quadtree (`ratio = 4`, `nk = 1`).
!<
!< Usage: `mpirun -np N exe/test_amr_regrid_roundtrip`.

use :: adam_common_library, only : realm_object
use :: adam_parameters,     only : AMR_PROLONGATION_CONSERVATIVE, AMR_PROLONGATION_LINEAR, TO_BE_DEREFINED, TO_BE_REFINED, &
                                   TO_NOT_TOUCH
use :: adam_mpih_global,    only : mpih
use :: mpi
use :: penf,                only : I4P, R8P, str

implicit none

real(R8P),    parameter :: TOLERANCE=1.e-13_R8P !< Round-off tolerance, relative to the field scale.
integer(I4P), parameter :: F_LINEAR=1_I4P       !< Linear field.
integer(I4P), parameter :: F_QUADRATIC=2_I4P    !< Quadratic field.
integer(I4P), parameter :: F_FRONT=3_I4P        !< Steep front from 1 to 0.01.
character(9), parameter :: F_NAME(3)=['linear   ', 'quadratic', 'front    '] !< Field names.
logical              :: test_passed          !< Aggregate pass flag.
type(realm_object), allocatable :: realm     !< Realm under test (one per configuration).
real(R8P), allocatable :: q(:,:,:,:,:)       !< Field (1, i, j, k, b).
integer(I4P)           :: ngc, ni, nj, nk    !< Grid dimensions of the configuration.
integer(I4P)           :: nd                 !< Refined axes, 3 octree, 2 quadtree.
integer(I4P)           :: field_kind         !< Field kind, F_*.
real(R8P)              :: scale              !< Field scale.
integer(I4P)           :: ierr               !< MPI error status.
integer(I4P)           :: fk                 !< Field kind counter.

call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)
test_passed = .true.
do fk=F_LINEAR, F_FRONT
   call check_configuration(label='octree,   linear,       ', ratio=8_I4P, nk_=8_I4P, kind_=fk, &
                            prolongation=AMR_PROLONGATION_LINEAR)
   call check_configuration(label='octree,   conservative, ', ratio=8_I4P, nk_=8_I4P, kind_=fk, &
                            prolongation=AMR_PROLONGATION_CONSERVATIVE)
   call check_configuration(label='quadtree, linear,       ', ratio=4_I4P, nk_=1_I4P, kind_=fk, &
                            prolongation=AMR_PROLONGATION_LINEAR)
   call check_configuration(label='quadtree, conservative, ', ratio=4_I4P, nk_=1_I4P, kind_=fk, &
                            prolongation=AMR_PROLONGATION_CONSERVATIVE)
enddo
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
   subroutine check_configuration(label, ratio, nk_, kind_, prolongation)
   !< Refine the origin block, derefine its children, check the field, the integral and the bookkeeping after each step.
   character(*), intent(in)  :: label          !< Configuration label.
   integer(I4P), intent(in)  :: ratio          !< Tree ratio.
   integer(I4P), intent(in)  :: nk_            !< Cells per block along z.
   integer(I4P), intent(in)  :: kind_          !< Field kind, F_*.
   integer(I4P), intent(in)  :: prolongation   !< Prolongation kind, AMR_PROLONGATION_*.
   character(:), allocatable :: filename       !< Input file name.
   real(R8P)                 :: dmin(3)        !< Domain origin.
   real(R8P)                 :: err_refine     !< Max |q - f| / scale after the refine.
   real(R8P)                 :: err_trip       !< Max |q - f| / scale after the round trip.
   real(R8P)                 :: total(0:2)     !< Integral sum(q dV): initial, refined, round trip.
   real(R8P)                 :: drift(1:2)     !< |integral - initial| / |initial|: refined, round trip.
   real(R8P)                 :: range_(0:1,2)  !< Min and max of q: initial, after the refine.
   integer(I4P)              :: leaves(0:2)    !< Blocks summed over the ranks: initial, refined, round trip.
   integer(I4P)              :: b              !< Block counter.
   logical                   :: conservative   !< Conservative prolongation.
   logical                   :: ok             !< Configuration pass flag.

   nk = nk_ ; field_kind = kind_
   nd = merge(3_I4P, 2_I4P, ratio == 8_I4P)
   conservative = prolongation == AMR_PROLONGATION_CONSERVATIVE
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
   select case(field_kind)
   case(F_LINEAR)    ; scale = 4._R8P
   case(F_QUADRATIC) ; scale = 10._R8P
   case(F_FRONT)     ; scale = 1._R8P
   endselect
   dmin = realm%adam%grid%domain_emin
   call fill(all_cells=.true.)
   leaves(0) = total_blocks()
   total(0) = integral()
   range_(0,:) = value_range()

   ! refine the block at the domain origin
   realm%adam%field%refinements_needed = [(TO_NOT_TOUCH, b=1, realm%adam%field%blocks_number)]
   do b=1, realm%adam%field%blocks_number
      if (all(abs(realm%adam%field%emin(:,b) - dmin) < 1.e-12_R8P)) realm%adam%field%refinements_needed(b) = TO_BE_REFINED
   enddo
   call realm%adam%amr_update(q=q, is_marked_by_field=.true., do_mpi_redistribute=.true., do_blocks_reorder=.false., &
                              prolongation=prolongation)
   leaves(1) = total_blocks()
   total(1) = integral()
   err_refine = max_error()
   range_(1,:) = value_range()

   ! derefine the children (the level-2 blocks), back to the uniform tree; ghosts refreshed exactly first
   call fill(all_cells=.false.)
   realm%adam%field%refinements_needed = [(TO_NOT_TOUCH, b=1, realm%adam%field%blocks_number)]
   do b=1, realm%adam%field%blocks_number
      if (realm%adam%field%coordinates(4,b) == 2_I4P) realm%adam%field%refinements_needed(b) = TO_BE_DEREFINED
   enddo
   call realm%adam%amr_update(q=q, is_marked_by_field=.true., do_mpi_redistribute=.true., do_blocks_reorder=.false., &
                              prolongation=prolongation)
   leaves(2) = total_blocks()
   total(2) = integral()
   err_trip = max_error()
   drift = abs(total(1:2) - total(0)) / abs(total(0))

   ok = leaves(1) == leaves(0) + (ratio - 1_I4P) .and. leaves(2) == leaves(0)
   if (field_kind == F_LINEAR) ok = ok .and. err_refine <= TOLERANCE
   if (conservative) then
      ok = ok .and. err_trip <= TOLERANCE .and. all(drift <= TOLERANCE)
      ! the front's data (ghosts included) lie in (0.01, 1): a child outside left the range of its stencil
      if (field_kind == F_FRONT) ok = ok .and. range_(1,1) >= 0.01_R8P .and. range_(1,2) <= 1._R8P
   else
      if (field_kind == F_LINEAR) ok = ok .and. err_trip <= TOLERANCE
      ! negative control: the linear prolongation is not conservative on curved data
      if (field_kind == F_QUADRATIC) ok = ok .and. drift(1) > 1.e3_R8P * TOLERANCE .and. err_trip > 1.e3_R8P * TOLERANCE
   endif
   ok = ok .and. err_refine == err_refine .and. err_trip == err_trip .and. all(drift == drift) ! finite
   test_passed = test_passed .and. ok
   if (mpih%myrank == 0) then
      print '(A)', '   '//label//F_NAME(field_kind)//': leaves '//trim(str(leaves(0), .true.))//' -> '//                 &
                   trim(str(leaves(1), .true.))//' -> '//trim(str(leaves(2), .true.))//', max |q - f| / scale refine '// &
                   trim(str(err_refine, .true.))//' round trip '//trim(str(err_trip, .true.))//', integral drift '//     &
                   trim(str(drift(1), .true.))//' '//trim(str(drift(2), .true.))//', range '//                          &
                   trim(str(range_(1,1)))//' '//trim(str(range_(1,2)))//' (initial '//               &
                   trim(str(range_(0,1)))//' '//trim(str(range_(0,2)))//'): '//merge('PASS', 'FAIL', ok)
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
                                     realm%adam%field%z_cell(k,b)], field_kind)
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
                                                        realm%adam%field%z_cell(k,b)], field_kind)) / scale)
            enddo
         enddo
      enddo
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, err, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   endfunction max_error

   function integral() result(total)
   !< Sum over the ranks of q dV on the interior cells, dV the cell size along the refined axes.
   real(R8P)    :: total      !< Integral.
   integer(I4P) :: b          !< Counter.

   total = 0._R8P
   do b=1, realm%adam%field%blocks_number
      total = total + sum(q(1,1:ni,1:nj,1:nk,b)) * product(realm%adam%field%dxyz(1:nd,b))
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, total, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
   endfunction integral

   function value_range() result(range_)
   !< Min and max over the ranks of q on the interior cells.
   real(R8P)    :: range_(2) !< Min and max.
   integer(I4P) :: b         !< Counter.

   range_ = [huge(1._R8P), -huge(1._R8P)]
   do b=1, realm%adam%field%blocks_number
      range_(1) = min(range_(1), minval(q(1,1:ni,1:nj,1:nk,b)))
      range_(2) = max(range_(2), maxval(q(1,1:ni,1:nj,1:nk,b)))
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, range_(1), 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE, range_(2), 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   endfunction value_range

   function total_blocks() result(n)
   !< Blocks summed over the ranks; must equal the tree leaves.
   integer(I4P) :: n !< Count.

   n = realm%adam%field%blocks_number
   call MPI_ALLREDUCE(MPI_IN_PLACE, n, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   endfunction total_blocks

   pure function exact(xyz, kind_) result(f)
   !< Exact field at a point: linear, quadratic, or a steep front across x + y + z = 0.45 (octree) or x + y = 0.3.
   real(R8P),    intent(in) :: xyz(3) !< Point.
   integer(I4P), intent(in) :: kind_  !< Field kind, F_*.
   real(R8P)                :: f      !< Value.

   select case(kind_)
   case(F_FRONT)
      if (nd == 3_I4P) then
         f = 0.505_R8P + 0.495_R8P * tanh((0.45_R8P - xyz(1) - xyz(2) - xyz(3)) / 0.03_R8P)
      else
         f = 0.505_R8P + 0.495_R8P * tanh((0.3_R8P - xyz(1) - xyz(2)) / 0.03_R8P)
      endif
   case default
      f = 1._R8P + 0.3_R8P * xyz(1) + 0.7_R8P * xyz(2) + 0.5_R8P * xyz(3)
      if (kind_ == F_QUADRATIC) f = f + 2._R8P * xyz(1)**2 + 1.5_R8P * xyz(2)**2 + xyz(1) * xyz(2)
   endselect
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
