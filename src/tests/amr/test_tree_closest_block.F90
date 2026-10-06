!< Unit test of the tree point lookup `tree%get_closest_block` (issue #54).
program test_tree_closest_block
!< Unit test of the tree point lookup `tree%get_closest_block` (issue #54).
!<
!< **Why this test exists.** The forest finds, on any rank, the peer cell of an inter-realm seam through
!< `tree%get_closest_block` (`build_seam_rows`, `register_inter_realm_seams`); the slices and the immersed boundary use it
!< too. Two defects were found while verifying #51: on a quadtree the lookup returned the wrong leaf (an inter-realm seam
!< stopped at initialization with a derived cell index past the block), and on a single-block tree it aborted (error
!< -111) instead of returning the root.
!<
!< **What it pins.** For binary trees (ratio 2), quadtrees (4) and octrees (8), uniformly refined to level 0, 1 and 2, with
!< `max_level` equal to the level and two levels above it (the lookup then walks up from a finest-level code that does not
!< exist), every leaf is probed at its centre, next to each of its faces and just past its lower corner: the returned code
!< must exist and decode (`morton_to_coordinates`, independent of the encoding under test) to the leaf's coordinates and
!< level. On a level-0 tree the leaf is the root, code -1.
!<
!< The realms are initialized from INI files the test writes (`test_tree_closest_block-r<ratio>-l<level>-m<max>.ini`) on
!< an anisotropic domain [0,1] x [0,2] x [0,3], so an axis mix-up cannot cancel out.
!<
!< Usage: `mpirun -np N exe/test_tree_closest_block`.

use :: adam_common_library, only : realm_object
use :: adam_mpih_global,    only : mpih
use :: mpi
use :: penf,                only : I4P, I8P, R8P, str

implicit none

integer(I4P), parameter :: RATIOS(3)=[2_I4P, 4_I4P, 8_I4P] !< Tree ratios: binary tree, quadtree, octree.
integer(I4P)            :: failures                         !< Failed probes, all configurations.
integer(I4P)            :: probes                           !< Probes, all configurations.
integer(I4P)            :: r, l, m                          !< Counters.

call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)
failures = 0_I4P
probes = 0_I4P
do r=1, size(RATIOS)
   do l=0, 2
      do m=l, l+2, 2
         call check_configuration(ratio=RATIOS(r), level=l, max_level=m)
      enddo
   enddo
enddo
if (mpih%myrank == 0) then
   if (failures == 0_I4P) then
      print '(A)', 'test_tree_closest_block: PASSED, '//trim(str(probes, .true.))//' probes'
   else
      print '(A)', 'test_tree_closest_block: FAILED, '//trim(str(failures, .true.))//' of '//trim(str(probes, .true.))// &
                   ' probes'
   endif
endif
call mpih%finalize
if (failures > 0_I4P) error stop 1

contains
   subroutine check_configuration(ratio, level, max_level)
   !< Initialize a realm with the given tree, refine it uniformly and probe every leaf.
   integer(I4P), intent(in)  :: ratio            !< Tree ratio.
   integer(I4P), intent(in)  :: level            !< Uniform refinement level.
   integer(I4P), intent(in)  :: max_level        !< Maximum refinement level.
   type(realm_object)        :: realm            !< Realm under test.
   character(:), allocatable :: filename         !< Input file name.
   real(R8P),    allocatable :: q(:,:,:,:,:)     !< Dummy field, refined with the tree.
   integer(I4P)              :: nb(3)            !< Leaves per direction.
   integer(I4P)              :: ijk(3), got(4)   !< Leaf coordinates, decoded coordinates and level.
   real(R8P)                 :: emin(3), emax(3) !< Leaf extents.
   real(R8P)                 :: point(3)         !< Probe.
   integer(I8P)              :: code             !< Returned code.
   integer(I4P)              :: config_failures  !< Failed probes of this configuration.
   integer(I4P)              :: i, j, k, p, a, s !< Counters.
   integer(I4P)              :: ierr             !< MPI error status.

   filename = 'test_tree_closest_block-r'//trim(str(ratio, .true.))//'-l'//trim(str(level, .true.))//'-m'// &
              trim(str(max_level, .true.))//'.ini'
   if (mpih%myrank == 0) call write_input(filename=filename, ratio=ratio, level=level, max_level=max_level)
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   call realm%initialize(filename=filename, memory_avail=1._R8P, nv=1_I4P)
   associate(ngc=>realm%adam%grid%ngc, ni=>realm%adam%grid%ni, nj=>realm%adam%grid%nj, nk=>realm%adam%grid%nk)
   allocate(q(1:1,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:realm%adam%field%nb))
   endassociate
   q = 0._R8P
   if (level > 0) call realm%adam%refine_uniform(refinement_levels=level, q=q, do_mpi_redistribute=.true., &
                                                 do_blocks_reorder=.false.)
   nb = 1_I4P
   nb(1) = 2_I4P**level
   if (ratio >= 4_I4P) nb(2) = 2_I4P**level
   if (ratio == 8_I4P) nb(3) = 2_I4P**level
   config_failures = 0_I4P
   do k=0, nb(3)-1
      do j=0, nb(2)-1
         do i=0, nb(1)-1
            ijk = [i, j, k]
            call realm%adam%grid%compute_metrics(coordinates=[ijk, level], emin=emin, emax=emax)
            do p=0, 7
               ! 0: centre; 1..6: next to the minus/plus face of each axis; 7: just past the lower corner
               point = 0.5_R8P * (emin + emax)
               if (p >= 1 .and. p <= 6) then
                  a = (p + 1) / 2
                  s = 2 * mod(p, 2) - 1
                  point(a) = point(a) + s * 0.499_R8P * (emax(a) - emin(a))
               elseif (p == 7) then
                  point = emin + 1.e-6_R8P * (emax - emin)
               endif
               code = realm%adam%tree%get_closest_block(grid=realm%adam%grid, point=point)
               call realm%adam%tree%morton_to_coordinates(code=code, i=got(1), j=got(2), k=got(3), l=got(4))
               probes = probes + 1_I4P
               if (.not.realm%adam%tree%has_code(code=code) .or. any(got /= [ijk, level])) then
                  config_failures = config_failures + 1_I4P
                  if (config_failures <= 3_I4P .and. mpih%myrank == 0) &
                     print '(A)', '   leaf '//trim(str(ijk))//' probe '//trim(str(p, .true.))//' point '// &
                                  trim(str(point))//': code '//trim(str(code))//' decodes to '//trim(str(got))
               endif
            enddo
         enddo
      enddo
   enddo
   failures = failures + config_failures
   if (mpih%myrank == 0) print '(A)', '   ratio '//trim(str(ratio, .true.))//', level '//trim(str(level, .true.))// &
                                      ', max_level '//trim(str(max_level, .true.))//': '//                       &
                                      trim(str(product(nb), .true.))//' leaves, '//                               &
                                      trim(str(config_failures, .true.))//' failed probes'
   endsubroutine check_configuration

   subroutine write_input(filename, ratio, level, max_level)
   !< Write the input of one configuration: the sections realm_object%initialize reads.
   character(*), intent(in) :: filename  !< Input file name.
   integer(I4P), intent(in) :: ratio     !< Tree ratio.
   integer(I4P), intent(in) :: level     !< Uniform refinement level.
   integer(I4P), intent(in) :: max_level !< Maximum refinement level.
   character(7)             :: null_y    !< Null y flag.
   character(7)             :: null_z    !< Null z flag.
   integer(I4P)             :: nj, nk    !< Cells along y, z.
   integer(I4P)             :: u         !< File unit.

   null_y = '.false.' ; nj = 8_I4P
   null_z = '.false.' ; nk = 8_I4P
   if (ratio <= 4_I4P) then
      null_z = '.true.'
      nk = 1_I4P
   endif
   if (ratio == 2_I4P) then
      null_y = '.true.'
      nj = 1_I4P
   endif
   open(newunit=u, file=filename, action='write', status='replace')
   write(u, '(A)') '[IO]'
   write(u, '(A)') 'output_basename        = test_tree_closest_block'
   write(u, '(A)') 'it_save                = 999999'
   write(u, '(A)') 'restart                = .false.'
   write(u, '(A)') 'restart_basename       = test_tree_closest_block-restart'
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
   write(u, '(A)') 'nj     = '//trim(str(nj, .true.))
   write(u, '(A)') 'nk     = '//trim(str(nk, .true.))
   write(u, '(A)') 'ngc    = 3'
   write(u, '(A)') 'emin_x = 0.0'
   write(u, '(A)') 'emin_y = 0.0'
   write(u, '(A)') 'emin_z = 0.0'
   write(u, '(A)') 'emax_x = 1.0'
   write(u, '(A)') 'emax_y = 2.0'
   write(u, '(A)') 'emax_z = 3.0'
   write(u, '(A)') 'null_x = .false.'
   write(u, '(A)') 'null_y = '//null_y
   write(u, '(A)') 'null_z = '//null_z
   write(u, '(A)') '[amr]'
   write(u, '(A)') 'max_level      = '//trim(str(max_level, .true.))
   write(u, '(A)') 'ratio          = '//trim(str(ratio, .true.))
   write(u, '(A)') 'iu_ref_levels  = '//trim(str(level, .true.))
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
endprogram test_tree_closest_block
