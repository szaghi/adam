!< ADAM, test of the analytical rectangle's distance function against the exact signed distance (issue #143).
program test_ib_rectangle_phi
!< ADAM, test of the analytical rectangle's distance function against the exact signed distance (issue #143).
!<
!< `ib_object%compute_phi` on an `analytical_rectangle` solid must give the exact signed distance from the rectangle's
!< boundary, positive inside the solid and negative outside: the convention of the sphere and of the circle, whose phi is
!< `-(|p - c| - R)`. The rectangle lies in the plane normal to its axis; `edge(1)`, `edge(2)` are its side lengths along
!< the first and the second of the two other axes in x, y, z order (axis x: y, z; axis y: x, z; axis z: x, y), and phi
!< does not depend on the coordinate along the axis.
!<
!< The cell coordinates are set by hand, not by a grid: `compute_phi` evaluates phi at `(x_cell(i), y_cell(j),
!< z_cell(k))`, so the two in-plane coordinate lists below give the tensor product of the probe abscissae: below, on and
!< just inside the lower edge, the centre, just inside, on and beyond the upper edge, along both in-plane directions.
!< That covers points inside near each face, inside at the centre, outside beside each face, outside past each corner
!< and on the boundary. The reference is computed independently of the kernel's formula: the minimum distance from the
!< four boundary segments, signed by an inside test. Two values worked by hand in the issue are checked explicitly.
!< Exit status 0 when every check passes, 1 otherwise.
use :: adam_field_object, only : field_object
use :: adam_grid_object,  only : grid_object
use :: adam_ib_object,    only : ib_object
use :: penf,              only : I4P, R8P

implicit none

real(R8P), parameter :: TOLERANCE = 1.e-12_R8P                                  !< Absolute tolerance on phi.
real(R8P), parameter :: CENTER(3) = [1.8_R8P, 0._R8P, 0.3_R8P]                  !< Centre, in-plane values reused.
real(R8P), parameter :: EDGE(2)   = [2.4_R8P, 0.4_R8P]                           !< Side lengths: [0.6, 3] x [-0.2, 0.2].
real(R8P), parameter :: P1(*)     = [-0.5_R8P, 0.5_R8P, 0.6_R8P, 0.65_R8P, 1.8_R8P, &
                                     2.9_R8P, 3._R8P, 3.2_R8P, 4.5_R8P]          !< First in-plane abscissae.
real(R8P), parameter :: P2(*)     = [-1._R8P, -0.25_R8P, -0.2_R8P, -0.15_R8P, 0._R8P, &
                                     0.1_R8P, 0.2_R8P, 0.3_R8P, 1.2_R8P]         !< Second in-plane abscissae.
real(R8P), parameter :: PA(2)     = [-7._R8P, 11._R8P]                           !< Abscissae along the axis.
integer(I4P)         :: failed                                                   !< Failed checks.
character(1)         :: axis                                                     !< Rectangle axis.
integer(I4P)         :: a                                                        !< Counter.

failed = 0
do a=1, 3
   axis = 'xyz'(a:a)
   call check_axis(axis=axis)
enddo
call check_issue_values

if (failed > 0) then
   print '(A,I0,A)', 'test_ib_rectangle_phi: FAILED, ', failed, ' check(s)'
   stop 1
endif
print '(A)', 'test_ib_rectangle_phi: all checks passed'

contains
   subroutine setup(axis, x, y, z, field, grid, ib)
   !< Build a one-block field holding the given cell coordinates and an IB of one rectangle with the test's geometry.
   character(1),       intent(in)  :: axis  !< Rectangle axis.
   real(R8P),          intent(in)  :: x(:)  !< x cell coordinates.
   real(R8P),          intent(in)  :: y(:)  !< y cell coordinates.
   real(R8P),          intent(in)  :: z(:)  !< z cell coordinates.
   type(field_object), intent(out) :: field !< Field holding the coordinates.
   type(grid_object),  intent(out) :: grid  !< Grid holding the sizes.
   type(ib_object),    intent(out) :: ib    !< IB, one rectangle.

   grid%ngc = 0
   grid%ni = size(x, kind=I4P)
   grid%nj = size(y, kind=I4P)
   grid%nk = size(z, kind=I4P)
   field%blocks_number = 1
   allocate(field%x_cell(1:grid%ni, 1:1), source=reshape(x, [grid%ni, 1_I4P]))
   allocate(field%y_cell(1:grid%nj, 1:1), source=reshape(y, [grid%nj, 1_I4P]))
   allocate(field%z_cell(1:grid%nk, 1:1), source=reshape(z, [grid%nk, 1_I4P]))
   ib%solids_number = 1
   allocate(ib%definition(1), ib%rectangle(1))
   ib%definition(1) = 'analytical_rectangle'
   ib%rectangle(1)%axis = axis
   ib%rectangle(1)%edge = EDGE
   select case(axis) ! the in-plane centre coordinates are CENTER(1:2) whatever the axis; the axial one is CENTER(3)
   case('x')
      ib%rectangle(1)%center = [CENTER(3), CENTER(1), CENTER(2)]
   case('y')
      ib%rectangle(1)%center = [CENTER(1), CENTER(3), CENTER(2)]
   case('z')
      ib%rectangle(1)%center = [CENTER(1), CENTER(2), CENTER(3)]
   endselect
   allocate(ib%phi(1:2, 1:grid%ni, 1:grid%nj, 1:grid%nk, 1:1))
   ib%phi = -1._R8P
   endsubroutine setup

   subroutine check_axis(axis)
   !< Check phi on the tensor product of the probe abscissae, in the plane normal to `axis`, at two axial positions.
   character(1), intent(in) :: axis       !< Rectangle axis.
   type(field_object)       :: field      !< Field holding the coordinates.
   type(grid_object)        :: grid       !< Grid holding the sizes.
   type(ib_object)          :: ib         !< IB, one rectangle.
   integer(I4P)             :: i, j, k    !< Counters.
   integer(I4P)             :: m, n       !< In-plane probe indexes.
   real(R8P)                :: phi        !< Computed phi.
   real(R8P)                :: phi_all    !< Computed phi, all-solids summary.
   real(R8P)                :: exact      !< Exact signed distance.

   select case(axis)
   case('x')
      call setup(axis=axis, x=PA, y=P1, z=P2, field=field, grid=grid, ib=ib)
   case('y')
      call setup(axis=axis, x=P1, y=PA, z=P2, field=field, grid=grid, ib=ib)
   case('z')
      call setup(axis=axis, x=P1, y=P2, z=PA, field=field, grid=grid, ib=ib)
   endselect
   call ib%compute_phi(field=field, grid=grid)
   do k=1, grid%nk
      do j=1, grid%nj
         do i=1, grid%ni
            select case(axis)
            case('x')
               m = j ; n = k
            case('y')
               m = i ; n = k
            case default
               m = i ; n = j
            endselect
            exact = exact_signed_distance(p=[P1(m), P2(n)])
            phi = ib%phi(1,i,j,k,1)
            phi_all = ib%phi(2,i,j,k,1)
            if (abs(phi - exact) > TOLERANCE .or. abs(phi_all - phi) > 0._R8P) then
               failed = failed + 1
               print '(A,A,A,2(F7.3,A),3(A,ES12.4))', 'FAIL axis ', axis, ': in-plane point (', P1(m), ',', P2(n), &
                     ')', ' phi ', phi, ' exact ', exact, ' summary ', phi_all
            endif
         enddo
      enddo
   enddo
   print '(A,A,A,I0,A)', 'axis ', axis, ': ', size(P1)*size(P2)*size(PA), ' points checked'
   endsubroutine check_axis

   subroutine check_issue_values
   !< The two values worked by hand in issue #143 (axis z): (0.5, 0.1) outside, 0.1 from the face x = 0.6, phi = -0.1;
   !< (0.65, 0.1) inside, 0.05 from the face x = 0.6, phi = +0.05.
   type(field_object) :: field     !< Field holding the coordinates.
   type(grid_object)  :: grid      !< Grid holding the sizes.
   type(ib_object)    :: ib        !< IB, one rectangle.
   real(R8P)          :: phi(2)    !< Computed phi.
   real(R8P), parameter :: EXPECTED(2) = [-0.1_R8P, 0.05_R8P] !< Expected phi.
   integer(I4P)       :: i         !< Counter.

   call setup(axis='z', x=[0.5_R8P, 0.65_R8P], y=[0.1_R8P], z=[0._R8P], field=field, grid=grid, ib=ib)
   call ib%compute_phi(field=field, grid=grid)
   phi = ib%phi(1,1:2,1,1,1)
   do i=1, 2
      if (abs(phi(i) - EXPECTED(i)) > TOLERANCE) then
         failed = failed + 1
         print '(A,F5.2,A,F8.4,A,F8.4)', 'FAIL issue #143 value at (', field%x_cell(i,1), ', 0.1): phi ', phi(i), &
               ' expected ', EXPECTED(i)
      endif
   enddo
   print '(A,2F9.4)', 'issue #143 values (0.5,0.1), (0.65,0.1): phi =', phi
   endsubroutine check_issue_values

   pure function exact_signed_distance(p) result(distance)
   !< Exact signed distance of the in-plane point `p` from the rectangle's boundary, positive inside: the minimum distance
   !< from the four boundary segments, signed by an inside test.
   real(R8P), intent(in) :: p(2)     !< In-plane point.
   real(R8P)             :: distance !< Signed distance.
   real(R8P)             :: lo(2)    !< Lower corner.
   real(R8P)             :: hi(2)    !< Upper corner.

   lo = CENTER(1:2) - EDGE * 0.5_R8P
   hi = CENTER(1:2) + EDGE * 0.5_R8P
   distance = min(segment_distance(p, [lo(1), lo(2)], [hi(1), lo(2)]), &
                  segment_distance(p, [hi(1), lo(2)], [hi(1), hi(2)]), &
                  segment_distance(p, [hi(1), hi(2)], [lo(1), hi(2)]), &
                  segment_distance(p, [lo(1), hi(2)], [lo(1), lo(2)]))
   if (any(p < lo) .or. any(p > hi)) distance = -distance
   endfunction exact_signed_distance

   pure function segment_distance(p, a, b) result(distance)
   !< Distance of the point `p` from the segment `[a, b]`.
   real(R8P), intent(in) :: p(2)     !< Point.
   real(R8P), intent(in) :: a(2)     !< Segment start.
   real(R8P), intent(in) :: b(2)     !< Segment end.
   real(R8P)             :: distance !< Distance.
   real(R8P)             :: t        !< Parameter of the projection, clipped to the segment.

   t = max(0._R8P, min(1._R8P, dot_product(p - a, b - a) / dot_product(b - a, b - a)))
   distance = norm2(p - a - t * (b - a))
   endfunction segment_distance
endprogram test_ib_rectangle_phi
