!< Standalone check of the PRISM PIC gathers against the deposits they pair with (no ADAM dependency).
!<
!< The deposits replicate the CPU backend (adam_prism_pic_object, adam_prism_cpu_object): B-spline charge with the optional
!< binomial filter, the Esirkepov directional charge change with its face prefix sum, the matched currents
!< (esirkepov_Q<N>, totorica: rho = P S, J_d = Q_d P_transverse F_d) and londrillo (A^-1 line solve with the checkerboard
!< constraint). Two gathers are compared:
!<   current  : the B-spline gather of field_weighting_model (order and filter as configured), the same for E and B;
!<   coherent : the adjoint of the deposit, E with the charge operator and B_c paired with the current component j_d it
!<              multiplies in v x B, gathered with the operator of J_d.
!< Metrics (unit charge, dx = 1, eps = 1):
!<   adjoint  : |sum_cells(deposit*F) - gather(F)| / sum|deposit*F| for a random F (rho and J_x, J_y, J_z);
!<   F_E      : electrostatic self-force, max over the sub-cell positions (periodic box, centered FD Poisson, F_ref = 1/(4 pi));
!<   F_B      : magnetostatic self-force v x B of the particle's own current (B = curl_h K J, K = (-lap_h)^-1).
!< Build and run: gfortran -O2 gather_consistency.F90 -o gather_consistency && ./gather_consistency
module gc_kernels
implicit none
integer, parameter :: R8P = selected_real_kind(15,307)
integer, parameter :: N = 32 !< Cells per axis.
real(R8P), parameter :: PI = acos(-1._R8P)
real(R8P), parameter :: FV1(4,4) = reshape([0.5_R8P,0._R8P,0._R8P,0._R8P,                            &
                                            7._R8P/12._R8P,-1._R8P/12._R8P,0._R8P,0._R8P,           &
                                            37._R8P/60._R8P,-8._R8P/60._R8P,1._R8P/60._R8P,0._R8P,   &
                                            533._R8P/840._R8P,-139._R8P/840._R8P,29._R8P/840._R8P,-3._R8P/840._R8P], [4,4])
real(R8P), parameter :: FD1(4,4) = reshape([0.5_R8P,0._R8P,0._R8P,0._R8P,                           &
                                            8._R8P/12._R8P,-1._R8P/12._R8P,0._R8P,0._R8P,           &
                                            45._R8P/60._R8P,-9._R8P/60._R8P,1._R8P/60._R8P,0._R8P,   &
                                            672._R8P/840._R8P,-168._R8P/840._R8P,32._R8P/840._R8P,-3._R8P/840._R8P], [4,4])

type :: config
   character(32) :: name = ''
   character(16) :: kind = 'poly'   !< poly (non-conserving), matched, londrillo.
   integer       :: order = 2       !< Particle B-spline order.
   logical       :: filter = .false. !< filter_deposition flag.
   integer       :: hs = 3          !< FD half stencil.
   real(R8P)     :: p(-7:7) = 0._R8P
   real(R8P)     :: qw(0:9) = 0._R8P
   integer       :: radius = 0, qfirst = 0, nq = 0
   real(R8P)     :: shift = 0._R8P  !< Odd esirkepov_Q<N>: base shape on the dual nodes.
   integer       :: gather_order = 2 !< field_weighting_model order of the current gather.
   real(R8P)     :: minv(N,N) = 0._R8P !< londrillo: inverse of the constrained face map.
endtype config

contains
   pure function bspl(order, r) result(w)
   integer,   intent(in) :: order
   real(R8P), intent(in) :: r
   real(R8P)             :: w, a
   a = abs(r) ; w = 0._R8P
   select case(order)
   case(0) ; if (r > -0.5_R8P .and. r <= 0.5_R8P) w = 1._R8P
   case(1) ; if (a <= 1._R8P) w = 1._R8P - a
   case(2)
      if (a <= 0.5_R8P) then ; w = 0.75_R8P - a*a
      elseif (a <= 1.5_R8P) then ; w = 0.5_R8P*(1.5_R8P-a)**2 ; endif
   case(3)
      if (a <= 1._R8P) then ; w = (4._R8P - 6._R8P*a*a + 3._R8P*a**3)/6._R8P
      elseif (a <= 2._R8P) then ; w = (2._R8P-a)**3/6._R8P ; endif
   case(4)
      if (a <= 0.5_R8P) then ; w = ((2.5_R8P-a)**4 - 5._R8P*(1.5_R8P-a)**4 + 10._R8P*(0.5_R8P-a)**4)/24._R8P
      elseif (a <= 1.5_R8P) then ; w = ((2.5_R8P-a)**4 - 5._R8P*(1.5_R8P-a)**4)/24._R8P
      elseif (a <= 2.5_R8P) then ; w = (2.5_R8P-a)**4/24._R8P ; endif
   case(5)
      if (a <= 1._R8P) then ; w = ((3._R8P-a)**5 - 6._R8P*(2._R8P-a)**5 + 15._R8P*(1._R8P-a)**5)/120._R8P
      elseif (a <= 2._R8P) then ; w = ((3._R8P-a)**5 - 6._R8P*(2._R8P-a)**5)/120._R8P
      elseif (a <= 3._R8P) then ; w = (3._R8P-a)**5/120._R8P ; endif
   case(6)
      if (a <= 0.5_R8P) then
         w = ((3.5_R8P-a)**6 - 7._R8P*(2.5_R8P-a)**6 + 21._R8P*(1.5_R8P-a)**6 - 35._R8P*(0.5_R8P-a)**6)/720._R8P
      elseif (a <= 1.5_R8P) then ; w = ((3.5_R8P-a)**6 - 7._R8P*(2.5_R8P-a)**6 + 21._R8P*(1.5_R8P-a)**6)/720._R8P
      elseif (a <= 2.5_R8P) then ; w = ((3.5_R8P-a)**6 - 7._R8P*(2.5_R8P-a)**6)/720._R8P
      elseif (a <= 3.5_R8P) then ; w = (3.5_R8P-a)**6/720._R8P ; endif
   endselect
   endfunction bspl

   pure function eff(order, r, filter) result(w)
   !< effective_bspline_weight of adam_prism_pic_object (filtered NGP included).
   integer,   intent(in) :: order
   real(R8P), intent(in) :: r
   logical,   intent(in) :: filter
   real(R8P)             :: w
   real(R8P), parameter  :: tol = 64._R8P*epsilon(1._R8P)
   if (filter) then
      if (order == 0) then
         w = 0._R8P
         if (r > -0.5_R8P+tol .and. r <= 0.5_R8P+tol) then ; w = 0.5_R8P
         elseif ((r > 0.5_R8P-tol .and. r <= 1.5_R8P+tol) .or. (r > -1.5_R8P-tol .and. r <= -0.5_R8P+tol)) then
            w = 0.25_R8P
         endif
      else
         w = 0.25_R8P*bspl(order,r+1._R8P) + 0.5_R8P*bspl(order,r) + 0.25_R8P*bspl(order,r-1._R8P)
      endif
   else
      w = bspl(order,r)
   endif
   endfunction eff

   pure function xc(i) result(x)
   integer, intent(in) :: i
   real(R8P)           :: x
   x = real(i,R8P) - 0.5_R8P
   endfunction xc

   subroutine set_matched(c, q)
   !< modified_kernels / totorica_kernels of PRISM (q = 0: totorica with c%hs).
   type(config), intent(inout) :: c
   integer,      intent(in)    :: q
   integer :: m
   c%kind = 'matched' ; c%p = 0._R8P ; c%qw = 0._R8P ; c%shift = 0._R8P
   select case(q)
   case(0)
      do m=1,c%hs
         c%p(1-m) = c%p(1-m) + 0.5_R8P*FV1(m,c%hs) ; c%p(m)   = c%p(m)   + 0.5_R8P*FV1(m,c%hs)
         c%p(-m)  = c%p(-m)  + 0.5_R8P*FV1(m,c%hs) ; c%p(m-1) = c%p(m-1) + 0.5_R8P*FV1(m,c%hs)
      enddo
      c%qw(0:1) = 0.5_R8P ; c%radius = c%hs ; c%qfirst = -1 ; c%nq = 2
   case(5)
      c%p(-5:4) = real([1,2,3,49,1985,1985,49,3,2,1],R8P)/4080._R8P
      c%qw(0:4) = real([1,10,46,10,1],R8P)/68._R8P
      c%radius = 6 ; c%qfirst = -3 ; c%nq = 5 ; c%shift = 0.5_R8P
   case(6)
      c%p(-5:5) = real([1,1,1,118,1378,2522,1378,118,1,1,1],R8P)/5520._R8P
      c%qw(0:5) = real([1,9,36,36,9,1],R8P)/92._R8P
      c%radius = 5 ; c%qfirst = -3 ; c%nq = 6
   case(8)
      c%p(-6:6) = real([1,1,1,152,1167,3327,4622,3327,1167,152,1,1,1],R8P)/13920._R8P
      c%qw(0:7) = real([1,9,36,70,70,36,9,1],R8P)/232._R8P
      c%radius = 6 ; c%qfirst = -4 ; c%nq = 8
   case(10)
      c%p(-7:7) = real([1,1,1,166,1111,3487,6567,8052,6567,3487,1111,166,1,1,1],R8P)/30720._R8P
      c%qw(0:9) = real([1,9,36,84,126,126,84,36,9,1],R8P)/512._R8P
      c%radius = 7 ; c%qfirst = -5 ; c%nq = 10
   endselect
   endsubroutine set_matched

   subroutine set_londrillo(c)
   !< Constrained face map of ensure_londrillo_current_solver_cache (odd mirror at the line ends), inverted densely.
   type(config), intent(inout) :: c
   real(R8P) :: a(N,N)
   integer   :: face, m, col
   real(R8P) :: sgn
   c%kind = 'londrillo' ; a = 0._R8P
   do face=1,N-1
      do m=1,c%hs
         call mirror(face+m, col, sgn)   ; a(face,col) = a(face,col) + FV1(m,c%hs)*sgn
         call mirror(face+1-m, col, sgn) ; a(face,col) = a(face,col) + FV1(m,c%hs)*sgn
      enddo
   enddo
   do col=1,N
      a(N,col) = (-1._R8P)**col
   enddo
   c%minv = inverse(a)
   contains
      subroutine mirror(idx, col, sgn)
      integer,   intent(in)  :: idx
      integer,   intent(out) :: col
      real(R8P), intent(out) :: sgn
      sgn = 1._R8P ; col = idx
      if (idx < 1) then ; col = 1-idx ; sgn = -1._R8P ; endif
      if (idx > N) then ; col = 2*N+1-idx ; sgn = -1._R8P ; endif
      endsubroutine mirror
   endsubroutine set_londrillo

   function inverse(a) result(b)
   real(R8P), intent(in) :: a(:,:)
   real(R8P)             :: b(size(a,1),size(a,1)), w(size(a,1),2*size(a,1)), t(2*size(a,1))
   integer               :: n0, i, k, p
   n0 = size(a,1) ; w = 0._R8P ; w(:,1:n0) = a
   do i=1,n0 ; w(i,n0+i) = 1._R8P ; enddo
   do k=1,n0
      p = k-1+maxloc(abs(w(k:n0,k)),dim=1)
      t = w(k,:) ; w(k,:) = w(p,:) ; w(p,:) = t
      w(k,:) = w(k,:)/w(k,k)
      do i=1,n0
         if (i /= k) w(i,:) = w(i,:) - w(i,k)*w(k,:)
      enddo
   enddo
   b = w(:,n0+1:)
   endfunction inverse

   ! ---------------------------------------------------------------------------------------------------- 1D line operators
   function filt_p(c, a) result(o)
   !< filter_modified_axis: o_i = sum_s p(s) a_{i+s} (zero extension).
   type(config), intent(in) :: c
   real(R8P),    intent(in) :: a(N)
   real(R8P)                :: o(N)
   integer :: i, s
   o = 0._R8P
   do i=1,N ; do s=-c%radius,c%radius
      if (i+s >= 1 .and. i+s <= N) o(i) = o(i) + c%p(s)*a(i+s)
   enddo ; enddo
   endfunction filt_p

   function filt_pt(c, a) result(o)
   !< Transpose of filt_p: o_m = sum_s p(s) a_{m-s}.
   type(config), intent(in) :: c
   real(R8P),    intent(in) :: a(N)
   real(R8P)                :: o(N)
   integer :: m, s
   o = 0._R8P
   do m=1,N ; do s=-c%radius,c%radius
      if (m-s >= 1 .and. m-s <= N) o(m) = o(m) + c%p(s)*a(m-s)
   enddo ; enddo
   endfunction filt_pt

   function current_line(c, h, dt) result(j)
   !< build_esirkepov_face_line + modified_current_line (matched) or the londrillo solve, one line along d.
   type(config), intent(in) :: c
   real(R8P),    intent(in) :: h(N), dt
   real(R8P)                :: j(N), flux(0:N), rhs(N)
   integer :: i, t, f
   flux(0) = 0._R8P
   do i=1,N ; flux(i) = flux(i-1) - h(i)/dt ; enddo
   j = 0._R8P
   select case(trim(c%kind))
   case('matched')
      do i=1,N ; do t=0,c%nq-1
         f = i+c%qfirst+t
         if (f >= 0 .and. f <= N) j(i) = j(i) + c%qw(t)*flux(f)
      enddo ; enddo
   case('londrillo')
      rhs(1:N-1) = flux(1:N-1) ; rhs(N) = 0._R8P
      j = matmul(c%minv, rhs)
   endselect
   endfunction current_line

   function line_to_faces_t(c, g) result(gf)
   !< Transpose of the faces -> centers map of current_line: gf(f), f = 0..N (flux index, face x_f + 1/2).
   type(config), intent(in) :: c
   real(R8P),    intent(in) :: g(N)
   real(R8P)                :: gf(0:N), y(N)
   integer :: i, t, f
   gf = 0._R8P
   select case(trim(c%kind))
   case('matched')
      do i=1,N ; do t=0,c%nq-1
         f = i+c%qfirst+t
         if (f >= 0 .and. f <= N) gf(f) = gf(f) + c%qw(t)*g(i)
      enddo ; enddo
   case('londrillo')
      y = matmul(transpose(c%minv), g)
      gf(1:N-1) = y(1:N-1)
   endselect
   endfunction line_to_faces_t

   ! --------------------------------------------------------------------------------------------------------- deposits
   subroutine shape_1d(c, x, w)
   !< Charge shape of one axis on the cell centers (before P for matched).
   type(config), intent(in)  :: c
   real(R8P),    intent(in)  :: x
   real(R8P),    intent(out) :: w(N)
   integer :: i
   do i=1,N
      if (trim(c%kind) == 'matched') then
         w(i) = bspl(c%order, x - c%shift - xc(i))
      else
         w(i) = eff(c%order, x - xc(i), c%filter)
      endif
   enddo
   endsubroutine shape_1d

   subroutine deposit_rho(c, xp, rho)
   type(config), intent(in)  :: c
   real(R8P),    intent(in)  :: xp(3)
   real(R8P),    intent(out) :: rho(N,N,N)
   real(R8P) :: w(N,3)
   integer   :: a, i, j, k
   do a=1,3
      call shape_1d(c, xp(a), w(:,a))
      if (trim(c%kind) == 'matched') w(:,a) = filt_p(c, w(:,a))
   enddo
   do k=1,N ; do j=1,N ; do i=1,N
      rho(i,j,k) = w(i,1)*w(j,2)*w(k,3)
   enddo ; enddo ; enddo
   endsubroutine deposit_rho

   subroutine deposit_j(c, xp, v, jc)
   !< Esirkepov over the symmetric displacement xp -/+ v*dt/2 (dt small), then the solver; poly: J = rho v.
   type(config), intent(in)  :: c
   real(R8P),    intent(in)  :: xp(3), v(3)
   real(R8P),    intent(out) :: jc(3,N,N,N)
   real(R8P), parameter :: dt = 1.e-4_R8P
   real(R8P) :: w0(N,3), w1(N,3), h(3,N,N,N), a0,a1,b0,b1,c0,c1,da,db,dc, rho(N,N,N), line(N)
   integer   :: a, i, j, k, d
   if (trim(c%kind) == 'poly') then
      call deposit_rho(c, xp, rho)
      do d=1,3 ; jc(d,:,:,:) = v(d)*rho ; enddo
      return
   endif
   do a=1,3
      call shape_1d(c, xp(a)-0.5_R8P*v(a)*dt, w0(:,a))
      call shape_1d(c, xp(a)+0.5_R8P*v(a)*dt, w1(:,a))
   enddo
   do k=1,N ; do j=1,N ; do i=1,N
      a0 = w0(i,1) ; a1 = w1(i,1) ; b0 = w0(j,2) ; b1 = w1(j,2) ; c0 = w0(k,3) ; c1 = w1(k,3)
      da = a1-a0 ; db = b1-b0 ; dc = c1-c0
      h(1,i,j,k) = da*(b0*c0 + 0.5_R8P*(db*c0+b0*dc) + db*dc/3._R8P)
      h(2,i,j,k) = db*(a0*c0 + 0.5_R8P*(da*c0+a0*dc) + da*dc/3._R8P)
      h(3,i,j,k) = dc*(a0*b0 + 0.5_R8P*(da*b0+a0*db) + da*db/3._R8P)
   enddo ; enddo ; enddo
   if (trim(c%kind) == 'matched') then ! P on the transverse axes
      do d=1,3
         do a=1,3
            if (a == d) cycle
            call apply_axis(h(d,:,:,:), a, .false.)
         enddo
      enddo
   endif
   do k=1,N ; do j=1,N
      line = current_line(c, h(1,:,j,k), dt) ; jc(1,:,j,k) = line
   enddo ; enddo
   do k=1,N ; do i=1,N
      line = current_line(c, h(2,i,:,k), dt) ; jc(2,i,:,k) = line
   enddo ; enddo
   do j=1,N ; do i=1,N
      line = current_line(c, h(3,i,j,:), dt) ; jc(3,i,j,:) = line
   enddo ; enddo
   contains
      subroutine apply_axis(f, ax, transposed)
      real(R8P), intent(inout) :: f(N,N,N)
      integer,   intent(in)    :: ax
      logical,   intent(in)    :: transposed
      integer :: p1, p2
      do p2=1,N ; do p1=1,N
         select case(ax)
         case(1) ; f(:,p1,p2) = filt_p(c, f(:,p1,p2))
         case(2) ; f(p1,:,p2) = filt_p(c, f(p1,:,p2))
         case(3) ; f(p1,p2,:) = filt_p(c, f(p1,p2,:))
         endselect
      enddo ; enddo
      if (transposed) continue
      endsubroutine apply_axis
   endsubroutine deposit_j

   ! ---------------------------------------------------------------------------------------------------------- gathers
   function gather_current(c, xp, f) result(g)
   !< bspline_field_weighting: order = field_weighting_model, filter = filter_deposition.
   type(config), intent(in) :: c
   real(R8P),    intent(in) :: xp(3), f(N,N,N)
   real(R8P)                :: g, w(N,3)
   integer :: a, i, j, k
   do a=1,3 ; do i=1,N
      w(i,a) = eff(c%gather_order, xp(a)-xc(i), c%filter)
   enddo ; enddo
   g = 0._R8P
   do k=1,N ; do j=1,N ; do i=1,N
      g = g + w(i,1)*w(j,2)*w(k,3)*f(i,j,k)
   enddo ; enddo ; enddo
   endfunction gather_current

   function gather_e(c, xp, f) result(g)
   !< Coherent E gather: transpose of the charge operator, then the charge shape.
   type(config), intent(in) :: c
   real(R8P),    intent(in) :: xp(3), f(N,N,N)
   real(R8P)                :: g, w(N,3), ft(N,N,N)
   integer :: a
   ft = f
   if (trim(c%kind) == 'matched') then
      do a=1,3 ; call axis_op(c, ft, a) ; enddo
   endif
   do a=1,3 ; call shape_1d(c, xp(a), w(:,a)) ; enddo
   g = contract(w(:,1), w(:,2), w(:,3), ft)
   endfunction gather_e

   function gather_j(c, d, xp, f) result(g)
   !< Coherent gather of a field multiplying j_d: transpose of the J_d operator, then the instantaneous Esirkepov shape
   !< (order-1 B-spline on the faces along d, the charge shape on the transverse axes).
   type(config), intent(in) :: c
   integer,      intent(in) :: d
   real(R8P),    intent(in) :: xp(3), f(N,N,N)
   real(R8P)                :: g, w(N,3), ft(N,N,N), wf(0:N), gf(0:N,N,N)
   integer :: a, i, j, k, p1, p2
   if (trim(c%kind) == 'poly') then
      g = gather_e(c, xp, f) ; return
   endif
   ft = f
   if (trim(c%kind) == 'matched') then
      do a=1,3
         if (a /= d) call axis_op(c, ft, a)
      enddo
   endif
   do a=1,3 ; call shape_1d(c, xp(a), w(:,a)) ; enddo
   do i=0,N
      if (trim(c%kind) == 'matched') then
         wf(i) = bspl(c%order-1, xp(d) - c%shift - (xc(i)+0.5_R8P))
      else
         wf(i) = eff(c%order-1, xp(d) - (xc(i)+0.5_R8P), c%filter)
      endif
   enddo
   ! lines along d to faces, the faces indexed first
   do p2=1,N ; do p1=1,N
      select case(d)
      case(1) ; gf(:,p1,p2) = line_to_faces_t(c, ft(:,p1,p2))
      case(2) ; gf(:,p1,p2) = line_to_faces_t(c, ft(p1,:,p2))
      case(3) ; gf(:,p1,p2) = line_to_faces_t(c, ft(p1,p2,:))
      endselect
   enddo ; enddo
   g = 0._R8P
   select case(d)
   case(1)
      do k=1,N ; do j=1,N ; do i=0,N ; g = g + wf(i)*w(j,2)*w(k,3)*gf(i,j,k) ; enddo ; enddo ; enddo
   case(2)
      do k=1,N ; do j=0,N ; do i=1,N ; g = g + w(i,1)*wf(j)*w(k,3)*gf(j,i,k) ; enddo ; enddo ; enddo
   case(3)
      do k=0,N ; do j=1,N ; do i=1,N ; g = g + w(i,1)*w(j,2)*wf(k)*gf(k,i,j) ; enddo ; enddo ; enddo
   endselect
   endfunction gather_j

   subroutine axis_op(c, f, ax)
   type(config), intent(in)    :: c
   real(R8P),    intent(inout) :: f(N,N,N)
   integer,      intent(in)    :: ax
   integer :: p1, p2
   do p2=1,N ; do p1=1,N
      select case(ax)
      case(1) ; f(:,p1,p2) = filt_pt(c, f(:,p1,p2))
      case(2) ; f(p1,:,p2) = filt_pt(c, f(p1,:,p2))
      case(3) ; f(p1,p2,:) = filt_pt(c, f(p1,p2,:))
      endselect
   enddo ; enddo
   endsubroutine axis_op

   pure function contract(wx, wy, wz, f) result(g)
   real(R8P), intent(in) :: wx(N), wy(N), wz(N), f(N,N,N)
   real(R8P)             :: g
   integer :: i, j, k
   g = 0._R8P
   do k=1,N ; do j=1,N ; do i=1,N ; g = g + wx(i)*wy(j)*wz(k)*f(i,j,k) ; enddo ; enddo ; enddo
   endfunction contract

   ! ------------------------------------------------------------------------------------------- periodic field solves
   subroutine dft3(f, inverse_flag)
   complex(R8P), intent(inout) :: f(N,N,N)
   logical,      intent(in)    :: inverse_flag
   complex(R8P) :: w(N,N), t(N)
   integer :: a, b, p1, p2
   real(R8P) :: sgn
   sgn = -1._R8P ; if (inverse_flag) sgn = 1._R8P
   do a=1,N ; do b=1,N
      w(a,b) = exp(cmplx(0._R8P, sgn*2._R8P*PI*real((a-1)*(b-1),R8P)/real(N,R8P), R8P))
   enddo ; enddo
   do p2=1,N ; do p1=1,N ; t = f(:,p1,p2) ; f(:,p1,p2) = matmul(w,t) ; enddo ; enddo
   do p2=1,N ; do p1=1,N ; t = f(p1,:,p2) ; f(p1,:,p2) = matmul(w,t) ; enddo ; enddo
   do p2=1,N ; do p1=1,N ; t = f(p1,p2,:) ; f(p1,p2,:) = matmul(w,t) ; enddo ; enddo
   if (inverse_flag) f = f/real(N,R8P)**3
   endsubroutine dft3

   pure function symbol(hs, m) result(d)
   !< Centered FD first derivative of half stencil hs, real part of its symbol / i.
   integer, intent(in) :: hs, m
   real(R8P)           :: d, k
   integer :: s
   k = 2._R8P*PI*real(m-1,R8P)/real(N,R8P) ; d = 0._R8P
   do s=1,hs ; d = d + 2._R8P*FD1(s,hs)*sin(real(s,R8P)*k) ; enddo
   endfunction symbol

   subroutine solve_e(hs, rho, e)
   !< E = -grad_h phi, -div_h grad_h phi = rho (null modes of the wide Laplacian dropped).
   integer,   intent(in)  :: hs
   real(R8P), intent(in)  :: rho(N,N,N)
   real(R8P), intent(out) :: e(3,N,N,N)
   complex(R8P) :: r(N,N,N), t(N,N,N)
   real(R8P) :: d(3), l
   integer :: i, j, k, a
   r = cmplx(rho, 0._R8P, R8P) ; call dft3(r, .false.)
   do a=1,3
      do k=1,N ; do j=1,N ; do i=1,N
         d = [symbol(hs,i), symbol(hs,j), symbol(hs,k)] ; l = sum(d*d)
         t(i,j,k) = (0._R8P,0._R8P)
         if (l > 1.e-12_R8P) t(i,j,k) = -cmplx(0._R8P, d(a), R8P)*r(i,j,k)/l
      enddo ; enddo ; enddo
      call dft3(t, .true.) ; e(a,:,:,:) = real(t, R8P)
   enddo
   endsubroutine solve_e

   subroutine solve_b(hs, jc, b)
   !< B = curl_h K J, K = (-lap_h)^-1 (magnetostatics of the particle's own current).
   integer,   intent(in)  :: hs
   real(R8P), intent(in)  :: jc(3,N,N,N)
   real(R8P), intent(out) :: b(3,N,N,N)
   complex(R8P) :: jh(3,N,N,N), t(N,N,N)
   real(R8P) :: d(3), l
   integer :: i, j, k, a
   do a=1,3
      t = cmplx(jc(a,:,:,:), 0._R8P, R8P) ; call dft3(t, .false.) ; jh(a,:,:,:) = t
   enddo
   do a=1,3
      do k=1,N ; do j=1,N ; do i=1,N
         d = [symbol(hs,i), symbol(hs,j), symbol(hs,k)] ; l = sum(d*d)
         t(i,j,k) = (0._R8P,0._R8P)
         if (l > 1.e-12_R8P) then
            select case(a)
            case(1) ; t(i,j,k) = cmplx(0._R8P,1._R8P,R8P)*(d(2)*jh(3,i,j,k) - d(3)*jh(2,i,j,k))/l
            case(2) ; t(i,j,k) = cmplx(0._R8P,1._R8P,R8P)*(d(3)*jh(1,i,j,k) - d(1)*jh(3,i,j,k))/l
            case(3) ; t(i,j,k) = cmplx(0._R8P,1._R8P,R8P)*(d(1)*jh(2,i,j,k) - d(2)*jh(1,i,j,k))/l
            endselect
         endif
      enddo ; enddo ; enddo
      call dft3(t, .true.) ; b(a,:,:,:) = real(t, R8P)
   enddo
   endsubroutine solve_b
endmodule gc_kernels

program gather_consistency
use gc_kernels
implicit none
type(config), allocatable :: cases(:)
real(R8P) :: frand(N,N,N), rho(N,N,N), jc(3,N,N,N), e(3,N,N,N), b(3,N,N,N), xp(3), v(3)
real(R8P) :: adj_cur(4), adj_coh(4), fe_cur, fe_coh, fb_cur, fb_coh, fb_rho, s_dep, s_abs, g, f(3), bt(3,3)
real(R8P), parameter :: F_REF = 1._R8P/(4._R8P*PI)
integer :: ic, d, ip, a, npos

call random_seed(put=[(12345+ic, ic=1,64)])
call random_number(frand) ; frand = frand - 0.5_R8P
call build_cases(cases)
v = [0.8_R8P, 0.5_R8P, -0.3_R8P]
npos = 8
write(*,'(a)') '# adjoint errors (rho, J_x, J_y, J_z) and self-forces / F_ref, current gather | coherent gather'
write(*,'(a28,2(a4,4a10,2a10),a4,a10)') 'case', ' |', 'rho','Jx','Jy','Jz','F_E','F_B', ' |', 'rho','Jx','Jy','Jz','F_E','F_B', &
                                         ' |', 'F_B(rho)'
do ic=1,size(cases)
   associate(c=>cases(ic))
   xp = [16.37_R8P, 15.81_R8P, 16.13_R8P]
   ! adjointness against a random field
   call deposit_rho(c, xp, rho)
   s_dep = sum(rho*frand) ; s_abs = sum(abs(rho*frand))
   adj_cur(1) = abs(s_dep - gather_current(c, xp, frand))/s_abs
   adj_coh(1) = abs(s_dep - gather_e(c, xp, frand))/s_abs
   call deposit_j(c, xp, v, jc)
   do d=1,3
      s_dep = sum(jc(d,:,:,:)*frand) ; s_abs = sum(abs(jc(d,:,:,:)*frand))
      adj_cur(1+d) = abs(s_dep - v(d)*gather_current(c, xp, frand))/s_abs
      adj_coh(1+d) = abs(s_dep - v(d)*gather_j(c, d, xp, frand))/s_abs
   enddo
   ! self-forces over sub-cell positions along x (y, z off-center)
   fe_cur = 0._R8P ; fe_coh = 0._R8P ; fb_cur = 0._R8P ; fb_coh = 0._R8P ; fb_rho = 0._R8P
   do ip=0,npos-1
      xp = [16._R8P + real(ip,R8P)/real(npos,R8P) + 0.03_R8P, 15.71_R8P, 16.29_R8P]
      call deposit_rho(c, xp, rho) ; call solve_e(c%hs, rho, e)
      do a=1,3
         fe_cur = max(fe_cur, abs(gather_current(c, xp, e(a,:,:,:))))
         fe_coh = max(fe_coh, abs(gather_e(c, xp, e(a,:,:,:))))
      enddo
      call deposit_j(c, xp, v, jc) ; call solve_b(c%hs, jc, b)
      do a=1,3 ; bt(a,1) = gather_current(c, xp, b(a,:,:,:)) ; enddo
      f = [v(2)*bt(3,1)-v(3)*bt(2,1), v(3)*bt(1,1)-v(1)*bt(3,1), v(1)*bt(2,1)-v(2)*bt(1,1)]
      fb_cur = max(fb_cur, maxval(abs(f)))
      ! bt(a,d): B_a gathered with the J_d operator
      do d=1,3 ; do a=1,3
         if (a /= d) bt(a,d) = gather_j(c, d, xp, b(a,:,:,:))
      enddo ; enddo
      f = [v(2)*bt(3,2)-v(3)*bt(2,3), v(3)*bt(1,3)-v(1)*bt(3,1), v(1)*bt(2,1)-v(2)*bt(1,2)]
      fb_coh = max(fb_coh, maxval(abs(f)))
      ! B gathered with the charge operator (the coherent E gather)
      do a=1,3 ; bt(a,1) = gather_e(c, xp, b(a,:,:,:)) ; enddo
      f = [v(2)*bt(3,1)-v(3)*bt(2,1), v(3)*bt(1,1)-v(1)*bt(3,1), v(1)*bt(2,1)-v(2)*bt(1,1)]
      fb_rho = max(fb_rho, maxval(abs(f)))
   enddo
   write(*,'(a28,2(a4,4es10.1,2es10.1),a4,es10.1)') trim(c%name), ' |', adj_cur, fe_cur/F_REF, fb_cur/F_REF, &
                                                         ' |', adj_coh, fe_coh/F_REF, fb_coh/F_REF, ' |', fb_rho/F_REF
   endassociate
enddo

contains
   subroutine build_cases(cs)
   type(config), allocatable, intent(out) :: cs(:)
   type(config) :: c
   integer :: o, q, h
   allocate(cs(0))
   ! non-conserving polynomial deposits, gather of the same order, filter off/on
   do o=0,6
      c = config() ; c%kind = 'poly' ; c%order = o ; c%gather_order = o ; c%hs = 3
      write(c%name,'(a,i0)') 'poly order ', o ; cs = [cs, c]
      c%filter = .true. ; write(c%name,'(a,i0,a)') 'poly order ', o, ' +filter' ; cs = [cs, c]
   enddo
   ! esirkepov_Q<N> on TSC, FD6
   do q=5,10
      if (q == 7 .or. q == 9) cycle
      c = config() ; c%order = 2 ; c%gather_order = 2 ; c%hs = 3 ; call set_matched(c, q)
      write(c%name,'(a,i0)') 'esirkepov_Q', q ; cs = [cs, c]
   enddo
   c = config() ; c%order = 2 ; c%gather_order = 2 ; c%hs = 3 ; call set_matched(c, 10) ; c%filter = .true.
   c%name = 'esirkepov_Q10 +filter flag' ; cs = [cs, c]
   ! totorica on TSC, FD2..FD8, and on cubic FD6
   do h=1,4
      c = config() ; c%order = 2 ; c%gather_order = 2 ; c%hs = h ; call set_matched(c, 0)
      write(c%name,'(a,i0)') 'totorica TSC FD', 2*h ; cs = [cs, c]
   enddo
   c = config() ; c%order = 3 ; c%gather_order = 3 ; c%hs = 3 ; call set_matched(c, 0)
   c%name = 'totorica cubic FD6' ; cs = [cs, c]
   ! londrillo on TSC FD6 (filter off/on) and on CIC FD4
   c = config() ; c%order = 2 ; c%gather_order = 2 ; c%hs = 3 ; call set_londrillo(c)
   c%name = 'londrillo TSC FD6' ; cs = [cs, c]
   c%filter = .true. ; c%name = 'londrillo TSC FD6 +filter' ; cs = [cs, c]
   c = config() ; c%order = 1 ; c%gather_order = 1 ; c%hs = 2 ; call set_londrillo(c)
   c%name = 'londrillo CIC FD4' ; cs = [cs, c]
   endsubroutine build_cases
endprogram gather_consistency
