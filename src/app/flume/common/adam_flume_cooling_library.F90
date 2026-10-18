!< ADAM, FLUME optically thin radiative cooling, Townsend's exact integration (issue #93, M11).
module adam_flume_cooling_library
!< ADAM, FLUME optically thin radiative cooling, Townsend's exact integration (issue #93, M11).
!<
!< Every routine is `pure`, takes explicit-size dummies and is tagged `!$acc routine seq` + `!$omp declare target`, as
!< in `adam_flume_euler_library`: the CPU loops and the FNL device kernels call the same source.
!<
!< The cooling law is `d(rho e)/dt = -n^2 Lambda(T)` at fixed density and momentum, `n = rho / m` (m the particle mass)
!< and `T = p / (rho R)`. `Lambda` is the piecewise power law through the nodes `(T_k, Lambda_k)`, `k = 1..N`:
!< `Lambda = Lambda_k (T / T_k)^alpha_k` on `[T_k, T_k+1]`, `alpha_k = ln(Lambda_k+1 / Lambda_k) / ln(T_k+1 / T_k)`;
!< above `T_N` the last segment continues; below `T_1` there is no cooling (the floor).
!<
!< Isochoric, `dT/dt = -(gamma - 1) n^2 Lambda(T) / (rho R)`. Townsend (2009, ApJS 181, 391, section 2) writes it
!< through the temporal evolution function
!<```
!< Y(T) = (Lambda_N / T_N) int_T^T_N dT' / Lambda(T'),
!< Y(T) = Y_k + (Lambda_N / Lambda_k) (T_k / T_N) [1 - (T / T_k)^(1 - alpha_k)] / (1 - alpha_k)   on segment k,
!<```
!< `Y_N = 0`, decreasing in `T` (the `alpha_k = 1` limit is the logarithm). Along the cooling `dY/dt` is the constant
!< `(gamma - 1) n^2 Lambda_N / (rho R T_N)`, so one step is exact for any `dt`: `Y_new = Y(T) + dt dY/dt`, then `T_new
!< = Y^-1(Y_new)` on the segment that holds `Y_new`, and `T_new = T_1` exactly when `Y_new >= Y_1`.
!<
!< Near `alpha_k = 1` both forms lose digits (`(1 - x^eps) / eps` with `eps` small); they are evaluated through
!< `expm1` and `log1p`, written here as series where the argument is small, so every `alpha_k` keeps full precision.

! third party modules
use :: penf, only : I4P, R8P

implicit none
private
public :: cooling_function
public :: cooling_evolution
public :: cooling_segment_y
public :: cooled_internal_energy
public :: cooled_temperature

real(R8P), parameter :: SERIES_THRESHOLD=1.e-2_R8P !< Below it, expm1 and log1p are summed as series (above, 1e-14 relative).

contains
   ! public procedures
   pure function cooling_function(n_nodes, t_nodes, l_nodes, alpha, temperature) result(lambda)
   !< Return `Lambda(T)`, zero below the first node.
   integer(I4P), intent(in) :: n_nodes          !< Nodes number, N >= 2.
   real(R8P),    intent(in) :: t_nodes(n_nodes) !< Node temperatures, increasing.
   real(R8P),    intent(in) :: l_nodes(n_nodes) !< Node values of Lambda.
   real(R8P),    intent(in) :: alpha(n_nodes)   !< Segment exponents, alpha(k) on [T_k, T_k+1] (alpha(N) unused).
   real(R8P),    intent(in) :: temperature      !< Temperature.
   real(R8P)                :: lambda           !< Lambda(T).
   integer(I4P)             :: k                !< Segment.
   !$acc routine seq
   !$omp declare target

   lambda = 0._R8P
   if (.not.(temperature > t_nodes(1))) return
   k = segment_of_temperature(n_nodes=n_nodes, t_nodes=t_nodes, temperature=temperature)
   lambda = l_nodes(k) * (temperature / t_nodes(k))**alpha(k)
   endfunction cooling_function

   pure function cooling_evolution(n_nodes, t_nodes, l_nodes, alpha, y_nodes, temperature) result(y)
   !< Return the temporal evolution function `Y(T)`, `T > T_1`.
   integer(I4P), intent(in) :: n_nodes          !< Nodes number, N >= 2.
   real(R8P),    intent(in) :: t_nodes(n_nodes) !< Node temperatures, increasing.
   real(R8P),    intent(in) :: l_nodes(n_nodes) !< Node values of Lambda.
   real(R8P),    intent(in) :: alpha(n_nodes)   !< Segment exponents.
   real(R8P),    intent(in) :: y_nodes(n_nodes) !< Y at the nodes, Y_N = 0.
   real(R8P),    intent(in) :: temperature      !< Temperature.
   real(R8P)                :: y                !< Y(T).
   integer(I4P)             :: k                !< Segment.
   !$acc routine seq
   !$omp declare target

   k = segment_of_temperature(n_nodes=n_nodes, t_nodes=t_nodes, temperature=temperature)
   y = y_nodes(k) + cooling_segment_y(n_nodes=n_nodes, t_nodes=t_nodes, l_nodes=l_nodes, alpha=alpha, k=k, &
                                      temperature=temperature)
   endfunction cooling_evolution

   pure function cooling_segment_y(n_nodes, t_nodes, l_nodes, alpha, k, temperature) result(y)
   !< Return `Y(T) - Y_k` along the power law of segment `k` (also beyond its end: the node values `Y_k` are built from
   !< it, and above `T_N` the last segment continues).
   integer(I4P), intent(in) :: n_nodes          !< Nodes number, N >= 2.
   real(R8P),    intent(in) :: t_nodes(n_nodes) !< Node temperatures, increasing.
   real(R8P),    intent(in) :: l_nodes(n_nodes) !< Node values of Lambda.
   real(R8P),    intent(in) :: alpha(n_nodes)   !< Segment exponents.
   integer(I4P), intent(in) :: k                !< Segment.
   real(R8P),    intent(in) :: temperature      !< Temperature.
   real(R8P)                :: y                !< Y(T) - Y_k.
   real(R8P)                :: c                !< Segment coefficient, (Lambda_N / Lambda_k) (T_k / T_N).
   !$acc routine seq
   !$omp declare target

   c = (l_nodes(n_nodes) / l_nodes(k)) * (t_nodes(k) / t_nodes(n_nodes))
   ! [1 - (T/T_k)^(1-alpha)] / (1-alpha) = -expm1((1-alpha) ln(T/T_k)) / (1-alpha)
   y = -c * expm1_over(eps=1._R8P - alpha(k), x=log(temperature / t_nodes(k)))
   endfunction cooling_segment_y

   pure function cooled_temperature(n_nodes, t_nodes, l_nodes, alpha, y_nodes, temperature, dy) result(t_new)
   !< Return the temperature after an isochoric cooling step that raises `Y` by `dy >= 0`; `T <= T_1` is not cooled.
   integer(I4P), intent(in) :: n_nodes          !< Nodes number, N >= 2.
   real(R8P),    intent(in) :: t_nodes(n_nodes) !< Node temperatures, increasing.
   real(R8P),    intent(in) :: l_nodes(n_nodes) !< Node values of Lambda.
   real(R8P),    intent(in) :: alpha(n_nodes)   !< Segment exponents.
   real(R8P),    intent(in) :: y_nodes(n_nodes) !< Y at the nodes, Y_N = 0.
   real(R8P),    intent(in) :: temperature      !< Temperature before the step.
   real(R8P),    intent(in) :: dy               !< Increment of Y over the step.
   real(R8P)                :: t_new            !< Temperature after the step.
   real(R8P)                :: y_new            !< Y after the step.
   real(R8P)                :: c                !< Segment coefficient.
   integer(I4P)             :: k                !< Segment.
   !$acc routine seq
   !$omp declare target

   t_new = temperature
   if (.not.(temperature > t_nodes(1))) return
   y_new = cooling_evolution(n_nodes=n_nodes, t_nodes=t_nodes, l_nodes=l_nodes, alpha=alpha, y_nodes=y_nodes, &
                             temperature=temperature) + dy
   if (y_new >= y_nodes(1)) then
      t_new = t_nodes(1)
      return
   endif
   ! the segment that holds y_new: y_nodes(k) > y_new >= y_nodes(k+1), the last one above T_N
   k = n_nodes - 1_I4P
   do while (k > 1_I4P)
      if (y_new < y_nodes(k)) exit
      k = k - 1_I4P
   enddo
   c = (l_nodes(n_nodes) / l_nodes(k)) * (t_nodes(k) / t_nodes(n_nodes))
   ! (T/T_k)^(1-alpha) = 1 - (1-alpha) (y_new - Y_k) / c, so ln(T/T_k) = log1p(-(1-alpha) z) / (1-alpha), z = (y-Y_k)/c
   t_new = t_nodes(k) * exp(log1p_over(eps=1._R8P - alpha(k), z=-(y_new - y_nodes(k)) / c))
   t_new = min(max(t_new, t_nodes(1)), temperature) ! round-off guard: a cooling step never heats nor undershoots
   endfunction cooled_temperature

   pure function cooled_internal_energy(n_nodes, t_nodes, l_nodes, alpha, y_nodes, gamma, R, dy_rate, dt, rho, e) &
      result(e_new)
   !< Return the internal energy density of a cell after an isochoric cooling step `dt`; a cell at or below the floor,
   !< or with `e <= 0`, keeps its value bit for bit (no round trip through the temperature).
   integer(I4P), intent(in) :: n_nodes          !< Nodes number, N >= 2.
   real(R8P),    intent(in) :: t_nodes(n_nodes) !< Node temperatures, increasing.
   real(R8P),    intent(in) :: l_nodes(n_nodes) !< Node values of Lambda.
   real(R8P),    intent(in) :: alpha(n_nodes)   !< Segment exponents.
   real(R8P),    intent(in) :: y_nodes(n_nodes) !< Y at the nodes, Y_N = 0.
   real(R8P),    intent(in) :: gamma            !< Specific heats ratio.
   real(R8P),    intent(in) :: R                !< Gas constant.
   real(R8P),    intent(in) :: dy_rate          !< (gamma - 1) Lambda_N / (m^2 R T_N).
   real(R8P),    intent(in) :: dt               !< Time step.
   real(R8P),    intent(in) :: rho              !< Density.
   real(R8P),    intent(in) :: e                !< Internal energy density, rho e.
   real(R8P)                :: e_new            !< Internal energy density after the step.
   real(R8P)                :: temperature      !< Temperature before the step.
   real(R8P)                :: t_new            !< Temperature after the step.
   !$acc routine seq
   !$omp declare target

   e_new = e
   temperature = (gamma - 1._R8P) * e / (rho * R)
   if (.not.(temperature > t_nodes(1))) return
   t_new = cooled_temperature(n_nodes=n_nodes, t_nodes=t_nodes, l_nodes=l_nodes, alpha=alpha, y_nodes=y_nodes, &
                              temperature=temperature, dy=dt * dy_rate * rho)
   if (t_new < temperature) e_new = min(rho * R * t_new / (gamma - 1._R8P), e)
   endfunction cooled_internal_energy

   ! private procedures
   pure function segment_of_temperature(n_nodes, t_nodes, temperature) result(k)
   !< Return the segment of a temperature `T > T_1`: `T_k <= T < T_k+1`, `N - 1` from `T_N - 1` up.
   integer(I4P), intent(in) :: n_nodes          !< Nodes number.
   real(R8P),    intent(in) :: t_nodes(n_nodes) !< Node temperatures, increasing.
   real(R8P),    intent(in) :: temperature      !< Temperature.
   integer(I4P)             :: k                !< Segment.
   !$acc routine seq
   !$omp declare target

   k = 1_I4P
   do while (k < n_nodes - 1_I4P)
      if (temperature < t_nodes(k+1)) exit
      k = k + 1_I4P
   enddo
   endfunction segment_of_temperature

   pure function expm1_over(eps, x) result(f)
   !< Return `(exp(eps x) - 1) / eps`, `x` if `eps = 0`, to full precision for any `eps`.
   real(R8P), intent(in) :: eps !< Exponent factor.
   real(R8P), intent(in) :: x   !< Argument.
   real(R8P)             :: f   !< Result.
   real(R8P)             :: a   !< eps x.
   !$acc routine seq
   !$omp declare target

   a = eps * x
   if (abs(a) < SERIES_THRESHOLD) then
      ! x (1 + a/2 + a^2/6 + ... + a^6/5040): the next term, a^7/40320 x, is below 3e-19 x
      f = x * (1._R8P + a * (1._R8P / 2._R8P + a * (1._R8P / 6._R8P + a * (1._R8P / 24._R8P + a * (1._R8P / 120._R8P &
              + a * (1._R8P / 720._R8P + a / 5040._R8P))))))
   else
      f = (exp(a) - 1._R8P) / eps
   endif
   endfunction expm1_over

   pure function log1p_over(eps, z) result(f)
   !< Return `log(1 + eps z) / eps`, `z` if `eps = 0`, to full precision for any `eps` (`1 + eps z > 0`).
   real(R8P), intent(in) :: eps !< Factor.
   real(R8P), intent(in) :: z   !< Argument.
   real(R8P)             :: f   !< Result.
   real(R8P)             :: a   !< eps z.
   !$acc routine seq
   !$omp declare target

   a = eps * z
   if (abs(a) < SERIES_THRESHOLD) then
      ! z (1 - a/2 + a^2/3 - ... - a^7/8): the next term, a^8/9 z, is below 2e-17 z
      f = z * (1._R8P - a * (1._R8P / 2._R8P - a * (1._R8P / 3._R8P - a * (1._R8P / 4._R8P - a * (1._R8P / 5._R8P &
              - a * (1._R8P / 6._R8P - a * (1._R8P / 7._R8P - a / 8._R8P)))))))
   else
      f = log(1._R8P + a) / eps
   endif
   endfunction log1p_over
endmodule adam_flume_cooling_library
