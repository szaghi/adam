!< ADAM, FLUME barotropic equation-of-state tables (issue #100, M13 P3).
module adam_flume_eos_tables_library
!< ADAM, FLUME barotropic equation-of-state tables (issue #100, M13 P3).
!<
!< A barotropic gas is a continuous piecewise polytrope, `P = K_k rho^G_k` on `rho_{k-1} <= rho < rho_k` (`rho_0 = 0`,
!< the last segment unbounded), `K_1 = c0^2 rho_ref^(1 - G_1)` and `K_{k+1} = K_k rho_k^(G_k - G_{k+1})` (P continuous),
!< `c^2 = dP/drho = G_k P / rho`: the law of src/tests/flume/verification/eos/eos_oracle.py (class Eos). c0 is the
!< isothermal speed sqrt(P / rho) at rho_ref, the sound speed there when G_1 = 1 (Mathew & Federrath 2021).
!<
!< The kernels take the EOS through one scalar, the argument named `gamma` (D-M13-1: one interface for every Euler
!< model): a table cannot travel that way, so each realm registers its table here and passes the table's id (as a
!< real). The registry is fixed-size, so the device holds it as module arrays (no allocatable descriptors in device
!< routines). Registering an identical table returns its id: realms with one law share one entry, a realm with another
!< law gets its own. A table is configuration, set once at initialization and never changed.
!< Every lookup is `pure` + `!$acc routine seq` + `!$omp declare target`, like the pointwise libraries.

! third party modules
use :: penf, only : I4P, R8P

implicit none
private
public :: EOS_SEGMENTS_MAX
public :: EOS_TABLES_MAX
public :: barotropic_c2
public :: barotropic_pressure
public :: eos_tables_k, eos_tables_breaks, eos_tables_exponents, eos_tables_segments
public :: register_barotropic_table

integer(I4P), parameter :: EOS_SEGMENTS_MAX = 8_I4P !< Segments of one table (fixed: device arrays).
integer(I4P), parameter :: EOS_TABLES_MAX = 4_I4P   !< Tables of one run (one per distinct realm law).

integer(I4P) :: eos_tables_number = 0_I4P                                          !< Registered tables.
integer(I4P) :: eos_tables_segments(EOS_TABLES_MAX) = 0_I4P                        !< Segments of each table.
real(R8P)    :: eos_tables_breaks(EOS_SEGMENTS_MAX, EOS_TABLES_MAX) = huge(1._R8P) !< Upper break of each segment.
real(R8P)    :: eos_tables_exponents(EOS_SEGMENTS_MAX, EOS_TABLES_MAX) = 1._R8P    !< G_k.
real(R8P)    :: eos_tables_k(EOS_SEGMENTS_MAX, EOS_TABLES_MAX) = 0._R8P            !< K_k.
!$acc declare create(eos_tables_segments, eos_tables_breaks, eos_tables_exponents, eos_tables_k)
!$omp declare target(eos_tables_segments, eos_tables_breaks, eos_tables_exponents, eos_tables_k)

contains
   ! public procedures
   function register_barotropic_table(c0, rho_ref, breaks, exponents) result(id)
   !< Register the table of `c0`, `rho_ref`, the breaks and the exponents (validated by the caller) and return its id;
   !< an identical table already registered returns that one's id. More than EOS_TABLES_MAX distinct tables: id 0
   !< (the caller stops). The device copy is the caller's (the FNL backend updates the arrays after initialization).
   real(R8P), intent(in) :: c0           !< Isothermal speed sqrt(P / rho) at rho_ref.
   real(R8P), intent(in) :: rho_ref      !< Reference density of the first segment.
   real(R8P), intent(in) :: breaks(:)    !< Densities of the breaks, increasing.
   real(R8P), intent(in) :: exponents(:) !< Exponent of each segment, size(breaks) + 1.
   integer(I4P)          :: id           !< Table id, 1..EOS_TABLES_MAX, 0 if the registry is full.
   real(R8P)             :: b(EOS_SEGMENTS_MAX), g(EOS_SEGMENTS_MAX), k(EOS_SEGMENTS_MAX) !< The candidate table.
   integer(I4P)          :: n, s, t      !< Segments, counters.

   n = size(exponents, kind=I4P)
   b = huge(1._R8P) ; g = 1._R8P ; k = 0._R8P
   b(1:n-1) = breaks(1:n-1)
   g(1:n) = exponents(1:n)
   k(1) = c0**2 * rho_ref**(1._R8P - g(1))
   do s=2, n
      k(s) = k(s-1) * b(s-1)**(g(s-1) - g(s))
   enddo
   do t=1, eos_tables_number
      if (eos_tables_segments(t) == n .and. all(eos_tables_breaks(:,t) == b) .and. &
          all(eos_tables_exponents(:,t) == g) .and. all(eos_tables_k(:,t) == k)) then
         id = t
         return
      endif
   enddo
   if (eos_tables_number == EOS_TABLES_MAX) then
      id = 0_I4P
      return
   endif
   eos_tables_number = eos_tables_number + 1_I4P
   id = eos_tables_number
   eos_tables_segments(id) = n
   eos_tables_breaks(:,id) = b
   eos_tables_exponents(:,id) = g
   eos_tables_k(:,id) = k
   endfunction register_barotropic_table

   pure function barotropic_pressure(id, rho) result(p)
   !< The pressure P(rho) of table `id` (a break belongs to the segment above it).
   integer(I4P), intent(in) :: id  !< Table id.
   real(R8P),    intent(in) :: rho !< Density.
   real(R8P)                :: p   !< Pressure.
   integer(I4P)             :: s   !< Segment.
   !$acc routine seq
   !$omp declare target

   s = segment(id=id, rho=rho)
   p = eos_tables_k(s,id) * rho**eos_tables_exponents(s,id)
   endfunction barotropic_pressure

   pure function barotropic_c2(id, rho) result(c2)
   !< The squared sound speed dP/drho = G P / rho of table `id`.
   integer(I4P), intent(in) :: id  !< Table id.
   real(R8P),    intent(in) :: rho !< Density.
   real(R8P)                :: c2  !< Squared sound speed.
   integer(I4P)             :: s   !< Segment.
   !$acc routine seq
   !$omp declare target

   s = segment(id=id, rho=rho)
   c2 = eos_tables_exponents(s,id) * eos_tables_k(s,id) * rho**(eos_tables_exponents(s,id) - 1._R8P)
   endfunction barotropic_c2

   ! private procedures
   pure function segment(id, rho) result(s)
   !< The segment of `rho` in table `id`: the first whose upper break exceeds rho (the last one is unbounded).
   integer(I4P), intent(in) :: id  !< Table id.
   real(R8P),    intent(in) :: rho !< Density.
   integer(I4P)             :: s   !< Segment.
   !$acc routine seq
   !$omp declare target

   s = 1_I4P
   do while (s < eos_tables_segments(id))
      if (rho < eos_tables_breaks(s,id)) exit
      s = s + 1_I4P
   enddo
   endfunction segment
endmodule adam_flume_eos_tables_library
