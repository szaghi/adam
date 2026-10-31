!< Print the eddy viscosities of the subgrid-scale formulas for the velocity gradients read on stdin (issue #119, T1).
program test_sgs
!< Print the eddy viscosities of the subgrid-scale formulas for the velocity gradients read on stdin (issue #119, T1).
!<
!< Each input line: `c dl(3) delta g(1:3,1:3)` (`g` in column order); each output line: Vreman, WALE, Sigma. sgs_oracle.py
!< compares them with its own computation.

use :: adam_turbulence_library, only : sgs_sigma, sgs_vreman, sgs_wale
use :: penf,                    only : I4P, R8P

implicit none
real(R8P)    :: c, dl(3), delta, g(3,3), nu(3)
integer(I4P) :: stat

do
   read(*, *, iostat=stat) c, dl, delta, g
   if (stat /= 0_I4P) exit
   call sgs_vreman(c=c, dl=dl, delta=delta, g=g, nu=nu(1))
   call sgs_wale(c=c, dl=dl, delta=delta, g=g, nu=nu(2))
   call sgs_sigma(c=c, dl=dl, delta=delta, g=g, nu=nu(3))
   write(*, '(3(ES24.16E3,1X))') nu
enddo
endprogram test_sgs
