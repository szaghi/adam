!< ADAM, global MPI handler singleton — single program-scope mpih_object instance.
module adam_mpih_global
!< ADAM, global MPI handler singleton — single program-scope mpih_object instance.
!<
!< `mpih` is an instance of [HERALD](https://github.com/szaghi/HERALD)'s `mpih_object` (issue #158): ranks and node
!< topology, the host memory budget of a rank, messages with the rank prefix, an `error_stop` that stops every rank and,
!< in the device builds, the binding of the rank to its device. It is ADAM's only program-scope singleton.

! external libraries
use :: herald, only : mpih_object
use :: penf,   only : I8P, R8P

implicit none
private
public :: host_memory_avail
public :: mpih

type(mpih_object), target :: mpih !< Program-scope MPI handler singleton.

contains
   function host_memory_avail() result(memory_avail)
   !< Return the host memory budget of this rank in the unit of `compute_blocks_number`: kB/1e6 (units of 1.024e9 bytes).
   !<
   !< The value of ADAM's former handler, bit for bit (issue #158, D1): HERALD's binding limit (bytes, the smaller of
   !< `MemTotal` and the cgroup limit) truncated to kB, divided by 1e6 and by the ranks of the node. HERALD's own
   !< `memory_avail` is in GiB (6.9 % smaller): the block capacity of the CPU path must not change with the wiring.
   !< Stops every rank if the budget is unknown (no `/proc`).
   real(R8P) :: memory_avail !< Host memory budget of this rank, kB/1e6.

   if (mpih%memory%limit <= 0_I8P) call mpih%error_stop(msg='host memory budget unknown: no /proc/meminfo nor cgroup limit')
   memory_avail = real(mpih%memory%limit/1024_I8P, R8P)/1e6/mpih%node_procs_number
   endfunction host_memory_avail
endmodule adam_mpih_global
