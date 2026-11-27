"""#161 1.3 instrumentation (night scratch, never committed): per-rank timers of the CPU ghost exchange.
Usage: python3 -I instrument.py TREE [--barriers]   (--barriers restores b8641a59's two MPI_Barrier, timed)"""
import sys, pathlib
tree = pathlib.Path(sys.argv[1]); barriers = '--barriers' in sys.argv
f = tree / 'src/lib/common/adam_field_object.F90'; s = f.read_text()
def sub(old, new):
    global s
    assert s.count(old) == 1, old
    s = s.replace(old, new)
sub("public :: restrict_seam_high_order\n",
    "public :: restrict_seam_high_order\npublic :: ghost_mpi_wtime\n"
    "real(R8P), save :: ghost_mpi_wtime(6) = 0._R8P !< #161 1.3: pack, post, waitall, barriers, unpack, calls.\n")
sub("   if (do_step(1)) then\n      req_send_recv = MPI_REQUEST_NULL\n",
    "   if (do_step(1)) then\n      ghost_mpi_wtime(6) = ghost_mpi_wtime(6) + 1._R8P\n      ghost_mpi_wtime(1) = ghost_mpi_wtime(1) - MPI_Wtime()\n      req_send_recv = MPI_REQUEST_NULL\n")
sub("   if (do_step(2)) then\n      ! receive\n",
    "   if (do_step(1)) ghost_mpi_wtime(1) = ghost_mpi_wtime(1) + MPI_Wtime()\n   if (do_step(2)) then\n      ghost_mpi_wtime(2) = ghost_mpi_wtime(2) - MPI_Wtime()\n      ! receive\n")
sub("   if (do_step(3)) then\n      ! No barrier",
    "   if (do_step(2)) ghost_mpi_wtime(2) = ghost_mpi_wtime(2) + MPI_Wtime()\n   if (do_step(3)) then\n      ! No barrier")
bar = ("      ghost_mpi_wtime(4) = ghost_mpi_wtime(4) - MPI_Wtime()\n      call MPI_Barrier(MPI_COMM_WORLD, error)\n"
       "      ghost_mpi_wtime(4) = ghost_mpi_wtime(4) + MPI_Wtime()\n") if barriers else ""
sub("      call MPI_WAITALL(procs_number * 2, req_send_recv, MPI_STATUSES_IGNORE, error)\n      if (allocated(maps%comm_map_recv_ghost_cell)) then\n",
    "      ghost_mpi_wtime(3) = ghost_mpi_wtime(3) - MPI_Wtime()\n"
    "      call MPI_WAITALL(procs_number * 2, req_send_recv, MPI_STATUSES_IGNORE, error)\n"
    "      ghost_mpi_wtime(3) = ghost_mpi_wtime(3) + MPI_Wtime()\n" + bar +
    "      ghost_mpi_wtime(5) = ghost_mpi_wtime(5) - MPI_Wtime()\n"
    "      if (allocated(maps%comm_map_recv_ghost_cell)) then\n")
sub("            q(v_recv,i_recv,j_recv,k_recv,b_recv) = recv_buffer_ghost(c_send)\n         enddo\n      endif\n   endif\n   endassociate\n   endsubroutine update_ghost_mpi",
    "            q(v_recv,i_recv,j_recv,k_recv,b_recv) = recv_buffer_ghost(c_send)\n         enddo\n      endif\n"
    "      ghost_mpi_wtime(5) = ghost_mpi_wtime(5) + MPI_Wtime()\n" + bar.replace('      ', '      ', 3) +
    "   endif\n   endassociate\n   endsubroutine update_ghost_mpi")
f.write_text(s)
g = tree / 'src/lib/common/adam_forest_object.F90'; t = g.read_text()
def subg(old, new):
    global t
    assert t.count(old) == 1, old
    t = t.replace(old, new)
subg("use :: adam_globals,              only : mpih\n", "use :: adam_globals,              only : mpih\nuse :: adam_field_object,         only : ghost_mpi_wtime\n")
subg("   if (mpih%myrank == 0) call mpih%print_message(trim(line))\n   write(line, '(A,4(F11.4,A),I0,A)') 'forest regrid timing",
     "   if (mpih%myrank == 0) call mpih%print_message(trim(line))\n   call report_ghost_mpi(steps=self%wtime(2))\n   write(line, '(A,4(F11.4,A),I0,A)') 'forest regrid timing")
subg("   endsubroutine report_timing\n",
"""   endsubroutine report_timing

   subroutine report_ghost_mpi(steps)
   !< #161 1.3 instrumentation: min / mean / max over the ranks of this rank's steps time and ghost exchange timers.
   real(R8P), intent(in) :: steps     !< This rank's steps time.
   real(R8P)             :: v(7)      !< steps, pack, post, waitall, barriers, unpack, steps minus waitall and barriers.
   real(R8P)             :: lo(7), hi(7), su(7)
   integer               :: ierr, n
   character(len=512)    :: line
   v(1) = steps ; v(2:6) = ghost_mpi_wtime(1:5) ; v(7) = steps - ghost_mpi_wtime(3) - ghost_mpi_wtime(4)
   call MPI_ALLREDUCE(v, lo, 7, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(v, hi, 7, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(v, su, 7, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, n, ierr)
   write(line, '(A,I0,A,7(A,3(F10.3,A)))') 'ghost mpi timing (#161, min/mean/max over ranks, s; ', &
      nint(ghost_mpi_wtime(6)), ' calls on rank 0): ', &
      'steps ', lo(1), ' ', su(1)/n, ' ', hi(1), '; ', 'pack ', lo(2), ' ', su(2)/n, ' ', hi(2), '; ', &
      'post ', lo(3), ' ', su(3)/n, ' ', hi(3), '; ', 'waitall ', lo(4), ' ', su(4)/n, ' ', hi(4), '; ', &
      'barriers ', lo(5), ' ', su(5)/n, ' ', hi(5), '; ', 'unpack ', lo(6), ' ', su(6)/n, ' ', hi(6), '; ', &
      'busy ', lo(7), ' ', su(7)/n, ' ', hi(7), ''
   if (mpih%myrank == 0) call mpih%print_message(trim(line))
   endsubroutine report_ghost_mpi
""")
g.write_text(t)
print('instrumented', tree, 'barriers' if barriers else 'no barriers')
