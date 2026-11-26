"""Apply the measurement-only phase-timer patch to a scratch ADAM tree (#133, 2026-11-24; updated for #156 L2, 2026-11-25: the ghost-exchange barriers are gone since b8641a59, and no FTZ switch, which nvfortran 26.3 does not resolve). NOT for commit."""

import re
import shutil
import sys
from pathlib import Path

T = Path(sys.argv[1])
shutil.copy(
    Path(__file__).with_name("adam_phase_timer.F90.txt"),
    T / "src/lib/common/adam_phase_timer.F90",
)


def region(s, start, end):
    a = s.index(start)
    b = s.index(end, a) + len(end)
    return a, b


def sub(s, start, end, old, new, count=1):
    a, b = region(s, start, end)
    r = s[a:b]
    n = r.count(old)
    assert n == count, f"{start[:50]!r}: {old[:60]!r} found {n} times, expected {count}"
    return s[:a] + r.replace(old, new) + s[b:]


def stmt_end(r, start):
    """index just after the Fortran statement beginning at `start` (follows & continuations)."""
    i = start
    while True:
        j = r.index("\n", i)
        line = r[i:j].split("!")[0].rstrip()
        if not line.endswith("&"):
            return j + 1
        i = j + 1


def wrap_stmt(s, start, end, head, ph, indent):
    a, b = region(s, start, end)
    r = s[a:b]
    assert r.count(head) == 1, f"{head[:50]!r} found {r.count(head)}"
    i = r.index(head)
    ls = r.rindex("\n", 0, i) + 1
    e = stmt_end(r, ls)
    r = (
        r[:ls]
        + f"{indent}call ph_tic({ph})\n"
        + r[ls:e]
        + f"{indent}call ph_toc({ph})\n"
        + r[e:]
    )
    return s[:a] + r + s[b:]


def add_use(s, after_regex):
    m = re.search(after_regex, s, re.M)
    assert m, after_regex
    e = s.index("\n", m.start()) + 1
    return s[:e] + "use :: adam_phase_timer\n" + s[e:]


# ---- forest
p = T / "src/lib/common/adam_forest_object.F90"
s = p.read_text()
s = add_use(s, r"^use :: adam_ini_library")
s = s.replace(
    "use :: adam_phase_timer\n",
    "use :: adam_phase_timer\n",
    1,
)
s = sub(
    s,
    "subroutine start_timing(self, realm)",
    "endsubroutine start_timing",
    "   call get_environment_variable('ADAM_SEAM_TIMING'",
    "   call ph_init()\n   call get_environment_variable('ADAM_SEAM_TIMING'",
)
s = sub(
    s,
    "subroutine compute_global_dt",
    "endsubroutine compute_global_dt",
    "   dt = huge(0._R8P)\n",
    "   dt = huge(0._R8P)\n   call ph_tic(PH_DT_LOCAL)\n",
)
s = sub(
    s,
    "subroutine compute_global_dt",
    "endsubroutine compute_global_dt",
    "   call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)\n",
    "   call ph_toc(PH_DT_LOCAL)\n   call ph_tic(PH_DT_ALLREDUCE)\n"
    "   call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)\n   call ph_toc(PH_DT_ALLREDUCE)\n",
)
loop_old = (
    "      if (self%timing) self%wtime(2) = self%wtime(2) - MPI_Wtime()\n"
    "      call self%evolve_one_step(realm=realm)\n"
    "      if (self%timing) self%wtime(2) = self%wtime(2) + MPI_Wtime()\n"
    "      call self%regrid(realm=realm)\n"
    "      call self%post_step(realm=realm)\n"
    "      call self%is_done(realm=realm, done=done)\n"
    "      if (done) exit\n   enddo\n"
)
loop_new = (
    "      call ph_tic(PH_STEP)\n"
    "      if (self%timing) self%wtime(2) = self%wtime(2) - MPI_Wtime()\n"
    "      call ph_tic(PH_EVOLVE)\n"
    "      call self%evolve_one_step(realm=realm)\n"
    "      call ph_toc(PH_EVOLVE)\n"
    "      if (self%timing) self%wtime(2) = self%wtime(2) + MPI_Wtime()\n"
    "      call ph_tic(PH_REGRID_DONE)\n      call self%regrid(realm=realm)\n      call ph_toc(PH_REGRID_DONE)\n"
    "      call ph_tic(PH_POST_STEP)\n      call self%post_step(realm=realm)\n      call ph_toc(PH_POST_STEP)\n"
    "      call ph_tic(PH_REGRID_DONE)\n      call self%is_done(realm=realm, done=done)\n      call ph_toc(PH_REGRID_DONE)\n"
    "      call ph_toc(PH_STEP)\n      call ph_end_step()\n"
    "      if (done) exit\n   enddo\n   call ph_report()\n"
)
n = s.count(loop_old)
assert n >= 1, "forest loop"
s = s.replace(loop_old, loop_new)
ftz = "   do\n      call ph_tic(PH_STEP)\n"  # no FTZ switch: nvfortran 26.3 does not resolve ieee_set_underflow_mode(gradual=)
m = s.count("   do\n      call ph_tic(PH_STEP)\n")
assert m == n, ("ftz anchor", m, n)
s = s.replace("   do\n      call ph_tic(PH_STEP)\n", ftz)
print("forest loops patched:", n)
p.write_text(s)

# ---- CPU field: update_ghost_mpi internals
p = T / "src/lib/common/adam_field_object.F90"
s = p.read_text()
s = add_use(s, r"^use :: adam_mpih_global")
R0, R1 = (
    "subroutine update_ghost_mpi(self, grid, maps, q, step, seam_guard)",
    "endsubroutine update_ghost_mpi",
)
s = sub(
    s,
    R0,
    R1,
    "   if (do_step(1)) then\n",
    "   call ph_tic(PH_MPI_PACK)\n   if (do_step(1)) then\n",
)
s = sub(
    s,
    R0,
    R1,
    "   if (do_step(2)) then\n",
    "   call ph_toc(PH_MPI_PACK)\n   call ph_tic(PH_MPI_POST)\n   if (do_step(2)) then\n",
)
s = sub(
    s,
    R0,
    R1,
    "   if (do_step(3)) then\n",
    "   call ph_toc(PH_MPI_POST)\n   if (do_step(3)) then\n",
)
s = sub(
    s,
    R0,
    R1,
    "      call MPI_WAITALL(procs_number * 2, req_send_recv, MPI_STATUSES_IGNORE, error)\n",
    "      call ph_tic(PH_MPI_WAITALL)\n"
    "      call MPI_WAITALL(procs_number * 2, req_send_recv, MPI_STATUSES_IGNORE, error)\n      call ph_toc(PH_MPI_WAITALL)\n"
    "      call ph_tic(PH_MPI_UNPACK)\n",
)
s = sub(
    s,
    R0,
    R1,
    "      endif\n   endif\n   endassociate\n",
    "      endif\n      call ph_toc(PH_MPI_UNPACK)\n   endif\n   endassociate\n",
)
p.write_text(s)

# ---- FNL field: update_ghost_mpi_gpu internals
p = T / "src/lib/fnl/adam_fnl_field_object.F90"
s = p.read_text()
s = add_use(s, r"^use :: adam_fnl_mpih_global")
R0, R1 = "subroutine update_ghost_mpi_gpu(", "endsubroutine update_ghost_mpi_gpu"
s = sub(
    s,
    R0,
    R1,
    "   if (do_step(1)) then\n",
    "   call ph_tic(PH_MPI_PACK)\n   if (do_step(1)) then\n",
)
s = sub(
    s,
    R0,
    R1,
    "   if (do_step(2)) then\n",
    "   call ph_toc(PH_MPI_PACK)\n   call ph_tic(PH_MPI_POST)\n   if (do_step(2)) then\n",
)
s = sub(
    s,
    R0,
    R1,
    "   if (do_step(3)) then\n",
    "   call ph_toc(PH_MPI_POST)\n   if (do_step(3)) then\n",
)
s = sub(
    s,
    R0,
    R1,
    "      call MPI_WAITALL(procs_number * 2, req_send_recv, MPI_STATUSES_IGNORE, error)\n",
    "      call ph_tic(PH_MPI_WAITALL)\n"
    "      call MPI_WAITALL(procs_number * 2, req_send_recv, MPI_STATUSES_IGNORE, error)\n      call ph_toc(PH_MPI_WAITALL)\n"
    "      call ph_tic(PH_MPI_UNPACK)\n",
)
s = sub(
    s,
    R0,
    R1,
    "   endif\n   endassociate\n",
    "      call ph_toc(PH_MPI_UNPACK)\n   endif\n   endassociate\n",
)
p.write_text(s)

# ---- FLUME CPU
p = T / "src/app/flume/cpu/adam_flume_cpu_object.F90"
s = p.read_text()
s = add_use(s, r"^use :: penf")
R0, R1 = "subroutine update_ghost(self, q)", "endsubroutine update_ghost\n"
s = sub(
    s,
    R0,
    R1,
    "   call self%adam%field%update_ghost_local(grid=self%adam%grid, maps=self%adam%maps, q=q)\n"
    "   call self%adam%field%update_ghost_mpi(grid=self%adam%grid, maps=self%adam%maps, q=q)\n"
    "   call self%set_boundary_conditions(q=q)\n",
    "   call ph_tic(PH_GHOST_LOCAL)\n   call self%adam%field%update_ghost_local(grid=self%adam%grid, maps=self%adam%maps, q=q)\n"
    "   call ph_toc(PH_GHOST_LOCAL)\n   call ph_tic(PH_GHOST_MPI)\n"
    "   call self%adam%field%update_ghost_mpi(grid=self%adam%grid, maps=self%adam%maps, q=q)\n   call ph_toc(PH_GHOST_MPI)\n"
    "   call ph_tic(PH_BC)\n   call self%set_boundary_conditions(q=q)\n",
)
s = sub(
    s,
    R0,
    R1,
    "   endsubroutine update_ghost\n",
    "   call ph_toc(PH_BC)\n   endsubroutine update_ghost\n",
)
R0, R1 = "subroutine compute_residuals_weno(", "endsubroutine compute_residuals_weno"
s = wrap_stmt(s, R0, R1, "call self%ib%evolve_eikonal(", "PH_IB", "         ")
s = wrap_stmt(s, R0, R1, "call self%ib%invert_eikonal(", "PH_IB", "      ")
s = sub(
    s,
    R0,
    R1,
    "   call self%compute_q_aux(q=q)\n   is_char = ",
    "   call ph_tic(PH_Q_AUX)\n   call self%compute_q_aux(q=q)\n   call ph_toc(PH_Q_AUX)\n   call ph_tic(PH_FACE_FLUXES)\n   is_char = ",
)
a, b = region(s, R0, R1)
r = s[a:b]
i = r.index("   call compute_flux_difference(nv=self%physics%nv,")
e = stmt_end(r, i)
r = (
    r[:i]
    + "   call ph_toc(PH_FACE_FLUXES)\n   call ph_tic(PH_FLUX_DIFF)\n"
    + r[i:e]
    + "   call ph_toc(PH_FLUX_DIFF)\n"
    + r[e:]
)
s = s[:a] + r + s[b:]
R0, R1 = "subroutine integrate_rk_ssp(self)", "endsubroutine integrate_rk_ssp"
for head, ind in [
    ("call self%ts%begin_step(", "   "),
    ("call self%ts%begin_stage(", "      "),
    ("call self%ts%end_stage(", "      "),
    ("call self%ts%finish_step(", "   "),
]:
    s = wrap_stmt(s, R0, R1, head, "PH_INTEGRATOR", ind)
s = wrap_stmt(s, R0, R1, "call self%save_residuals", "PH_RESIDUALS_SAVE", "   ")
s = wrap_stmt(
    s,
    "subroutine advance_one_step_forest(",
    "endsubroutine advance_one_step_forest",
    "call self%time%print_progress(",
    "PH_PROGRESS",
    "   ",
)
p.write_text(s)

# ---- FLUME FNL
p = T / "src/app/flume/fnl/adam_flume_fnl_object.F90"
s = p.read_text()
s = add_use(s, r"^use :: penf")
R0, R1 = "subroutine update_ghost(self, q_gpu)", "endsubroutine update_ghost\n"
s = wrap_stmt(
    s, R0, R1, "call self%field_fnl%update_ghost_local_gpu(", "PH_GHOST_LOCAL", "   "
)
s = wrap_stmt(
    s, R0, R1, "call self%field_fnl%update_ghost_mpi_gpu(", "PH_GHOST_MPI", "   "
)
s = sub(
    s,
    R0,
    R1,
    "   call self%set_boundary_conditions(q_gpu=q_gpu)\n",
    "   call ph_tic(PH_BC)\n   call self%set_boundary_conditions(q_gpu=q_gpu)\n",
)
s = sub(
    s,
    R0,
    R1,
    "   endsubroutine update_ghost\n",
    "   call ph_toc(PH_BC)\n   endsubroutine update_ghost\n",
)
R0, R1 = (
    "subroutine compute_residuals_weno_dev(",
    "endsubroutine compute_residuals_weno_dev",
)
s = wrap_stmt(s, R0, R1, "call self%ib_fnl%evolve_eikonal(", "PH_IB", "         ")
s = wrap_stmt(s, R0, R1, "call self%ib_fnl%invert_eikonal(", "PH_IB", "      ")
s = sub(
    s,
    R0,
    R1,
    "   call self%compute_q_aux(q_gpu=q_gpu)\n   is_char = ",
    "   call ph_tic(PH_Q_AUX)\n   call self%compute_q_aux(q_gpu=q_gpu)\n   call ph_toc(PH_Q_AUX)\n   call ph_tic(PH_FACE_FLUXES)\n   is_char = ",
)
old = "   if (self%ib%solids_number > 0_I4P) then\n      call compute_flux_difference_ib_dev("
s = sub(
    s,
    R0,
    R1,
    old,
    "   call ph_toc(PH_FACE_FLUXES)\n   call ph_tic(PH_FLUX_DIFF)\n" + old,
)
a, b = region(s, R0, R1)
r = s[a:b]
i = r.index("      call compute_flux_difference_dev(")
e = stmt_end(r, i)
assert r[e:].startswith("   endif\n"), r[e : e + 40]
r = r[:e] + "   endif\n   call ph_toc(PH_FLUX_DIFF)\n" + r[e + len("   endif\n") :]
s = s[:a] + r + s[b:]
R0, R1 = "subroutine integrate_rk_ssp_dev(self)", "endsubroutine integrate_rk_ssp_dev"
for head, ind in [
    ("call ts_begin_step(self)", "   "),
    ("call ts_begin_stage(self, s=s)", "      "),
    ("call self%ts%set_rhs(", "      "),
    ("call ts_end_stage(self, s=s)", "      "),
    ("call ts_finish_step(self)", "   "),
]:
    s = wrap_stmt(s, R0, R1, head, "PH_INTEGRATOR", ind)
s = wrap_stmt(s, R0, R1, "call self%save_residuals", "PH_RESIDUALS_SAVE", "   ")
s = wrap_stmt(
    s,
    "subroutine advance_one_step_forest(",
    "endsubroutine advance_one_step_forest",
    "call self%time%print_progress(",
    "PH_PROGRESS",
    "   ",
)
p.write_text(s)
print("phase-timer patch applied to", T)
