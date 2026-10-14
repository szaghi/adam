#!/usr/bin/env python3
"""Write the FLUME self-gravity Jeans cases of issue #81 (M8-P3): `sine-wave` with `wave_mode = jeans`.

A uniform gas, rho_0 = 1, p_0 = 3/5, gamma = 5/3 (sound speed c = 1), every active face periodic, with a Jeans mode of
wave vector k = 2 pi (1, 1) (`wave_angle` 45 degrees, `wavelength` 1/sqrt(2), so k^2 = 8 pi^2) under the self-gravity
of the FFT solver; the dispersion relation is w^2 = c^2 k^2 - 4 pi G rho_0:

    jeans2d   stable, G = pi/2 (4 pi G rho_0 = 2 pi^2, lambda = lambda_J / 2): w = pi sqrt(6), a right-running wave of
              phase speed w/k = sqrt(3)/2, amplitude 1e-6, to one period t = 2/sqrt(6) (SV-3, 2-D quadtree);
    jeans3d   jeans2d on the unit cube, an octree, z periodic (SV-3, 3-D: the wave in the x-y plane, z uniform);
    jeansu2d  unstable, G = 9 pi/2 (4 pi G rho_0 = 18 pi^2, lambda = 1.5 lambda_J): the growing mode, sigma = pi sqrt(10),
              amplitude 1e-6, to t = 0.5 (SV-4: the growth rate and the mirror symmetry);
    jeansc2d  jeansu2d at amplitude 1e-2 to t = 0.4, into the nonlinear collapse (the density contrast reaches about 0.5):
              SV-4, the drift of E + W and of the momentum against dt.

The `jeans` mode is `rho_0 (1 + A sin(2 pi (x + y)))`, even about the point (1/8, 1/8). The unstable cases shift the
box to [1/8, 9/8]^2, so that point is the box corner: the point reflection about it, an exact symmetry of the
equations (density even, velocity odd), maps cell (i, j) to (N + 1 - i, N + 1 - j) (`jeans_oracle.py symmetry`).

Usage:
    make_jeans.py <out.ini> <case> [--cells N] [--cfl C] [--time-max T] [--it-save N] [--amplitude A] [--levels L]

`--cells` sets the cells per block side (default 8); the 2-D cases have 4x4 blocks (two uniform levels), the 3-D one
2x2x2 (one level), so the grid is 4 N (2-D) or 2 N (3-D) cells per side.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

CASES = ("jeans2d", "jeans3d", "jeansu2d", "jeansc2d")
G_STABLE = math.pi / 2.0  # 4 pi G rho_0 = 2 pi^2: lambda = lambda_J / 2
G_UNSTABLE = 4.5 * math.pi  # 4 pi G rho_0 = 18 pi^2: lambda = 1.5 lambda_J
GAMMA = 5.0 / 3.0
P0 = 1.0 / GAMMA


def case_ini(case: str, cells: int, cfl: float, time_max: float | None, it_save: int, amplitude: float | None,
             levels: int | None) -> str:
    """Return the INI of a case."""
    is3d = case == "jeans3d"
    unstable = case in ("jeansu2d", "jeansc2d")
    lo = 0.125 if unstable else 0.0
    t_end = time_max if time_max is not None else {"jeansu2d": 0.5, "jeansc2d": 0.4}.get(case, 2.0 / math.sqrt(6.0))
    amp = amplitude if amplitude is not None else (1.0e-2 if case == "jeansc2d" else 1.0e-6)
    lev = levels if levels is not None else (1 if is3d else 2)
    lines = [
        f"; FLUME self-gravity case {case} (verification/selfgravity/make_jeans.py, issue #81 M8-P3)",
        "[IO]", f"output_basename = {case}", f"it_save = {it_save}", "restart = .false.",
        f"restart_basename = {case}-restart", "restart_save = 0", "residuals_save = 1", "save_memory_status = .false.",
        "save_residual_fields = .false.", "save_auxiliary_fields = .false.", "save_curl_fields = .false.",
        "save_divergence_fields = .false.", "save_gradient_fields = .false.", "save_laplacian_fields = .false.",
        "[grid]", f"ni = {cells}", f"nj = {cells}", f"nk = {cells if is3d else 1}", "ngc = 4",
        f"emin_x = {lo!r}", f"emin_y = {lo!r}", f"emin_z = {0.0!r}", f"emax_x = {lo + 1.0!r}", f"emax_y = {lo + 1.0!r}",
        f"emax_z = {1.0!r}",
        "null_x = .false.", "null_y = .false.", f"null_z = {'.false.' if is3d else '.true.'}",
        "[amr]", f"max_level = {lev}", f"ratio = {8 if is3d else 4}", f"iu_ref_levels = {lev}",
        "i_prune = 0", "j_prune = 0", "k_prune = 0", "l_prune = -1", "frequency = 0", "iters = 1", "markers_number = 0",
        "[field]", "nv = 5", "[runge_kutta]", "scheme = runge-kutta-ssp-33",
        "[weno]", "scheme = weno-u-5", "ror_number = 0", "ror_threshold = 0.9", "ror_vars_number = 0",
        "enable_ror_stats = .false.", "ib_reduction_extent = 0", "ib_reduced_order = 2",
        "[linear-algebra]", "smoothing = gauss-seidel", "iterations_init = 3", "iterations_coarse = 10",
        "iterations_fine = 3", "iterations = 10", "tolerance = 1.e-20",
        "[fdv]", "fdv_scheme = fd", "fdv_order = 2", "[solids]", "number = 0", "n_eikonal = 0",
        "[slices]", "slices_number = 0",
        "[physics]", "physical_model = euler", f"gamma = {GAMMA!r}",
        "[numerics]", "scheme_space = weno", "reconstruction_variables = characteristic", "reflux = .true.",
    ]
    for face in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max"):
        bc = "extrapolation" if face[0] == "z" and not is3d else "periodic"
        lines += [f"[bc_{face}]", f"type = {bc}"]
    lines += ["[time]", "it_max = -1", f"time_max = {t_end!r}", f"CFL = {cfl!r}",
              "[diagnostics]", "conservation_history_save = 1",
              "[gravity]", "self_gravity = .true.", f"G = {G_UNSTABLE if unstable else G_STABLE!r}",
              "[initial_conditions]", "type = sine-wave", "amr_iterations = 0", "wave_mode = jeans", "wave_angle = 45.0",
              f"wave_amplitude = {amp!r}", f"wavelength = {1.0 / math.sqrt(2.0)!r}",
              "[initial_conditions_region_1]", "r = 1.0", "u = 0.0", "v = 0.0", "w = 0.0", f"p = {P0!r}"]
    return "\n".join(lines) + "\n"


def main() -> int:
    """Parse the command line and write the input."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out", type=Path)
    parser.add_argument("case", choices=CASES)
    parser.add_argument("--cells", type=int, default=8)
    parser.add_argument("--cfl", type=float, default=0.4)
    parser.add_argument("--time-max", type=float, default=None)
    parser.add_argument("--it-save", type=int, default=1000000)
    parser.add_argument("--amplitude", type=float, default=None)
    parser.add_argument("--levels", type=int, default=None)
    args = parser.parse_args()
    args.out.write_text(case_ini(args.case, args.cells, args.cfl, args.time_max, args.it_save, args.amplitude,
                                 args.levels))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
