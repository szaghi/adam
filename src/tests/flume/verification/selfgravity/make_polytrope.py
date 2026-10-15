#!/usr/bin/env python3
"""Write the FLUME self-gravity polytrope cases of issue #81 (M8-P4): `polytrope` with isolated boundaries.

An n = 1 polytrope star, radius a = 1/2, central density rho_c = 1, G = 1, centred in the cube [-1, 1]^3 (an octree,
uniform), under the isolated self-gravity of the FFT solver (`poisson_bc = isolated`, zero padding with the lattice
Green's function), inviscid walls at every face (a closed box: no mass or energy crosses it, so E + W has no
boundary term); gamma = 2, the adiabatic exponent of the n = 1 polytrope p = K rho^2,
K = 2 G a^2 / pi. The atmosphere has rho_1 = 1e-2 rho_c at the surface (`--atmosphere`), at the central temperature
(p_1 = K rho_c rho_1), isothermal and hydrostatic in the star's field outside it (the `polytrope` initial
conditions). It is not lighter on purpose: at 1e-4, the gas just outside the surface (a cell from a star density of
~0.1) is driven to |u| ~ 0.4 within ten steps and its pressure goes negative (32^3), which only the cell positivity
limiter (`--limiter`) survives; and that limiter's gravity energy correction (M7-G1, option B) breaks the exact
E + W exchange wherever it acts (2e-3 relative by t = 1, against 1e-9 without it). At 1e-2 the cases run without it.
The mean density of the star is 3 rho_c / pi^2, so its free-fall time is t_ff = sqrt(3 pi / (32 G rho_mean)) = 0.98.

    star   the equilibrium, `star_pressure_factor` 1, to t = 4 (SV-7: the spurious velocity, the drift of the central
           density, the drift of E + W);
    starc  the pressure halved (`star_pressure_factor` 0.5): the star collapses, bounces and rings, to t = 1.5 (SV-8:
           the drift of E + W and of the momentum against dt, while E and W exchange).

Usage:
    make_polytrope.py <out.ini> <case> [--cells N] [--levels L] [--cfl C] [--time-max T] [--it-save N]

`--cells` sets the cells per block side (default 8), `--levels` the uniform levels (default 2): the grid is
N 2^L cells per side (32^3 by default).
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

CASES = ("star", "starc")
G = 1.0
RADIUS = 0.5
RHO_C = 1.0
GAMMA = 2.0
RHO_ATM = 1.0e-2
K = 2.0 * G * RADIUS**2 / math.pi


def case_ini(case: str, cells: int, levels: int, cfl: float, time_max: float | None, it_save: int,
             atmosphere: float, limiter: bool) -> str:
    """Return the INI of a case."""
    t_end = time_max if time_max is not None else {"star": 4.0, "starc": 1.5}[case]
    factor = 1.0 if case == "star" else 0.5
    lines = [
        f"; FLUME self-gravity case {case} (verification/selfgravity/make_polytrope.py, issue #81 M8-P4)",
        "[IO]", f"output_basename = {case}", f"it_save = {it_save}", "restart = .false.",
        f"restart_basename = {case}-restart", "restart_save = 0", "residuals_save = 1", "save_memory_status = .false.",
        "save_residual_fields = .false.", "save_auxiliary_fields = .false.", "save_curl_fields = .false.",
        "save_divergence_fields = .false.", "save_gradient_fields = .false.", "save_laplacian_fields = .false.",
        "[grid]", f"ni = {cells}", f"nj = {cells}", f"nk = {cells}", "ngc = 4",
        "emin_x = -1.0", "emin_y = -1.0", "emin_z = -1.0", "emax_x = 1.0", "emax_y = 1.0", "emax_z = 1.0",
        "null_x = .false.", "null_y = .false.", "null_z = .false.",
        "[amr]", f"max_level = {levels}", "ratio = 8", f"iu_ref_levels = {levels}",
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
        f"positivity_limiter = {'cell' if limiter else 'none'}",
    ]
    for face in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max"):
        lines += [f"[bc_{face}]", "type = wall-inviscid"]
    lines += ["[time]", "it_max = -1", f"time_max = {t_end!r}", f"CFL = {cfl!r}",
              "[diagnostics]", "conservation_history_save = 1",
              "[gravity]", "self_gravity = .true.", f"G = {G!r}", "poisson_bc = isolated",
              "[initial_conditions]", "type = polytrope", "amr_iterations = 0", f"star_radius = {RADIUS!r}",
              f"star_density = {RHO_C!r}", f"star_pressure_factor = {factor!r}",
              "[initial_conditions_region_1]", f"r = {atmosphere * RHO_C!r}", "u = 0.0", "v = 0.0", "w = 0.0",
              f"p = {K * RHO_C * atmosphere * RHO_C!r}"]
    return "\n".join(lines) + "\n"


def main() -> int:
    """Parse the command line and write the input."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out", type=Path)
    parser.add_argument("case", choices=CASES)
    parser.add_argument("--cells", type=int, default=8)
    parser.add_argument("--levels", type=int, default=2)
    parser.add_argument("--cfl", type=float, default=0.4)
    parser.add_argument("--time-max", type=float, default=None)
    parser.add_argument("--it-save", type=int, default=1000000)
    parser.add_argument("--atmosphere", type=float, default=RHO_ATM, help="surface density of the atmosphere / rho_c")
    parser.add_argument("--limiter", action="store_true", help="the cell positivity limiter")
    args = parser.parse_args()
    args.out.write_text(case_ini(args.case, args.cells, args.levels, args.cfl, args.time_max, args.it_save,
                                 args.atmosphere, args.limiter))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
