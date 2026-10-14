#!/usr/bin/env python3
"""Write the FLUME gravity cases of issue #73 (M7): states under the constant gravity of `[gravity]`.

M7-G0 cases (the inputs): one negligible step (`CFL = 1e-30`) with the ghosts poisoned (`[diagnostics] ghost_poison`),
so the post-step checkpoint holds the initial atmosphere and its ghosts filled in the order of a Runge-Kutta stage;
`atmosphere_oracle.py profile` checks both:

    iso2d        isothermal, g = (0, -1), [0,1]^2 quadtree 4x4 blocks of 8^2, x periodic, y hydrostatic (Euler);
    poly2d       polytropic, G = 1.2, otherwise as iso2d;
    tilted2d     isothermal, g = (-0.6, -0.8), every face hydrostatic (the realm corners too), one refined block;
    tilted3d     isothermal, g = (-0.48, -0.64, -0.6), [0,1]^3 octree 4x4x4 blocks of 8^3, every face hydrostatic;
    iso2d-mhd    iso2d with EGLM and a field oblique to gravity, B = (0.5, 0.7, 0.2) (even at hydrostatic faces);
    froude2d     iso2d written with the number form, froude = 1, axis = -y.

M7-G1 cases (the force), a uniform state rho = 1, p = 1 at rest under g = (0, -1), run at CFL 0.4:

    freefall2d   every face periodic: the state falls freely, rho u = rho g t, E = E_0 + rho g^2 t^2 / 2 (50 steps);
    box2d        x periodic, y wall-inviscid: the gas falls onto the bottom wall and sloshes, a closed box whose mass and
                 total energy E + rho Phi (Phi = -g . x) are conserved (200 steps);
    box2d-lim    box2d with the positivity limiter and p = 0.05 (300 steps): the gas falls away from the top wall
                 toward vacuum, so the limiter acts (gravity in its factor, option B);
    box2d-stiff  box2d-lim with p = 1e-3: g dt is about 4.5 times the sound speed, the stiff regime;
    box2d-mhd    box2d with EGLM and a horizontal field B = (0.5, 0, 0) (the limiter on, as the EGLM regression cases).

M7-G2a case (the order of the energy source):

    wave2d       every face periodic, rho = 1 + 0.2 sin(2 pi (x + y)) at u_0 = (1, 0.5), p = 1, under g = (-0.6, -0.8):
                 the density wave falls freely (exact: advected by u_0 t + g t^2 / 2, u = u_0 + g t, p uniform), so the
                 errors at time_max give the order of the scheme with the gravity sources (`atmosphere_oracle.py wave`).

M7-G2b cases (the well-balanced mode), atmospheres at rest under a strong stratification, run at CFL 0.4 to
--time-max (default 0.5), x periodic, y hydrostatic:

    strat2d      isothermal, g = (0, -10): ten scale heights in the box (density ratio 2.2e4);
    polystrat2d  polytropic, G = 1.2, g = (0, -4) (density ratio 243);
    strat2d-mhdy strat2d with EGLM and a uniform field along g, B = (0, 0.05, 0) (beta 7 at the top);
    strat2d-mhdx strat2d with EGLM and a uniform field across g, B = (0.05, 0, 0);
    strat2d-amr  strat2d with a refined box [0.25, 0.75]^2 (2:1 seams across the stratification).

M7-G3 cases (Rayleigh-Taylor, the `rayleigh-taylor` initial condition), x periodic, y hydrostatic, g = (0, -1), two
isothermal layers 2 over 1 (rho_1 = 1, p_1 = 1 at the interface height 0.5), a solenoidal single mode of wavelength
0.5 seeded in the velocity, Gaussian envelope of width 0.1 about the interface, CFL 0.4:

    rt2d         interface width 0.05, seed 1e-5, to t = 4, on [0, 0.5] x [0.25, 0.75] (one wavelength, the mode
                 down to exp(-pi) at the y faces): the linear stage, sigma = 1.523 by rt_oracle.py (GV-4);
    rt2d-nl      interface width 0.01, seed 1e-2, to t = 5 (sigma = 1.869): the nonlinear single mode (bubbles,
                 spikes, roll-up);
    rt2d-mhd     rt2d-nl with EGLM and a horizontal field B = (0.1, 0, 0), half the critical field of the mode
                 (sqrt(drho g / (2 k)) = 0.20 for a sharp incompressible interface): tension slows the growth and
                 suppresses the small scales.

The MHD cases clean with EGLM, `glm_ch` = 3, `glm_alpha` = 0.18 and `glm_damping_length` = 1; `--divergence-control`,
`--glm-ch`, `--glm-alpha` and `--glm-damping-length` override them (issue #83: the mirror symmetry of rt2d-mhd under
each cleaning). `--balanced` sets [gravity] equilibrium = atmosphere; `--pulse A` adds the relative pressure pulse
p (1 + A exp(-|x - (0.5, 0.5, 0.5)|^2 / 0.05^2)) (the 2-D cells sit at z = 0.5).

Two uniform levels, so the blocks spread over the ranks (one level leaves every block on rank 0). Every case has four
ghost cells, the stencil of the fourth-order gravity energy source with WENO-5 (the face beyond each block end).

Usage:
    make_atmosphere.py <out.ini> <case> [--steps N] [--cfl C] [--cells N] [--time-max T] [--weights js|si]
                       [--riemann SOLVER] [--balanced] [--pulse A] [--it-save N] [--seed A]
                       [--divergence-control none|glm|eglm] [--glm-ch C] [--glm-alpha A] [--glm-damping-length L]

`--riemann SOLVER` selects the scheme `weno-riemann` with that solver and the sixth-order flux correction (the default
is the flux-splitting `weno`).
"""

from __future__ import annotations

import argparse
from pathlib import Path

G0_CASES = ("iso2d", "poly2d", "tilted2d", "tilted3d", "iso2d-mhd", "froude2d")
G1_CASES = ("freefall2d", "box2d", "box2d-lim", "box2d-stiff", "box2d-mhd")
G2_CASES = ("wave2d",)
G2B_CASES = ("strat2d", "polystrat2d", "strat2d-mhdy", "strat2d-mhdx", "strat2d-amr")
G3_CASES = ("rt2d", "rt2d-nl", "rt2d-mhd")
CASES = G0_CASES + G1_CASES + G2_CASES + G2B_CASES + G3_CASES
GRAVITY = {"iso2d": (0.0, -1.0, 0.0), "poly2d": (0.0, -1.0, 0.0), "tilted2d": (-0.6, -0.8, 0.0),
           "tilted3d": (-0.48, -0.64, -0.6), "iso2d-mhd": (0.0, -1.0, 0.0), "wave2d": (-0.6, -0.8, 0.0),
           "strat2d": (0.0, -10.0, 0.0), "polystrat2d": (0.0, -4.0, 0.0), "strat2d-mhdy": (0.0, -10.0, 0.0),
           "strat2d-mhdx": (0.0, -10.0, 0.0), "strat2d-amr": (0.0, -10.0, 0.0), "rt2d": (0.0, -1.0, 0.0),
           "rt2d-nl": (0.0, -1.0, 0.0), "rt2d-mhd": (0.0, -1.0, 0.0)}
G1_STEPS = {"freefall2d": 50, "box2d": 200, "box2d-lim": 300, "box2d-stiff": 200, "box2d-mhd": 200, "wave2d": 100}


def case_ini(case: str, steps: int | None, cfl: float | None, cells: int, time_max: float | None,
             weights: str | None, riemann: str | None = None, balanced: bool = False, pulse: float = 0.0,
             it_save: int | None = None, seed: float | None = None, cleaning: dict[str, str] | None = None) -> str:
    """Return the INI of a case; `cleaning` overrides the [mhd] divergence-control keys of the MHD cases."""
    g3 = case in G3_CASES
    g2b = case in G2B_CASES or g3
    if g2b and time_max is None and steps is None:
        time_max = {"rt2d": 4.0, "rt2d-nl": 5.0, "rt2d-mhd": 5.0}.get(case, 0.5)
    if g2b and cfl is None:
        cfl = 0.4
    g1 = case in G1_CASES or case in G2_CASES
    steps = steps if steps is not None else (G1_STEPS[case] if g1 else 1)
    cfl = cfl if cfl is not None else (0.4 if g1 else 1.0e-30)
    is3d = case == "tilted3d"
    mhd = case.endswith("-mhd") or "-mhd" in case
    tilted = case.startswith("tilted")
    clean = {"divergence_control": "eglm", "glm_ch": "3.0", "glm_alpha": "0.18", "glm_damping_length": "1.0"}
    clean |= cleaning or {}
    box = {"tilted2d": (0.0, 0.5, 0.0, 0.5, 0.0, 1.0), "strat2d-amr": (0.25, 0.75, 0.25, 0.75, 0.0, 1.0)}.get(case)
    limiter = case in ("box2d-lim", "box2d-stiff", "box2d-mhd")
    # rt2d: one wavelength wide, half a unit tall about the interface (the mode is down to exp(-pi) at the faces)
    lo_y, hi_x, hi_y = (0.25, 0.5, 0.75) if case == "rt2d" else (0.0, 1.0, 1.0)
    lines = [
        f"; FLUME gravity case {case} (verification/gravity/make_atmosphere.py, issue #73 M7)",
        "[IO]", f"output_basename = {case}", f"it_save = {it_save or (10 if g1 else 1000000)}", "restart = .false.",
        f"restart_basename = {case}-restart", "restart_save = 0", "residuals_save = 1", "save_memory_status = .false.",
        "save_residual_fields = .false.", "save_auxiliary_fields = .false.", "save_curl_fields = .false.",
        "save_divergence_fields = .false.", "save_gradient_fields = .false.", "save_laplacian_fields = .false.",
        "[grid]", f"ni = {cells}", f"nj = {cells}", f"nk = {cells if is3d else 1}", "ngc = 4",
        "emin_x = 0.0", f"emin_y = {lo_y!r}", "emin_z = 0.0", f"emax_x = {hi_x!r}", f"emax_y = {hi_y!r}",
        "emax_z = 1.0",
        "null_x = .false.", "null_y = .false.", f"null_z = {'.false.' if is3d else '.true.'}",
        "[amr]", f"max_level = {3 if box else 2}", f"ratio = {8 if is3d else 4}", "iu_ref_levels = 2",
        "i_prune = 0", "j_prune = 0", "k_prune = 0", "l_prune = -1", "frequency = 0", "iters = 1",
        f"markers_number = {1 if box else 0}",
    ]
    if box:
        lines += ["[amr_marker_1]", "mode = 1", "geo_type = primitive-box", "delta_type = max", "delta_fine = 0.0",
                  "delta_coarse = 0.0", f"box_xmin = {box[0]!r}", f"box_ymin = {box[2]!r}", f"box_zmin = {box[4]!r}",
                  f"box_xmax = {box[1]!r}", f"box_ymax = {box[3]!r}", f"box_zmax = {box[5]!r}", "target_level = 3"]
    lines += [
        "[field]", f"nv = {(8 if clean['divergence_control'] == 'none' else 9) if mhd else 5}", "[runge_kutta]",
        "scheme = runge-kutta-ssp-33",
        "[weno]", "scheme = weno-u-5", "ror_number = 0", "ror_threshold = 0.9", "ror_vars_number = 0",
        "enable_ror_stats = .false.", "ib_reduction_extent = 0", "ib_reduced_order = 2",
    ]
    if weights:
        lines += [f"weights = {weights}"]
    lines += [
        "[linear-algebra]", "smoothing = gauss-seidel", "iterations_init = 3", "iterations_coarse = 10",
        "iterations_fine = 3", "iterations = 10", "tolerance = 1.e-20",
        "[fdv]", "fdv_scheme = fd", "fdv_order = 2", "[solids]", "number = 0", "n_eikonal = 0",
        "[slices]", "slices_number = 0",
        "[physics]", f"physical_model = {'mhd-ideal' if mhd else 'euler'}", "gamma = 1.4",
        "[numerics]", f"scheme_space = {'weno-riemann' if riemann else 'weno'}",
        "reconstruction_variables = characteristic", "reflux = .true.",
    ]
    if riemann:
        lines += [f"riemann_solver = {riemann}", "flux_correction = 6th", "flux_correction_sensor = weno"]
    if limiter:
        lines += ["positivity_limiter = cell"]
    if mhd:
        lines += ["[mhd]", f"divergence_control = {clean['divergence_control']}", "divb_tol = 0.0",
                  "divb_error = .false.", "rho_floor = 0.0", "p_floor = 0.0", f"glm_ch = {clean['glm_ch']}",
                  f"glm_alpha = {clean['glm_alpha']}", f"glm_damping_length = {clean['glm_damping_length']}",
                  "glm_ch_check = warning"]
    for face in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max"):
        axis = face[0]
        if axis == "z" and not is3d:
            bc = "extrapolation"
        elif case in ("freefall2d", "wave2d") or (axis == "x" and not tilted):
            bc = "periodic"
        elif g1 and not g2b:
            bc = "wall-inviscid"
        else:
            bc = "hydrostatic"
        lines += [f"[bc_{face}]", f"type = {bc}"]
    lines += ["[time]", f"it_max = {steps if time_max is None else -1}",
              f"time_max = {1.0e30 if time_max is None else time_max!r}", f"CFL = {cfl!r}",
              "[diagnostics]", "ghost_poison = .true.", "conservation_history_save = 1"]
    if case == "froude2d":
        lines += ["[gravity]", "froude = 1.0", "axis = -y"]
    else:
        lines += ["[gravity]"] + [f"g_{x} = {v!r}" for x, v in zip("xyz", GRAVITY.get(case, (0.0, -1.0, 0.0)),
                                                                     strict=True) if v != 0.0]
    if balanced:
        lines += ["equilibrium = atmosphere"]
    if g3:
        lines += ["[initial_conditions]", "type = rayleigh-taylor", "amr_iterations = 1", "rt_height = 0.5",
                  f"rt_width = {0.05 if case == 'rt2d' else 0.01!r}", "rt_density_ratio = 2.0", "rt_wavelength = 0.5",
                  f"rt_amplitude = {seed if seed is not None else 1.0e-5 if case == 'rt2d' else 1.0e-2!r}",
                  "rt_seed_width = 0.1"]
    elif g2b:
        ic_type = "polytropic-atmosphere" if case == "polystrat2d" else "isothermal-atmosphere"
        lines += ["[initial_conditions]", f"type = {ic_type}", "amr_iterations = 1"]
        if case == "polystrat2d":
            lines += ["polytropic_gamma = 1.2"]
        if pulse:
            lines += [f"pulse_amplitude = {pulse!r}", "pulse_x = 0.5", "pulse_y = 0.5", "pulse_z = 0.5",
                      "pulse_width = 0.05"]
    elif case == "wave2d":
        lines += ["[initial_conditions]", "type = sine-wave", "amr_iterations = 1", "wave_mode = entropy",
                  "wave_angle = 45.0", "wave_amplitude = 0.2", "wavelength = 0.7071067811865476"]
    elif g1:
        lines += ["[initial_conditions]", "type = uniform", "amr_iterations = 1", "s = 0.0"]
    else:
        ic_type = "polytropic-atmosphere" if case == "poly2d" else "isothermal-atmosphere"
        lines += ["[initial_conditions]", f"type = {ic_type}", "amr_iterations = 1"]
        if case == "poly2d":
            lines += ["polytropic_gamma = 1.2"]
    state = {"r": 1.0, "u": 1.0 if case == "wave2d" else 0.0, "v": 0.5 if case == "wave2d" else 0.0, "w": 0.0, "p": {"box2d-lim": 0.05, "box2d-stiff": 1.0e-3}.get(case, 1.0)}
    if mhd:
        if case == "strat2d-mhdy":
            state |= {"bx": 0.0, "by": 0.05, "bz": 0.0}
        elif case == "rt2d-mhd":
            state |= {"bx": 0.1, "by": 0.0, "bz": 0.0}
        elif case == "strat2d-mhdx":
            state |= {"bx": 0.05, "by": 0.0, "bz": 0.0}
        elif g1 or g2b:
            state |= {"bx": 0.5, "by": 0.0, "bz": 0.0}
        else:
            state |= {"bx": 0.5, "by": 0.7, "bz": 0.2}
    lines += ["[initial_conditions_region_1]"] + [f"{k} = {v!r}" for k, v in state.items()]
    return "\n".join(lines) + "\n"


def main() -> int:
    """Parse the command line and write the case."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out", type=Path)
    parser.add_argument("case", choices=CASES)
    parser.add_argument("--steps", type=int, default=None)
    parser.add_argument("--cfl", type=float, default=None)
    parser.add_argument("--cells", type=int, default=8, help="cells per block side (the grid is 4 blocks a side)")
    parser.add_argument("--time-max", type=float, default=None)
    parser.add_argument("--weights", choices=("js", "si"), default=None)
    parser.add_argument("--riemann", choices=("llf", "hll", "hllc", "hlld"), default=None)
    parser.add_argument("--balanced", action="store_true", help="[gravity] equilibrium = atmosphere")
    parser.add_argument("--pulse", type=float, default=0.0, help="relative pressure pulse amplitude (G2b cases)")
    parser.add_argument("--it-save", type=int, default=None, help="checkpoint every N steps (default: the case's)")
    parser.add_argument("--seed", type=float, default=None, help="rt_amplitude of the G3 cases (0: the unseeded layers)")
    parser.add_argument("--divergence-control", choices=("none", "glm", "eglm"), default=None, help="MHD cases")
    parser.add_argument("--glm-ch", type=float, default=None, help="MHD cases: the cleaning speed")
    parser.add_argument("--glm-alpha", type=float, default=None, help="MHD cases: the damping coefficient")
    parser.add_argument("--glm-damping-length", default=None, help="MHD cases: a length or min-cell")
    args = parser.parse_args()
    cleaning = {k: str(v) for k, v in (("divergence_control", args.divergence_control), ("glm_ch", args.glm_ch),
                                       ("glm_alpha", args.glm_alpha), ("glm_damping_length", args.glm_damping_length))
                if v is not None}
    args.out.write_text(case_ini(args.case, args.steps, args.cfl, args.cells, args.time_max, args.weights,
                                 args.riemann, args.balanced, args.pulse, args.it_save, args.seed, cleaning))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
