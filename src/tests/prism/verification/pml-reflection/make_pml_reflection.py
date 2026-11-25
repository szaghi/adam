#!/usr/bin/env python3
"""Write a PRISM input of the 1-D Gaussian pulse at normal incidence on an absorbing boundary (issue #145, item 3).

Why: the paper's absorbing-layer test (a pulse leaving a 1-D domain through a 20-cell classical PML backed by a
perfect conductor) checks one point. This writes the inputs of the sweep that measures the layer: the reflected
amplitude against the layer thickness N and the damping rate gamma_max, the CFS variant, and the boundaries without a
layer (Silver-Mueller, Neumann, PEC). The interior is the same in every input, [-0.08, 0.08] m at dx = 1 mm, so the
measurement window and the pulse path are the same; the layer of N cells is added outside it, on both x faces, with a
PEC wall behind it.

The pulse is the sm_comparison one of the co-author's inputs (src/tests/prism/fnl/fd/BC_tests_plane_wave/sm_comparison):
`gaussian_pulse` along +x, centre 0, sigma = 4 mm, B0 = 10 T (B_z, and E_y = c B_z), cut at 5 sigma; sixth-order
centred differences, SSP(5,4), CFL 0.5; y and z periodic, 6 cells each (a 1-D problem on a 3-D grid).

Usage:
    make_pml_reflection.py <out.ini> --boundary pml|cfs|silver-muller|neumann|pec|far [--cells N]
                           [--gamma-max G] [--exponent n] [--alfa-max A] [--k-max K] [--it-save S]

--cells N is the layer thickness in cells (pml, cfs); `far` puts a PEC wall 0.40 m beyond the interior, so no
reflection comes back within the measurement: the run measures the noise floor of the window (the dispersive trail of
the pulse). The run lasts until the pulse reflected by the wall behind the layer has crossed the window
(make_pml_reflection.t_end).
"""

from __future__ import annotations

import argparse
from pathlib import Path

C0 = 299792458.0
DX = 1.0e-3  # cell size [m]
HALF = 0.08  # half width of the interior [m]
WINDOW = 0.04  # measurement window [-WINDOW, WINDOW] [m]
SIGMA = 0.004  # pulse standard deviation [m]
FAR = 0.40  # distance of the wall of the `far` input beyond the interior [m]


def t_end(layer: float, kappa_max: float = 1.0) -> float:
    """Return the final time: the pulse, centred at 0, reflected by a wall at HALF + layer, has its centre at -WINDOW.
    A CFS layer slows the pulse by its stretching kappa = 1 + (kappa_max - 1) (d/L)^2: its round trip takes
    2 (kappa_max - 1) L / (3 c) more."""
    return (2.0 * (HALF + layer) + WINDOW + 2.0 * (kappa_max - 1.0) * layer / 3.0) / C0


def main() -> None:
    """Write the input."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out", type=Path)
    parser.add_argument(
        "--boundary",
        required=True,
        choices=("pml", "cfs", "silver-muller", "neumann", "pec", "far"),
    )
    parser.add_argument(
        "--cells", type=int, default=20, help="layer thickness in cells (pml, cfs)"
    )
    parser.add_argument("--gamma-max", default="7.8E10", help="peak damping rate [1/s]")
    parser.add_argument("--exponent", default="2.0", help="grading exponent n (pml)")
    parser.add_argument("--alfa-max", default="0.0", help="CFS frequency shift [1/s]")
    parser.add_argument("--k-max", default="1.0", help="CFS stretching kappa_max")
    parser.add_argument("--it-save", default="4", help="checkpoint period [steps]")
    args = parser.parse_args()
    layered = args.boundary in ("pml", "cfs")
    cells = args.cells if layered else 0
    extra = cells * DX if layered else FAR if args.boundary == "far" else 0.0
    ni = round(2.0 * (HALF + extra) / DX)
    kappa = float(args.k_max) if args.boundary == "cfs" else 1.0
    wall = {"silver-muller": "Silver_Muller", "neumann": "Neumann"}.get(
        args.boundary, "PEC"
    )
    if layered:
        common = (
            f"width          = {(cells - 0.25) * DX:.6E} ; {cells} cells (the count is ceiling(width/dx))\n"
            "x_minus_layer  = .true.\nx_plus_layer   = .true.\ny_minus_layer  = .false.\ny_plus_layer   = .false.\n"
            "z_minus_layer  = .false.\nz_plus_layer   = .false.\n"
            f"gamma_max      = {args.gamma_max}\n"
        )
        if args.boundary == "pml":
            pml = (
                f"PML_type       = classic\n{common}gamma_exponent = {args.exponent}\n"
            )
        else:
            pml = f"PML_type       = CFS\n{common}alfa_max       = {args.alfa_max}\nk_max          = {args.k_max}\n"
    else:
        pml = "PML_type       = none\n"
    text = f"""; 1-D Gaussian pulse at normal incidence on `{args.boundary}`{f", {cells} cells" if layered else ""} (issue #145 item 3,
; written by make_pml_reflection.py). Interior [-{HALF}, {HALF}] m, dx = {DX} m.
[numerics]
reconstruction_variables = conservative
scheme_time              = runge_kutta
scheme_space             = fd_centered
constrained_transport    = no
divergence_correction    = no

[fdv]
fdv_scheme = fd
fdv_order  = 6

[physics]
physical_model = EM

[runge_kutta]
scheme = runge-kutta-ssp-54

[IO]
output_basename        = pulse
it_save                = {args.it_save}
restart                = .false.
restart_basename       = pulse-restart
restart_save           = 10000000
residuals_save         = 10
save_memory_status     = .false.
save_residual_fields   = .false.
save_curl_fields       = .false.
save_divergence_fields = .false.

[time]
time_max = {t_end(extra if args.boundary != "far" else 0.0, kappa):.9E}
it_max   = -1
CFL      = 0.5

[grid]
ni     = {ni}
nj     = 6
nk     = 6
ngc    = 3
emin_x = {-(HALF + extra):.6f}
emin_y = -0.003
emin_z = -0.003
emax_x = {HALF + extra:.6f}
emax_y = 0.003
emax_z = 0.003
null_x = .false.
null_y = .false.
null_z = .false.

[initial_conditions]
amr_iterations  = 0
type            = gaussian_pulse
k0              = 0.0
r0              = 0.0
pulse_direction = +x
sigma           = {SIGMA}
B0              = 10.0

[amr]
max_level      = 0
ratio          = 8
iu_ref_levels  = 0
l_prune        = -1
frequency      = 50000
iters          = 1
markers_number = 0

[bc_x_min]
type = {wall}
[bc_x_max]
type = {wall}
[bc_y_min]
type = periodic
[bc_y_max]
type = periodic
[bc_z_min]
type = periodic
[bc_z_max]
type = periodic

[fWLayer]
width = 0.0
x_minus_layer = no
x_plus_layer  = no
y_minus_layer = no
y_plus_layer  = no
z_minus_layer = no
z_plus_layer  = no

[PML]
{pml}
[coils_input]
circular_coils_number = 0
rectangular_coils_number = 0
solenoid_coils_number = 0
helicon_coils_number = 0

[external_fields]
external_fields_applied = none

[solids]
number = 0

[slices]
slices_number = 0
"""
    args.out.write_text(text)


if __name__ == "__main__":
    main()
