#!/usr/bin/env python3
"""Generate PRISM input.ini for the Tesi_Vasile 'Final' geometry (Ansys Maxwell reference).

Geometry (conductor centre-lines, from Final.aedtresults):
  rectangular: 4 planar loops 40 (transverse) x 50 (axial) mm, planes x,y = +-30.5 mm, I = 100 A
               Y coils (By) phase 0, X coils (Bx) phase -90
               Ansys convention: By = -B0 cos(wt), Bx = +B0 sin(wt)
  saddle     : 120 deg arcs at z = +-25 mm, axial wires 50 mm
               X pair R = 25 mm, I = 43.95 A, phase -90 ; Y pair R = 28 mm, I = 50.32 A, phase 0
               Ansys convention: By = +B0 cos(wt), Bx = +B0 sin(wt)
sigma = 2 * dh (finest level).
"""
import argparse, math, re, pathlib

p = argparse.ArgumentParser()
p.add_argument('--template', required=True)
p.add_argument('--outdir', required=True)
p.add_argument('--case', choices=['rect', 'saddle'], required=True)
p.add_argument('--hx', type=float, required=True, help='half extent x=y [m]')
p.add_argument('--hz', type=float, required=True, help='half extent z [m]')
p.add_argument('--ni', type=int, required=True, help='cells per block along x=y')
p.add_argument('--nk', type=int, required=True, help='cells per block along z')
p.add_argument('--levels', type=int, required=True, help='max_level = iu_ref_levels')
p.add_argument('--freq', type=float, default=2.5e6)
p.add_argument('--time_delay', type=float, default=2.0e-8)
p.add_argument('--periods', type=float, default=1.0, help='periods simulated after the ramp')
p.add_argument('--it_max', type=int, default=0, help='override it_max (smoke tests)')
p.add_argument('--it_save', type=int, default=0)
p.add_argument('--cfl', type=float, default=0.499666666666666)
p.add_argument('--basename', required=True)
a = p.parse_args()

dx0 = 2 * a.hx / a.ni
dz0 = 2 * a.hz / a.nk
assert abs(dx0 - dz0) < 1e-12, f'non-cubic cells: dx0={dx0} dz0={dz0}'
dh = dx0 / 2 ** a.levels
sigma = 2 * dh
# empirical dt of PRISM FD6/RK-SSP54 runs: dt = 1.66670847797719e-12 s at dh = 2 mm, CFL 0.4997
dt = 1.66670847797719e-12 * (dh / 2e-3) * (a.cfl / 0.499666666666666)
t_end = a.time_delay + a.periods / a.freq
it_max = a.it_max or math.ceil(t_end / dt)
T_it = 1 / a.freq / dt
it_save = a.it_save or max(1, round(T_it / 8))   # 8 snapshots per period

txt = pathlib.Path(a.template).read_text()

def setkey(section, key, val):
    """Line-based replace of `key = value` inside [section], keeping any trailing ; comment."""
    global txt
    lines = txt.split('\n'); cur = None; done = False
    for n, ln in enumerate(lines):
        h = re.match(r'^\[(.+?)\]\s*$', ln)
        if h:
            cur = h.group(1); continue
        if cur == section and not done:
            m = re.match(r'^(' + re.escape(key) + r'\s*=\s*)([^;]*?)(\s*;.*)?$', ln)
            if m:
                lines[n] = m.group(1) + str(val) + (m.group(3) or ''); done = True
    assert done, (section, key)
    txt = '\n'.join(lines)

setkey('IO', 'output_basename', a.basename)
setkey('IO', 'restart_basename', a.basename + '-restart')
setkey('IO', 'it_save', it_save)
setkey('grms', 'radius', 0.015)   # Ansys Volume_controllo: cylinder D30 x L50
setkey('grms', 'length', 0.05)
setkey('time', 'time_max', f'{t_end * 1.01:.6e}')
setkey('time', 'it_max', it_max)
setkey('time', 'CFL', a.cfl)
for k, v in (('ni', a.ni), ('nj', a.ni), ('nk', a.nk),
             ('emin_x', -a.hx), ('emin_y', -a.hx), ('emin_z', -a.hz),
             ('emax_x', a.hx), ('emax_y', a.hx), ('emax_z', a.hz)):
    setkey('grid', k, v)
setkey('amr', 'max_level', a.levels)
setkey('amr', 'iu_ref_levels', a.levels)
for s in ('bc_x_min', 'bc_x_max', 'bc_y_min', 'bc_y_max', 'bc_z_min', 'bc_z_max'):
    setkey(s, 'type', 'PEC')

# ---- coils: drop old coil sections, write new ones
i0 = txt.index('[coils_input]')
i1 = txt.index('[external_fields]')
f = f'{a.freq:.10g}'
coils = ['; Coils input (geometry from Ansys Maxwell Final.aedtresults)', '[coils_input]',
         'circular_coils_number = 0']
if a.case == 'rect':
    coils += ['rectangular_coils_number = 4', 'solenoid_coils_number = 0', 'helicon_coils_number = 0',
              f'time_delay = {a.time_delay:.6e}']
    # (centre, normal, phase): Y coils phase 0, X coils phase -90 (Ansys Final)
    for n, (xc, yc, nrm, ph) in enumerate([(0, 0.0305, '+y', 0), (0, -0.0305, '+y', 0),
                                           (0.0305, 0, '+x', -90), (-0.0305, 0, '+x', -90)], 1):
        coils += [f'[coils_input_coil_{n}]', 'coil_type = rectangular', 'current_type = AC_current',
                  'lx = 0.040', 'ly = 0.050', f'sigma = {sigma:.6g}',
                  f'x_center = {xc}', f'y_center = {yc}', 'z_center = 0', f'normal = {nrm}',
                  'Amplitude = 100', f'Frequency = {f}', f'Phase = {ph}']
else:
    coils += ['rectangular_coils_number = 0', 'solenoid_coils_number = 0', 'helicon_coils_number = 4',
              f'time_delay = {a.time_delay:.6e}']
    L = 0.025
    def sad(R, phi_wire_minus_z, phi_other):
        # polyline: axial wire at phi_wire_minus_z from +L to -L (current -z, as Ansys terminal),
        # arc at -L to phi_other, axial wire up to +L, arc back
        c = lambda ph: (R * math.cos(math.radians(ph)), R * math.sin(math.radians(ph)))
        (x1, y1), (x2, y2) = c(phi_wire_minus_z), c(phi_other)
        return [(x1, y1, L), (x1, y1, -L), (x2, y2, -L), (x2, y2, L), (x1, y1, L)]
    # (R, phi of -z wire, phi of +z wire, amplitude, phase)  -> Ansys Final terminals
    spec = [(0.025, -60, 60, 43.95, -90),    # Saddle_X_pos
            (0.025, 240, 120, 43.95, -90),   # Saddle_X_neg
            (0.028, 30, 150, 50.32, 0),      # Saddle_Y_pos
            (0.028, 330, 210, 50.32, 0)]     # Saddle_Y_neg (terminal +z at 210 deg == -z at 330 deg)
    for n, (R, pm, pp, I, ph) in enumerate(spec, 1):
        pts = sad(R, pm, pp)
        coils += [f'[coils_input_coil_{n}]', 'coil_type = helicon', 'current_type = AC_current',
                  f'sigma = {sigma:.6g}', 'x_center = 0.000', 'y_center = 0.000', 'z_center = 0.000',
                  'normal = +z', f'r_coil = {R}', f'Amplitude = {I}', f'Frequency = {f}', f'Phase = {ph}',
                  f'N_points = {len(pts)}']
        for k, (x, y, z) in enumerate(pts, 1):
            coils += [f'x_point_{k} = {x:.13f}', f'y_point_{k} = {y:.13f}', f'z_point_{k} = {z:.13f}']
txt = txt[:i0] + '\n'.join(coils) + '\n\n' + txt[i1:]

out = pathlib.Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
(out / 'input.ini').write_text(txt)
nb = 8 ** a.levels
cells = nb * a.ni * a.ni * a.nk
print(f'{a.case}: dh={dh*1e3:.4g} mm sigma={sigma*1e3:.4g} mm grid={a.ni*2**a.levels}x{a.ni*2**a.levels}x{a.nk*2**a.levels}'
      f' cells={cells:.3e} blocks={nb} dt~{dt:.3e} s it/period={T_it:.0f} it_max={it_max} it_save={it_save}')
