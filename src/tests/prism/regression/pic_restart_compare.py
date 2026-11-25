"""Restart leg of pic_check.sh (issue #129): a PIC run stopped at it_max/2 and restarted, against the uninterrupted run.

Usage: pic_restart_compare.py <uninterrupted dir> <restarted dir> <basename> <half> <it_max> <rtol fields> <rtol particles> [--no-particles]

Compares, with the largest absolute difference relative to the largest magnitude of each pair (rtol 0: bit-identical;
the particles variable by variable, the first dt and the particles with <rtol particles>, the fields with <rtol fields>):
- the first step after the restart: dt and omega_max dt of rank 0, against step half+1 of the uninterrupted run (the
  plasma-frequency limit is printed also where the CFL sets dt, so a wrong restart of it shows);
- the final .pbd of every rank: header (width, particles, history), particles and plasma-frequency maxima;
- the fields of the it_max checkpoint, div_* left out (0 after a restart, #129 R5).
--no-particles (a case without particles, pml-restart): no .pbd compared.
Exit 1 on any difference above rtol.
"""
import glob
import os
import re
import sys

import h5py
import numpy as np

no_particles = '--no-particles' in sys.argv
sys.argv = [a for a in sys.argv if a != '--no-particles']
wd, rd, base = sys.argv[1], sys.argv[2], sys.argv[3]
half, it_max, rtol, prtol = int(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])
bad = []


def rel(a, b):
    """Largest |a - b| over the largest |a|, |b| (inf if the shapes differ)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.shape != b.shape:
        return np.inf
    s = max(np.abs(a).max(initial=0.), np.abs(b).max(initial=0.))
    return 0. if s == 0. else float(np.abs(a - b).max()/s)


def steps(fn):
    return [l for l in open(fn, errors='replace') if l.startswith('[mpi-00000]time step:')]


def numbers(line):
    return [float(x) for x in re.findall(r'[-+]?\d+\.\d+E[-+]\d+', line)]


# first step after the restart
u, r = steps(f'{wd}/run.log'), steps(f'{rd}/run2.log')
if len(u) < it_max or not r:
    bad.append('dt lines')
else:
    d = rel(numbers(u[half]), numbers(r[0]))
    print(f'>> first dt after the restart: {r[0].split(":", 1)[1].strip()} (relative difference {d:.1e})')
    if d > prtol:
        bad.append('first dt')

# final particles and plasma-frequency maxima, rank by rank
pbd = [] if no_particles else sorted(glob.glob(f'{wd}/{base}-restart-proc*.pbd'))
if not pbd and not no_particles:
    bad.append('no .pbd')
def particles(fn):
    """Particles [n, width] of a .pbd (q_pic is stored as [width, capacity], variable fastest) and the trailing maxima."""
    width, n, _ = np.fromfile(fn, np.int32, 3)
    raw = np.fromfile(fn, np.float64, offset=12)
    cap = (raw.size - 3)//width if (raw.size - 3) % width == 0 else raw.size//width
    return raw[:width*cap].reshape(cap, width)[:n], raw[width*cap:]


dmax, worst = 0., ''
for fu in pbd:
    fr = f'{rd}/{os.path.basename(fu)}'
    if not os.path.exists(fr) or not np.array_equal(np.fromfile(fu, np.int32, 3), np.fromfile(fr, np.int32, 3)):
        bad.append(f'{os.path.basename(fu)} header')
        continue
    (pu, mu), (pr, mr) = particles(fu), particles(fr)
    # every variable on its own scale (positions, velocities, charge, mass differ by many orders of magnitude)
    for v in range(pu.shape[1]):
        d = rel(pu[:, v], pr[:, v])
        if d > dmax:
            dmax, worst = d, f'variable {v + 1}'
    for v in range(max(mu.size, mr.size)):
        d = rel(mu[v:v + 1], mr[v:v + 1])
        if d > dmax:
            dmax, worst = d, f'plasma-frequency maximum {v + 1}'
if not no_particles:
    print(f'>> particles (.pbd, {len(pbd)} ranks, per variable): max relative difference {dmax:.1e} {worst}')
if dmax > prtol:
    bad.append('particles')

# fields of the last checkpoint
dmax, worst = 0., ''
h5 = sorted(glob.glob(f'{wd}/{base}-{it_max:09d}-proc*.h5'))
if not h5:
    bad.append(f'no it {it_max} checkpoint')
for fu in h5:
    fr = f'{rd}/{os.path.basename(fu)}'
    if not os.path.exists(fr):
        bad.append(f'{os.path.basename(fu)} missing')
        continue
    with h5py.File(fu, 'r') as a, h5py.File(fr, 'r') as b:
        for k in a.keys():
            if k.rsplit('-', 1)[1].startswith('div_'):
                continue
            d = rel(a[k][()], b[k][()]) if k in b else np.inf
            if d > dmax:
                dmax, worst = d, k
print(f'>> fields at it {it_max}: max relative difference {dmax:.1e} {worst}')
if dmax > rtol:
    bad.append('fields')

if bad:
    print('  restart mismatch: ' + ', '.join(bad))
sys.exit(1 if bad else 0)
