#!/usr/bin/env python3
"""Check the one-sided Lax-Friedrichs partial states of the positivity limiter on first-order Rusanov states.

Along the first-order Rusanov run of lf_proto.py (GLM, the blast), count at every step the faces whose one-sided
partial states q_L - 2 D lambda (F_LF - f(q_L)) or q_R + 2 D lambda (F_LF - f(q_R)) have p <= 0, while the full
first-order update of every cell stays admissible (lf_proto.py): if they exist, the limiter's sufficient condition
is too strict for MHD (the Lax-Friedrichs splitting property fails, Wu 2018), not the scheme.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from lf_proto import GAMMA, flux, pressure  # noqa: E402

n, b0, ch, cfl, model = 128, 100.0 / np.sqrt(4.0 * np.pi), 60.0, 0.25, "glm"
h = 1.0 / n
x = (np.arange(n) + 0.5) * h
xx, yy = np.meshgrid(x, x, indexing="ij")
u = np.zeros((9, n, n))
u[0] = 1.0
u[5] = u[6] = b0 / np.sqrt(2.0)
u[4] = np.where((xx - 0.5) ** 2 + (yy - 0.5) ** 2 < 0.01, 1000.0, 0.1) / (GAMMA - 1.0) + 0.5 * (u[5] ** 2 + u[6] ** 2)
t, worst = 0.0, 0
for it in range(1, 616):
    fs, speeds = zip(*(flux(u, d, model, ch) for d in (0, 1)), strict=True)
    alpha = [np.maximum(s, ch) for s in speeds]
    dt = min(cfl / np.max(alpha[0] / h + alpha[1] / h), 0.01 - t)
    du = np.zeros_like(u)
    bad = 0
    for d in (0, 1):
        ur, fr, ar = np.roll(u, -1, axis=1 + d), np.roll(fs[d], -1, axis=1 + d), np.roll(alpha[d], -1, axis=d)
        a = np.maximum(alpha[d], ar)
        face = 0.5 * (fs[d] + fr) - 0.5 * a * (ur - u)
        step = 2 * 2 * dt / h
        vl = u - step * (face - fs[d])
        vr = ur + step * (face - fr)
        bad += int(np.sum(pressure(vl, model) <= 0.0) + np.sum(pressure(vr, model) <= 0.0))
        du -= (face - np.roll(face, 1, axis=1 + d)) / h
    u = u + dt * du
    t += dt
    worst = max(worst, bad)
    if it in (1, 2, 5, 10, 50, 100, 300, 615) or t >= 0.01 - 1e-15:
        print(f"step {it:4d} t {t:.3e}: inadmissible one-sided partial states {bad:5d}, "
              f"full-update min p {pressure(u, model).min():.3e}")
    if t >= 0.01 - 1e-15:
        break
print(f"max inadmissible partial states per step {worst}")
