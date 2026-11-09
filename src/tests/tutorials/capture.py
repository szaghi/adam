#!/usr/bin/env python3
"""Run the tutorials and regenerate the program output the documentation shows.

Every block of program output in the tutorials (docs/tutorials/, docs/start/first-run.md, the home page's "In a
nutshell") is a file under docs/snippets/tutorials/, written by this script from a real run and included by the page
(`<<< @/snippets/tutorials/<name>.txt` for a text block, `<!--@include: @/snippets/tutorials/<name>.md-->` for a table).
The numbers the prose quotes are checked as claims: the page must contain the string formatted from the run.

    capture.py                 run every CPU case, rewrite the snippets, check the claims
    capture.py --check         run, compare with the committed snippets (numbers within --rtol/--atol), check the claims
    capture.py --only sod,lax  only these runs (and the snippets and claims that read them)
    capture.py --gpu           also the FNL cases (PRISM pulse, RIPPLE r4): their round-off digits are the machine's
    capture.py --gpu --only pulse-fnl,ripple-r4 --with-existing   the FNL cases alone, after a CPU capture

Exit status 0 when every run behaved as expected (two runs must stop: an error message is their output), every
snippet matched (--check) and every claim held; 1 otherwise, with the differences printed.

Executables come from --exe-dir (default <repo>/exe): adam_flume_cpu, adam_prism_cpu, adam_ripple_cpu,
test_ripple_crest_criterion, and with --gpu adam_prism_fnl, adam_ripple_fnl. The runs go to --work (default
<repo>/exe/tutorials-work); RIPPLE's check.sh writes into src/tests/ripple/work (git-ignored). The FNL runs on WSL need
the GPU environment of the cookbook (UCX_RNDV_THRESH=inf and the rest) exported by the caller.
"""

import argparse
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
DOCS = REPO / "docs"
SNIPPETS = DOCS / "snippets" / "tutorials"
PYTHON = os.environ.get("PYTHON", str(REPO / "exe" / ".regression-venv" / "bin" / "python"))
MPIRUN = os.environ.get("MPIRUN", "mpirun")

SOD = "src/tests/flume/regression/sod-x/input.ini"
PULSE = "src/tests/prism/regression/periodic-pulse-blocks/input.ini"
LAX = {"IO": {"output_basename": "lax"},
       "initial_conditions_region_1": {"r": "0.445", "u": "0.698", "p": "3.528"},
       "initial_conditions_region_2": {"r": "0.5", "u": "0.0", "p": "0.571"},
       "time": {"time_max": "0.13"}}
HLLC = {"numerics": {"scheme_space": "weno-riemann", "reconstruction_variables": "characteristic",
                     "riemann_solver": "hllc", "flux_correction": "6th", "flux_correction_sensor": "weno",
                     "reflux": ".true."}}
LIMITER = {"time": {"it_max": "-1", "time_max": "0.01"}}
FACES = {"initial_conditions": {"r0": "0.0"}}
PML = {"PML": {"PML_type": "classic", "width": "0.04", "x_minus_layer": ".true.", "x_plus_layer": ".true.",
               "y_minus_layer": ".false.", "y_plus_layer": ".false.", "z_minus_layer": ".false.",
               "z_plus_layer": ".false.", "gamma_max": "7.8E10", "gamma_exponent": "2.0"},
       "time": {"it_max": "-1", "time_max": "6.6712814328639916E-10"}}


def merge(*edits):
    out = {}
    for e in edits:
        for sec, keys in e.items():
            out.setdefault(sec, {}).update(keys)
    return out


def faces(kind):
    return merge(FACES, {"bc_x_min": {"type": kind}, "bc_x_max": {"type": kind}})


# name: (executable, input, edits, ranks, must stop with an error, needs a GPU)
RUNS = {
    "sod": ("adam_flume_cpu", SOD, {}, 2, False, False),
    "sod-u3": ("adam_flume_cpu", SOD, {"weno": {"scheme": "weno-u-3"}}, 2, False, False),
    "sod-u7": ("adam_flume_cpu", SOD, {"weno": {"scheme": "weno-u-7"}, "grid": {"ngc": "4"}}, 2, False, False),
    "sod-u7-ngc3": ("adam_flume_cpu", SOD, {"weno": {"scheme": "weno-u-7"}}, 2, True, False),
    "lax": ("adam_flume_cpu", SOD, LAX, 2, False, False),
    "lax-hllc": ("adam_flume_cpu", SOD, merge(LAX, HLLC), 2, False, False),
    "blast-amr": ("adam_flume_cpu", "src/tests/flume/regression/blast-amr-limiter/input.ini", {}, 2, False, False),
    "limiter-none": ("adam_flume_cpu", "src/tests/flume/regression/blast-limiter/input.ini",
                     merge(LIMITER, {"numerics": {"positivity_limiter": "none"}}), 2, True, False),
    "limiter-cell": ("adam_flume_cpu", "src/tests/flume/regression/blast-limiter/input.ini",
                     merge(LIMITER, {"numerics": {"positivity_limiter": "cell"}}), 2, False, False),
    "pulse": ("adam_prism_cpu", PULSE, {}, 2, False, False),
    "pulse-centred": ("adam_prism_cpu", PULSE, FACES, 2, False, False),
    "pulse-pec": ("adam_prism_cpu", PULSE, faces("PEC"), 2, False, False),
    "pulse-sm": ("adam_prism_cpu", PULSE, faces("Silver_Muller"), 2, False, False),
    "pulse-pml": ("adam_prism_cpu", PULSE, merge(faces("Silver_Muller"), PML), 2, False, False),
    "pulse-fnl": ("adam_prism_fnl", PULSE, {}, 2, False, True),
    "ripple-periodic": ("adam_ripple_cpu", "src/tests/ripple/periodic.ini", {}, 2, False, False),
    "ripple-r1": ("check.sh", "r1", {}, 0, False, False),
    "ripple-r2r3": ("check.sh", "r2,r3", {}, 0, False, False),
    "ripple-r6": ("check.sh", "r6", {}, 0, False, False),
    "ripple-r4": ("check.sh", "r4", {}, 0, False, True),
}


# ---------------------------------------------------------------------------------------------------- input and runs
def edit_ini(text, edits):
    """Set `key = value` in `[section]` (replacing the line, or adding it after the header; a missing section is
    appended). Comments and the other lines are kept."""
    lines = text.splitlines()
    for sec, keys in edits.items():
        head = next((i for i, l in enumerate(lines) if l.strip().lower() == f"[{sec.lower()}]"), None)
        if head is None:
            lines += ["", f"[{sec}]"] + [f"{k} = {v}" for k, v in keys.items()]
            continue
        end = next((i for i in range(head + 1, len(lines)) if lines[i].lstrip().startswith("[")), len(lines))
        for k, v in keys.items():
            pat = re.compile(rf"^\s*{re.escape(k)}\s*=", re.IGNORECASE)
            hit = next((i for i in range(head + 1, end) if pat.match(lines[i])), None)
            if hit is None:
                lines.insert(head + 1, f"{k} = {v}")
                end += 1
            else:
                lines[hit] = f"{k} = {v}"
    return "\n".join(lines) + "\n"


class Run:
    def __init__(self, name, work, exe_dir):
        self.name = name
        self.exe, self.input, self.edits, self.np, self.must_fail, self.gpu = RUNS[name]
        self.dir = work / name
        self.exe_dir = exe_dir
        self.status = None
        self.stdout = ""

    def execute(self):
        if self.exe == "check.sh":
            env = dict(os.environ, RIPPLE_CPU_EXE=str(self.exe_dir / "adam_ripple_cpu"),
                       RIPPLE_FNL_EXE=str(self.exe_dir / "adam_ripple_fnl"),
                       RIPPLE_CRIT_EXE=str(self.exe_dir / "test_ripple_crest_criterion"), PYTHON=PYTHON)
            p = subprocess.run([str(REPO / "src/tests/ripple/check.sh"), "--legs", self.input], cwd=REPO, env=env,
                               capture_output=True, text=True)
            self.status, self.stdout = p.returncode, p.stdout
            self.dir.mkdir(parents=True, exist_ok=True)
            (self.dir / "log.txt").write_text(p.stdout + p.stderr)
            return
        if self.dir.exists():
            shutil.rmtree(self.dir)
        self.dir.mkdir(parents=True)
        (self.dir / "input.ini").write_text(edit_ini((REPO / self.input).read_text(), self.edits))
        with open(self.dir / "log.txt", "w") as log:
            self.status = subprocess.run([MPIRUN, "-np", str(self.np), str(self.exe_dir / self.exe), "input.ini"],
                                         cwd=self.dir, stdout=log, stderr=subprocess.STDOUT).returncode

    def ok(self):
        return (self.status != 0) if self.must_fail else (self.status == 0)

    def log(self):
        return (self.dir / "log.txt").read_text(errors="replace").splitlines()

    def first(self, text):
        """The first log line holding `text`, rank 0's when it printed one (the ranks print in no fixed order)."""
        hits = [l for l in self.log() if text in l]
        return next((l for l in hits if l.startswith("[mpi-00000]")), hits[0])

    def file(self, suffix):
        hits = sorted(self.dir.glob(f"*{suffix}"))
        if not hits:
            raise FileNotFoundError(f"{self.name}: no *{suffix}")
        return hits[0]

    def rows(self, suffix):
        """Data rows (first field an integer) of a history file, as lists of tokens."""
        out = []
        for line in self.file(suffix).read_text().splitlines():
            tok = line.split()
            if tok and re.fullmatch(r"[+-]?\d+", tok[0]):
                out.append(tok)
        return out

    def header(self, suffix):
        return next(l for l in self.file(suffix).read_text().splitlines() if l.strip() and
                    not re.fullmatch(r"[+-]?\d+", l.split()[0]))

    def run_py(self, *args):
        p = subprocess.run([PYTHON, *map(str, args)], cwd=self.dir, capture_output=True, text=True)
        return p.stdout.splitlines()


# ----------------------------------------------------------------------------------------------------- formatting
def num(tok):
    return float(tok.replace("D", "E"))


def tex(x, digits=2):
    """1.2e-06 -> 1.2\\cdot10^{-6}"""
    m, e = f"{x:.{digits - 1}e}".split("e")
    return f"{m}\\cdot10^{{{int(e)}}}"


def oracle(run, *extra):
    """The oracle's verdict line (the run directory, printed as a path prefix, is dropped)."""
    line = run.run_py(REPO / "src/tests/flume/verification/sod/sod_oracle.py", "input.ini", ".", *extra)[-1]
    return line.replace(str(run.dir), "")


def l1(run):
    return num(re.search(r"L1\(rho\) = (\S+)", oracle(run)).group(1))


def compare(run, a, b, ref=None):
    args = [REPO / "src/tests/tutorials/compare_saves.py", ref or ".", "periodic_pulse_blocks", a, b]
    if ref:
        args = [REPO / "src/tests/tutorials/compare_saves.py", ref, "periodic_pulse_blocks", a, b, "."]
    return [l for l in run.run_py(*args) if re.match(r"(Bz|Dy)\s", l)]


def energy(run, steps):
    rows = {int(r[0]): r for r in run.rows("-energy_history.dat")}
    return [(rows[s], num(rows[s][3]) + num(rows[s][4])) for s in steps]


def energy_total(run, it):
    return energy(run, [it])[0][1]


def watts(w):
    if w >= 1:
        return f"{w:.2f}"
    if w >= 0.1:
        return f"{w:.3f}"
    return f"${tex(w)}$"


def check_lines(run, start, stop=None):
    """The lines of a check.sh output from the first matching `start` to `stop` included (the last line otherwise)."""
    lines = run.stdout.splitlines()
    i = next(k for k, l in enumerate(lines) if re.search(start, l))
    j = next((k for k in range(i, len(lines)) if stop and re.search(stop, lines[k])), len(lines) - 1)
    return lines[i:j + 1]


def xdmf_blocks(path):
    """Blocks and distinct 1/dx of a FLUME save: one `block_*` Grid per block in the XDMF index, its spacing in the
    `-dxdydz` dataset of the rank's HDF5 file."""
    import h5py
    text = path.read_text()
    inv = set()
    for f, d in re.findall(r"(\S+\.h5):(\S+-dxdydz)", text):
        with h5py.File(path.parent / f, "r") as h5:
            inv.add(round(1.0 / float(h5[d][-1])))          # ORIGIN_DXDYDZ is (dz, dy, dx)
    return len(re.findall(r'<Grid Name="block_', text)), sorted(inv)


# --------------------------------------------------------------------------------------------- snippets and claims
def snippets(R):
    """name -> (runs it reads, function returning the text)."""
    def sod_log():
        log = [l for l in R["sod"].log() if l.startswith("[mpi-00000]")]
        first = next(i for i, l in enumerate(log) if re.fullmatch(r"\[mpi-00000\]t:\s+1", l))
        last = max(i for i, l in enumerate(log) if re.fullmatch(r"\[mpi-00000\]t:\s+\d+", l))
        return ["…"] + log[first:first + 5] + ["…"] + log[last:]

    def sod_files():
        names = sorted(p.name for p in R["sod"].dir.iterdir() if p.name not in ("input.ini", "log.txt"))
        saves = sorted({re.match(r"(.*-\d{9})", n).group(1) for n in names if re.match(r".*-\d{9}", n)})
        lines = ["  ".join(sorted(n for n in names if n.startswith(s))[::-1][:1] +
                           sorted(n for n in names if n.startswith(s) and n.endswith(".h5"))) for s in saves]
        return lines + ["  ".join(n for n in names if n.endswith(".dat"))]

    def sod_history():
        rows = R["sod"].rows("-conservation_history.dat")
        return [f"{r[0]:>10}  " + "  ".join(r[1:4]) + "  ...  " + r[-1] for r in (rows[0], rows[-1])]

    def sod_orders():
        rows = [("`weno-u-3`", 3, "sod-u3"), ("`weno-u-5` (the input)", 3, "sod"), ("`weno-u-7`", 4, "sod-u7")]
        out = ["| `[weno] scheme` | `ngc` | $L_1(\\rho)$ against the exact solution, 200 cells, $t = 0.2$ |",
               "|---|---|---|"]
        return out + [f"| {s} | {g} | ${tex(l1(R[r]), 7)}$ |" for s, g, r in rows]

    def blast_history():
        run = R["blast-amr"]
        rows = run.rows("-conservation_history.dat")
        head = " ".join(run.header("-conservation_history.dat").split()[:3]) + " ..."
        w = max(len(rows[0][0]), len(rows[-1][0])) + 1
        return [head] + [f"{r[0]:<{w}} {r[1]} {r[2]}" for r in (rows[0], rows[-1])]

    def pulse_energy():
        rows = R["pulse"].rows("-energy_history.dat")
        w = max(len(rows[0][0]), len(rows[-1][0]))
        return [f"{r[0]:<{w}} " + " ".join(r[1:5]) + " ..." for r in (rows[0], rows[-1])]

    def table(run, steps, time_row=False):
        e = energy(R[run], steps)
        out = ["| step | " + " | ".join(map(str, steps)) + " |", "|---" * (len(steps) + 1) + "|"]
        if time_row:
            out.append("| $t$ (s) | " + " | ".join(f"${tex(num(r[2]), 3)}$" for r, _ in e) + " |")
        return out + ["| $W_D + W_B$ (J) | " + " | ".join(watts(w) for _, w in e) + " |"]

    def pml_log():
        return [R["pulse-pml"].first("C_gamma:"), R["pulse-pml"].first("time step:")]

    def ripple_periodic():
        lines = R["ripple-periodic"].file("-history.dat").read_text().splitlines()
        return [lines[0], lines[1], lines[-1]]

    def gpu_energy():
        cpu, fnl = R["pulse"].rows("-energy_history.dat")[-1], R["pulse-fnl"].rows("-energy_history.dat")[-1]
        return [f"{lab}  " + " ".join(r[:5]) + " ..." for lab, r in (("CPU", cpu), ("FNL", fnl))]

    def error_line(run):
        return [R[run].first("error stop")]

    return {
        "flume-sod-log.txt": (["sod"], sod_log),
        "flume-sod-end.txt": (["sod"], lambda: sod_log()[-2:]),
        "flume-sod-files.txt": (["sod"], sod_files),
        "flume-sod-history.txt": (["sod"], sod_history),
        "flume-sod-oracle.txt": (["sod"], lambda: [oracle(R["sod"], "--l1-max", "3.31e-3")]),
        "flume-sod-orders.md": (["sod", "sod-u3", "sod-u7"], sod_orders),
        "flume-sod-ngc.txt": (["sod-u7-ngc3"], lambda: error_line("sod-u7-ngc3")),
        "flume-lax-splitting.txt": (["lax"], lambda: [oracle(R["lax"])]),
        "flume-lax-hllc.txt": (["lax-hllc"], lambda: [oracle(R["lax-hllc"])]),
        "flume-blast-amr-history.txt": (["blast-amr"], blast_history),
        "flume-limiter-none.txt": (["limiter-none"], lambda: error_line("limiter-none")),
        "flume-limiter-cell.txt": (["limiter-cell"], lambda: [" ".join(R["limiter-cell"].rows("-residuals.dat")[-1][:3])
                                                              + " ..."]),
        "prism-pulse-energy.txt": (["pulse"], pulse_energy),
        "prism-pulse-return.txt": (["pulse"], lambda: compare(R["pulse"], 0, 192)),
        "prism-pulse-centred.txt": (["pulse-centred"], lambda: compare(R["pulse-centred"], 0, 192)),
        "prism-faces-pec.md": (["pulse-pec"], lambda: table("pulse-pec", [10, 50, 100, 150, 192])),
        "prism-faces-sm.md": (["pulse-sm"], lambda: table("pulse-sm", [10, 50, 100, 150, 192])),
        "prism-pml-log.txt": (["pulse-pml"], pml_log),
        "prism-pml-energy.md": (["pulse-pml"], lambda: table("pulse-pml", [10, 110, 210, 310, 410,
                                                                           int(R["pulse-pml"].rows("-energy_history.dat")[-1][0])], True)),
        "prism-gpu-energy.txt": (["pulse", "pulse-fnl"], gpu_energy),
        "prism-gpu-compare.txt": (["pulse", "pulse-fnl"], lambda: compare(R["pulse-fnl"], 192, 192,
                                                                          ref=str(R["pulse"].dir))),
        "ripple-periodic-history.txt": (["ripple-periodic"], ripple_periodic),
        "ripple-r1.txt": (["ripple-r1"], lambda: check_lines(R["ripple-r1"], r"^\s*cells", r"PASS \[r1\]")),
        "ripple-r2r3.txt": (["ripple-r2r3"], lambda: check_lines(R["ripple-r2r3"], r"^\s*cells", r"PASS \[r3\]")),
        "ripple-r6.txt": (["ripple-r6"], lambda: check_lines(R["ripple-r6"], r"refined base blocks",
                                                             r"PASS \[r6 unregistered\]")),
        "ripple-r4.txt": (["ripple-r4"], lambda: check_lines(R["ripple-r4"], r"blocks, largest")),
    }


def claims(R):
    """(page, runs, function returning (ok, what)). A claim with a string is a literal the page must contain; a claim
    on round-off (a conservation drift) checks a bound only, since the digits are the compiler's."""
    def has(page, text):
        return (DOCS / page).read_text().find(text) >= 0, f"{page}: '{text}'"

    def bound(what, x, b):
        return x <= b, f"{what} = {x:.3e} <= {b:g}"

    def steps(run, suffix):
        return int(R[run].rows(suffix)[-1][0])

    def drift(run, suffix, col):
        rows = R[run].rows(suffix)
        return max(abs(num(r[col]) - num(rows[0][col])) for r in rows)

    def sod_mass():
        rows = R["sod"].rows("-conservation_history.dat")
        return num(rows[-1][2]), num(rows[-1][3])

    def blast():
        n0, dx0 = xdmf_blocks(sorted(R["blast-amr"].dir.glob("*-000000000.xdmf"))[0])
        n1, dx1 = xdmf_blocks(sorted(R["blast-amr"].dir.glob("*.xdmf"))[-1])
        return n0, n1, dx0, dx1

    pulse_e = lambda r: abs(energy_total(R[r], steps(r, "-energy_history.dat")) / energy_total(R[r], 10) - 1)
    first_dt = lambda r: num(R[r].first("time step:").split()[2])
    return [
        (["sod"], lambda: has("start/first-run.md", f"${tex(l1(R['sod']), 4)}$")),
        (["sod"], lambda: has("tutorials/flume/sod.md", f"${sod_mass()[0]:.4f}$")),
        (["sod"], lambda: bound("Sod mass drift", abs(sod_mass()[0] - 0.5625), 1e-14)),
        (["sod"], lambda: has("tutorials/flume/sod.md", f"${sod_mass()[1]:.2f} = (p_L - p_R)")),
        (["sod"], lambda: has("tutorials/flume/sod.md", f"step {steps('sod', '-conservation_history.dat')}")),
        (["lax"], lambda: has("tutorials/flume/lax.md", f"${tex(l1(R['lax']), 4)}$")),
        (["lax-hllc"], lambda: has("tutorials/flume/lax.md", f"${tex(l1(R['lax-hllc']), 4)}$")),
        (["lax", "lax-hllc"], lambda: has("tutorials/flume/lax.md",
                                          f"by {100 * (l1(R['lax-hllc']) / l1(R['lax']) - 1):.1f}%")),
        (["blast-amr"], lambda: has("tutorials/flume/blast-amr.md",
                                    f"$t = {tex(num(R['blast-amr'].rows('-conservation_history.dat')[-1][1]), 4)}$")),
        (["blast-amr"], lambda: has("tutorials/flume/blast-amr.md", f"both hold {blast()[0]} blocks")),
        (["blast-amr"], lambda: (blast()[0] == blast()[1] and blast()[2] == blast()[3] == [64, 128],
                                 f"blast-amr blocks {blast()[:2]}, 1/dx {blast()[2:]}")),
        (["blast-amr"], lambda: bound("blast-amr mass drift", drift("blast-amr", "-conservation_history.dat", 2), 1e-14)),
        (["limiter-cell"], lambda: has("tutorials/flume/limiter.md", f"in {steps('limiter-cell', '-residuals.dat')} steps")),
        (["pulse"], lambda: has("tutorials/prism/pulse.md", f"${tex(pulse_e('pulse'))}$ relative")),
        (["pulse"], lambda: (all(num(t) == 0 for t in R["pulse"].rows("-divergence_history.dat")[-1][3:5]),
                             "pulse: last divergence row all zero")),
        (["pulse-pec"], lambda: has("tutorials/prism/faces.md", f"${tex(pulse_e('pulse-pec'), 1)}$ over the run")),
        (["pulse-sm"], lambda: has("tutorials/prism/faces.md",
                                   f"${tex(energy_total(R['pulse-sm'], 192) / energy_total(R['pulse-sm'], 10))}$ of the initial")),
        (["pulse-pml"], lambda: has("tutorials/prism/pml.md", f"${tex(first_dt('pulse-pml'), 3)}$ s, against")),
        (["pulse", "pulse-pml"], lambda: has("tutorials/prism/pml.md", f"against ${tex(first_dt('pulse'), 3)}$ s")),
        (["pulse-pml"], lambda: has("tutorials/prism/pml.md",
                                    f"takes {steps('pulse-pml', '-energy_history.dat')} steps")),
        (["pulse-pml"], lambda: has("tutorials/prism/pml.md", "leaves $" + tex(
            energy_total(R["pulse-pml"], steps("pulse-pml", "-energy_history.dat")) / energy_total(R["pulse-pml"], 10))
            + "$ of the energy")),
        (["ripple-periodic"], lambda: has("tutorials/build-your-own/core.md",
                                          f"{int(R['ripple-periodic'].rows('-history.dat')[-1][0])} steps")),
    ]


# ------------------------------------------------------------------------------------------------------- compare
NUM = re.compile(r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][+-]?\d+)?")


def same(a, b, rtol, atol):
    """Equal text, numbers within rtol/atol (the round-off of another compiler)."""
    if NUM.sub("#", a) != NUM.sub("#", b):
        return False
    for x, y in zip(NUM.findall(a), NUM.findall(b)):
        x, y = num(x), num(y)
        if abs(x - y) > max(atol, rtol * max(abs(x), abs(y))):
            return False
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="compare with the committed snippets, write nothing")
    ap.add_argument("--only", help="comma-separated run names (default: all CPU runs)")
    ap.add_argument("--gpu", action="store_true", help="also the FNL runs")
    ap.add_argument("--exe-dir", type=Path, default=REPO / "exe")
    ap.add_argument("--work", type=Path, default=REPO / "exe" / "tutorials-work")
    ap.add_argument("--reuse", action="store_true", help="do not run; read the outputs already in --work")
    ap.add_argument("--with-existing", action="store_true",
                    help="also read the runs not selected by --only whose outputs are already in --work (e.g. --gpu "
                         "--only pulse-fnl,ripple-r4 --with-existing, after a CPU capture)")
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--atol", type=float, default=1e-12)
    ap.add_argument("--list", action="store_true", help="list the runs and snippets, and exit")
    args = ap.parse_args()

    names = [n for n in RUNS if (args.gpu or not RUNS[n][5])]
    if args.only:
        names = [n for n in args.only.split(",") if n in RUNS]
    R = {n: Run(n, args.work.resolve(), args.exe_dir.resolve()) for n in RUNS}
    if args.list:
        for n in RUNS:
            print(f"run {n:16s} {RUNS[n][0]:30s} {RUNS[n][1]}")
        for s, (runs, _) in snippets(R).items():
            print(f"snippet {s:30s} <- {','.join(runs)}")
        return 0

    bad = 0
    for n in names:
        run = R[n]
        if args.reuse:
            run.status = 1 if run.must_fail else 0
            if run.exe == "check.sh":
                run.stdout = (run.dir / "log.txt").read_text()
            continue
        print(f"run {n} ...", flush=True)
        run.execute()
        if not run.ok():
            bad += 1
            print(f"  FAIL {n}: exit {run.status} (expected {'non-zero' if run.must_fail else 0}); log {run.dir}/log.txt")

    done = set(names)
    if args.with_existing:
        for n, run in R.items():
            if n not in done and (run.dir / "log.txt").exists():
                run.status = 1 if run.must_fail else 0
                if run.exe == "check.sh":
                    run.stdout = (run.dir / "log.txt").read_text()
                done.add(n)
    out = Path(tempfile.mkdtemp()) if args.check else SNIPPETS
    out.mkdir(parents=True, exist_ok=True)
    for name, (runs, make) in snippets(R).items():
        if not set(runs) <= done:
            continue
        try:
            text = "\n".join(make()) + "\n"
        except Exception as exc:                                    # a missing file or line: the run changed
            bad += 1
            print(f"  FAIL snippet {name}: {type(exc).__name__}: {exc}")
            continue
        if args.check:
            ref = SNIPPETS / name
            old = ref.read_text().splitlines() if ref.exists() else []
            new = text.splitlines()
            if len(old) != len(new) or not all(same(a, b, args.rtol, args.atol) for a, b in zip(old, new)):
                bad += 1
                print(f"  DIFF {name}\n    committed: " + "\n               ".join(old) +
                      "\n    this run:  " + "\n               ".join(new))
        else:
            (out / name).write_text(text)
    for runs, check in claims(R):
        if not set(runs) <= done:
            continue
        try:
            ok, what = check()
        except Exception as exc:
            ok, what = False, f"{type(exc).__name__}: {exc}"
        if not ok:
            bad += 1
            print(f"  FAIL claim {what}")
    print(f"{len(names)} runs; {'all checks hold' if not bad else f'{bad} failure(s)'}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
