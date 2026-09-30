#!/usr/bin/env python3
"""Split a single-realm input along x at the domain centre into a 2-realm forest (issue #41, M2-P7c, MV-14).

Why: the multi-realm verification of #37 (verification/multirealm) glues two realms at x = 0.5 through a mirror seam
filled at every Runge-Kutta stage (beta cadence), where the seam is a block interface like any other: the union of the
realms must reproduce the single-realm run bit for bit. This script derives that pair of realms from any single-realm
input whose x axis carries an even number of cells per block (each realm keeps the block layout and halves ni), so
the MHD cases reuse their generators: realm 1 is [emin_x, centre], realm 2 [centre, emax_x]; the outer faces keep
their boundary conditions and the forest turns the inner ones into the seam.

Writes <out-dir>/<name>.ini (the forest manifest), <name>-r1.ini and <name>-r2.ini; --set overrides an option of
both realms (e.g. time.it_max=100, IO.restart=.true.).

Usage:
    make_split.py <single.ini> <out-dir> <name> [--set SECTION.KEY=VALUE ...]
"""

from __future__ import annotations

import argparse
import configparser
from pathlib import Path

MANIFEST = """[forest]
realms_number = 2

[realm.1]
ini = {name}-r1.ini

[realm.2]
ini = {name}-r2.ini

[forest.topology]
inter_realm_faces_number = 1

[forest.topology.face_1]
realm_a          = 1
face_a           = +x
realm_b          = 2
face_b           = -x
coupling         = mirror
coupling_cadence = stage_coincident
"""


def main() -> None:
    """Write the manifest and the two realm inputs."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("single", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("name")
    parser.add_argument("--set", action="append", default=[], metavar="SECTION.KEY=VALUE")
    args = parser.parse_args()
    for r in (1, 2):
        ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
        ini.optionxform = str
        ini.read(args.single)
        ni = int(ini["grid"]["ni"])
        if ni % 2:
            raise SystemExit(f"make_split: ni = {ni} is odd, the realms cannot halve it")
        xmin, xmax = float(ini["grid"]["emin_x"]), float(ini["grid"]["emax_x"])
        centre = repr(0.5 * (xmin + xmax))
        ini["grid"].update({"ni": str(ni // 2), "emax_x" if r == 1 else "emin_x": centre})
        ini["IO"].update({"output_basename": f"{args.name}-r{r}", "restart_basename": f"{args.name}-r{r}-restart"})
        for item in args.set:
            key, val = item.split("=", 1)
            section, option = key.split(".", 1)
            ini[section][option] = val
        with open(args.out / f"{args.name}-r{r}.ini", "w") as f:
            ini.write(f)
    (args.out / f"{args.name}.ini").write_text(MANIFEST.format(name=args.name))


if __name__ == "__main__":
    main()
