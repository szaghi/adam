#!/usr/bin/env python3
"""Rewrite options of a FLUME input in place (stdlib only).

Why: the speedup scripts generate an input, then switch a few options per run (output cadence, step count) without
regenerating it. The removed matrix generator carried this as its `set` subcommand (`make_inputs.py set`, #133
phase 0); it lives here on its own, for profile-couette3d/profile.sh now and the matrix harness of phase 3 later.

Each item is `section.key=value`; the section is everything before the last dot, so a dotted section name works. A
missing section is added; a missing key is added to its section. Comments are not kept (configparser drops them), and
keys keep their case. Stdlib only: medooza's python3 has no third-party packages.

Usage:
    ini_set.py <input.ini> <section.key=value> [<section.key=value> ...]
"""

from __future__ import annotations

import configparser
import sys
from pathlib import Path


def set_options(ini: Path, items: list[str]) -> None:
    """Rewrite `section.key=value` options of an input in place.

    Parameters
    ----------
    ini : Path
        The input file, rewritten in place (its comments are not kept).
    items : list[str]
        The options, each `section.key=value`.

    Raises
    ------
    SystemExit
        If the input does not exist or an item is not of the form `section.key=value`.
    """
    if not ini.is_file():
        raise SystemExit(f"ini_set: no input '{ini}'")
    parser = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None, strict=False)
    parser.optionxform = str  # keep the case of the keys
    parser.read(ini)
    for item in items:
        key, sep, value = item.partition("=")
        section, dot, option = key.rpartition(".")
        if not (sep and dot and section and option):
            raise SystemExit(f"ini_set: '{item}' is not section.key=value")
        if not parser.has_section(section):
            parser.add_section(section)
        parser[section][option] = value
    with ini.open("w") as f:
        parser.write(f)


def main() -> None:
    """Parse the command line and rewrite the input."""
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    set_options(Path(sys.argv[1]), sys.argv[2:])


if __name__ == "__main__":
    main()
