#!/usr/bin/env python3
"""Turn the regression sod-x input into the Shu-Osher problem on a given scheme_space (issue #119, T2.7). Arguments: input, output, scheme_space."""

from __future__ import annotations

import re
import sys


def section(text: str, head: str, subs: dict[str, str]) -> str:
    """Set keys in one section."""
    start = text.index(head)
    end = text.find("\n[", start + len(head))
    end = len(text) if end < 0 else end
    block = text[start:end]
    for key, value in subs.items():
        block, n = re.subn(rf"^{key}(\s*)=.*$", rf"{key}\g<1>= {value}", block, flags=re.M)
        if n != 1:
            sys.exit(f"make_shu_osher: key {key} not found in {head}")
    return text[:start] + block + text[end:]


def main() -> int:
    """Write the input."""
    source, target, scheme = sys.argv[1:4]
    text = open(source).read()
    text = section(text, "[grid]", {"ni": "100", "emin_x": "-5.0", "emax_x": "5.0"})
    text = section(text, "[numerics]", {"scheme_space": scheme})
    text = section(text, "[initial_conditions]", {"type": "shu-osher"})
    text = text.replace("[initial_conditions]\n", "[initial_conditions]\naxis           = x\ninterface      = -4.0\nrho_amplitude  = 0.2\nrho_wavenumber = 5.0\n", 1)
    text = section(text, "[initial_conditions_region_1]", {"r": "3.857143", "u": "2.629369", "p": "10.33333"})
    text = section(text, "[initial_conditions_region_2]", {"r": "1.0", "u": "0.0", "p": "1.0"})
    text = section(text, "[time]", {"time_max": "1.8"})
    open(target, "w").write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
