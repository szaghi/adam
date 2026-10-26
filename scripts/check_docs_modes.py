#!/usr/bin/env python3
"""Check that every build mode the documentation names exists in the fobos files.

Why this script exists: the docs audit of issue #60 found 38 references to build modes that did not exist, because a
command written in a page is not compiled. This check reads every `fobis <command> ... --mode X` in the documentation
(docs/**/*.md, README.md and the source READMEs) and fails when X is not a mode `fobis build --lmodes` lists.
Placeholders (`X`, `<mode>`, `$MODE`, names with `<...>`) are skipped. Standard library only; CI runs it (issue #60, P8).

Usage:
    check_docs_modes.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

from make_build_modes_page import ROOT, list_modes

# a fobis command with a --mode option, on one line (inline code or a shell block)
MODE_RE = re.compile(r"\bfobis\s+[a-z]+\b[^`\n]*?--mode[ =]([A-Za-z0-9_.<>$-]+)")


def documents() -> list[Path]:
    """Return the Markdown files that document the build: the site pages, the root README, the source READMEs."""
    skip = ("node_modules", ".vitepress", "third_party")
    pages = [p for p in (ROOT / "docs").rglob("*.md") if not any(s in p.parts for s in skip)]
    pages = [p for p in pages if p.relative_to(ROOT / "docs").parts[0] != "api"]
    readmes = [p for p in (ROOT / "src").rglob("README.md") if not any(s in p.parts for s in skip)]
    return sorted([*pages, ROOT / "README.md", *readmes])


def is_placeholder(name: str) -> bool:
    """Return True for a placeholder rather than a mode name."""
    return name == "X" or any(c in name for c in "<>$")


def main() -> int:
    """Report every documented mode that fobis does not define; exit 1 if there is any."""
    modes = set(list_modes(ROOT / "fobos"))
    if not modes:
        sys.exit("check_docs_modes: fobis listed no mode")
    missing, checked = [], 0
    for page in documents():
        for number, line in enumerate(page.read_text(errors="replace").splitlines(), start=1):
            for name in MODE_RE.findall(line):
                if is_placeholder(name):
                    continue
                checked += 1
                if name not in modes:
                    missing.append(f"{page.relative_to(ROOT)}:{number}: --mode {name}")
    for item in missing:
        print(f"FAIL: {item}")
    print(f"check_docs_modes: {checked} mode references, {len(missing)} not defined in fobis")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
