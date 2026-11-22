#!/usr/bin/env bash
# No WENO, TENO or hybrid operator written in ADAM's sources: the space operators are SETH's (#139). The self-test first
# (each kind planted must be caught), then the check over src/app and src/lib with its waivers (the files the migration
# has not reached yet; the list shrinks with each step of #139 and is empty after S7).
cd "$(dirname "$0")/../.." || exit 2
python3 scripts/check_space_separation.py --self-test && exec python3 scripts/check_space_separation.py
