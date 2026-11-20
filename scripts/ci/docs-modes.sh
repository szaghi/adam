#!/usr/bin/env bash
# Every `fobis ... --mode X` in the documentation must be a mode the fobos files define (issue #60, P8).
cd "$(dirname "$0")/../.." && exec python3 scripts/check_docs_modes.py
