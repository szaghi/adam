#!/usr/bin/env bash
# No Runge-Kutta combination written in src/app or src/lib: the time integrators are SELAH's (#132, 6.2). The check's
# self-test first (each planted kind caught, a plain dt product not), then the sources.
cd "$(dirname "$0")/../.." || exit 2
python3 scripts/check_time_separation.py --self-test && exec python3 scripts/check_time_separation.py
