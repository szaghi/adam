#!/usr/bin/env bash
# FLUME speedup matrix: analyse every measurement present (any leg, any np) into results/summary.{md,csv}, a result.md
# per measurement, and the Results section of every case README. Stdlib python only, so it also runs on medooza.
#
# Usage: ./analyze.sh [--ref LABEL]   LABEL: the reference of the backend speedups (default: the cpu label with the
#                                      most ranks), e.g. --ref cpu-np32
set -euo pipefail
exec python3 "$(cd "$(dirname "$0")" && pwd)/analyze.py" "$@"
