#!/usr/bin/env bash
# pic-periodic: PIC regression legs (CPU np 1 for the first 4 steps, FNL np 1 and np 2), see ../pic_check.sh.
# Usage: ./check.sh <cpu|fnl> [--capture]   (FNL_NP=1 on WSL; np 2 on the cluster)
exec "$(cd "$(dirname "$0")" && pwd)/../pic_check.sh" "$(cd "$(dirname "$0")" && pwd)" pic_periodic 4 "$@"
