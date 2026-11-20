# Sourced by the ADAM CI legs (scripts/ci/*.sh), run by fortran-toolchains' matrix.yml in the linux-gnu image from the
# repository root: the dependencies fetched, the image's Python (/opt/venv: h5py, numpy) as the regression venv
# (exe/.regression-venv, which the regression scripts would otherwise create with pip), and ARTIFACTS, the directory
# matrix.yml uploads whatever the outcome.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
fobis fetch
mkdir -p exe
if [[ -x /opt/venv/bin/python && ! -e exe/.regression-venv ]]; then ln -s /opt/venv exe/.regression-venv; fi
ARTIFACTS=ci-artifacts
mkdir -p "$ARTIFACTS"
collect() { # copy the files matching the globs into $ARTIFACTS, keeping their paths
   local f
   shopt -s globstar nullglob
   for f in "$@"; do mkdir -p "$ARTIFACTS/$(dirname "$f")" && cp "$f" "$ARTIFACTS/$f"; done
   shopt -u globstar nullglob
}
