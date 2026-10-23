#!/usr/bin/env bash
# FLUME verification of self-gravity (issue #81, M8). SV-2 (P2): the FFT Poisson solver on uniform periodic grids;
# SV-6 (P4): the same solver with isolated boundaries; SV-9 (P5a): the multigrid solver on AMR grids.
#
# Why: self-gravity conserves E + rho Phi and momentum exactly only if Phi solves the discrete Poisson equation
# L_h Phi = 4 pi G (rho - mean rho) to round-off (#81 P0). The library solver (adam_poisson_fft_object, FFTW; its FNL
# twin, cuFFT on rank 0's GPU) divides each Fourier mode by the eigenvalue of L_h, gathers the blocks to rank 0 and
# scatters Phi back with every ghost layer filled by the periodic wrap. The unit test test_poisson_fft holds it to:
#   - the discrete solution (two oblique modes on [0,1] x [0,2] x [0,3], each an eigenfunction of L_h), every cell and
#     ghost, at round-off (1e-12 relative);
#   - the discrete residual from the scattered ghosts (1e-10 relative);
#   - the continuum solution, second order (>= 1.9 on the finest pair, N = 32 -> 64);
# on a 3-D octree, a 2-D quadtree and a single block (ranks without blocks when np > 1). Each rank count of the list
# must pass; a non-periodic grid must be refused (exit 1, the message naming the direction).
# Measured (CPU and FNL): discrete <= 2.3e-15, residual <= 4.3e-14, orders 1.98 (3-D) and 1.98 (2-D); the printed
# numbers do not depend on the rank count.
#
# SV-6 (test_poisson_fft --isolated): Hockney-Eastwood zero padding with the lattice Green's function of L_h
# (adam_lattice_green_library, scripts/lattice_green.py), phi -> 0 at infinity, on the cube [-1,1]^3 (16^3 to 64^3):
#   - an n = 1 polytrope (smooth) and a uniform sphere (cell averages), radius 1/2, against their exact potentials,
#     every cell and ghost: least-squares order >= 1.9 (polytrope) and >= 1.2 (sphere: its surface cuts cells);
#   - the residual of L_h, 1e-10 relative (the kernel's own residual summed over the source, not round-off);
#   - refused, each with its message: a periodic grid, non-cubic cells, a 2-D grid.
# Measured (CPU): orders 2.03 (polytrope) and 1.35 (sphere), residual <= 2.1e-11 at 64^3.
#
# SV-9 (test_poisson_mg, P5a): adam_poisson_mg_object, FCG preconditioned by a V-cycle, on the composite Laplacian
# with the mimetic mean flux across 2:1 faces (symmetric and consistent), SV-2's source and box:
#   - SV-9a, uniform grids: the multigrid equals the FFT solver (1e-8 relative), in at most one iteration;
#   - SV-9b, the central box refined once (octree and quadtree, base levels 2 and 3): FCG residual below the
#     tolerance, symmetry defect of the assembled rows <= 1e-12, every face ghost filled, Phi order >= 1.7 against
#     the exact potential, max truncation error growing at most 1.5x under refinement (bounded: consistent);
#   - SV-12 (issue #87, P1), Dirichlet faces (no active direction periodic), the exact potential as the ghost values:
#     SV-12a uniform, at most one iteration (the coarsest level's odd extension is the whole solve), Phi order >= 1.7;
#     SV-12b as SV-9b with Dirichlet rows;
#   - refused: mixed periodic and non-periodic directions without bc (with bc they are accepted since issue #84,
#     E1b), a periodic axis in bc that the grid bounds, an axis with a ghost-Dirichlet face and a face kind.
# Measured (CPU, np 1-3, rank-independent to round-off): SV-9a 3e-15; SV-9b symmetry <= 8e-18, Phi order 1.92 (3-D
# and 2-D), truncation growth 1.005, 2:1 face gravity order 1.03 (3-D) and 1.11 (2-D), 24-39 FCG iterations;
# SV-12a 1 iteration, order 1.95; SV-12b symmetry <= 6.7e-18, Phi order 1.91 (3-D) and 1.93 (2-D), 22-36 iterations.
#
# SV-13 (test_poisson_mg --isolated, issue #87, P2): isolated boundaries on AMR by the symmetric James step (two
# Dirichlet solves around the ghost-to-ghost lattice Green's function, a padded FFT on the coarsest box widened by the
# ghost layer), SV-6's stars on the cube [-1,1]^3:
#   - SV-13a, uniform 16^3 and 32^3: the multigrid equals the FFT solver's isolated mode (1e-9 relative), each solve
#     in at most one iteration; SV-13b, the star's box [-1/2,1/2]^3 refined once over base levels 2 and 3: Phi order
#     >= 1.7 (polytrope) and >= 1.2 (uniform sphere); both: every face ghost filled, the solution map symmetric in
#     situ (|x . V K y - y . V K x| / |x . V K y| <= 1e-9, FCG tolerance 1e-12);
#   - refused: a leaf above the coarsest level on the boundary.
# Measured (CPU, np 1-3): SV-13a |Phi_mg - Phi_fft| <= 2.6e-13, 1 + 1 iterations, order 2.14; SV-13b orders 2.01
# (polytrope) and 1.86 (sphere), 23-25 iterations per solve; symmetry <= 1.7e-11; boundary ghosts within 1e-4 of the
# exact potential (2.5e-5 at base level 3).
#
# CPU = FNL: run this script once per build with --dump DIR, then
#   poisson_oracle.py compare DIR_CPU DIR_FNL      (max |phi_CPU - phi_FNL| / max |phi|, 1e-13)
#
# bc (test_flail_bc, issue #84, E1a and E2a): the FFT solver with a boundary kind per face (periodic, Dirichlet or
# Neumann at the face, ghost-Dirichlet), the compact operator and the wide D(G) of order 2, 4, 6, 8 (PRISM's elliptic
# operator, E2a): per-axis discrete eigenfunctions (the wide eigenvalue computed in real space, D1 applied twice),
# interior and ghost layers to depth p at round-off, the energy, the D/D order (compact >= 1.9, wide p = 4 >= 3.8);
# refused: a periodic face against a Dirichlet one, a ghost-Dirichlet face against a face kind, the wide operator
# with isolated boundaries.
#
# bcmg (test_poisson_mg --bc, issue #84, E1b): the multigrid with a boundary kind per face (periodic, Dirichlet or
# Neumann at the face, ghost-Dirichlet; mixed across and within axes):
#   - uniform grids: equal to the FFT solver with the same kinds (bc leg, E1a), 1e-8 relative, at most one iteration;
#   - static AMR, the corner octant [0, L/2] refined (leaves of both levels on the low faces), a manufactured potential
#     prod_d cos(2 pi x_d / L_d + theta_d) whose boundary values are the ghost values of each kind: residual,
#     symmetry <= 1e-12, every face ghost filled, Phi order >= 1.7 (up to a constant without a Dirichlet face), the
#     domain-face ghosts' order >= 1.7 against the potential at the ghost centres, truncation growth <= 1.5.
# Measured (CPU, np 1-3, iterations rank-independent): uniform <= 1.2e-14, 1 iteration; AMR orders of Phi 1.95
# (P/P D/N N/D), 2.00 (D/D N/N G/G), 1.95 (P/P N/N N/N), 2.00 (2-D N/D P/P), of the ghosts 1.95-2.03, truncation
# growth 1.00-1.02, symmetry <= 2.5e-17, 21-40 FCG iterations. Each mutation alone fails: the Neumann weight in B
# negated (orders 0.02-0.14), the Neumann ghost's h g negated (ghost orders ~1), the Dirichlet ghost 2 g + phi_C
# (ghost orders ~0.03).
#
# flail (test_poisson_mg --flail, issue #84, E1c): the FLAIL facade with boundary kinds, flail%setup(..., bc) and
# solve against the solver it chooses called directly, bitwise on every cell and ghost (FFT on a uniform grid, the
# multigrid on AMR with inhomogeneous boundary values); a vector with per-component kinds (a PEC face normal to z)
# as three FLAIL instances, solved interleaved, bitwise equal to three solved alone; refused: boundary values on the
# FFT path, isolated with bc.
#
# Usage: ./poisson.sh [--np-list "1 2 3"] [--legs LIST] [--dump DIR] [--keep]
#        LIST: sv2, sv6, sv9, sv13, bc, bcmg, flail or all
#        (default all)
#
# POISSON_EXE overrides the executable under test, e.g. POISSON_EXE=$REPO/exe/test_poisson_fft_fnl ./poisson.sh
# (built by `fobis build --mode test-poisson-fft-gnu` or `test-poisson-fft-fnl-nvf`); POISSON_MG_EXE the SV-9/13 one
# (`test-poisson-mg-gnu`, CPU only; --dump covers SV-2 and SV-6). The caller owns the matching
# environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${POISSON_EXE:-$REPO_ROOT/exe/test_poisson_fft}"
MG_EXE="${POISSON_MG_EXE:-$REPO_ROOT/exe/test_poisson_mg}"
BC_EXE="${POISSON_BC_EXE:-$REPO_ROOT/exe/test_flail_bc}"
NP_LIST="1 2 3"
source "$VERIF_DIR/legs.sh"
legs_declare "sv2 sv6 sv9 sv13 bc bcmg flail"
DUMP=""
KEEP=0

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np-list) NP_LIST="$2" ; shift 2 ;;
      --legs)    legs_select "$2" ; shift 2 ;;
      --dump)    DUMP="$2" ; shift 2 ;;
      --keep)    KEEP=1 ; shift ;;
      *)         echo "poisson.sh: unknown argument '$1' (accepted: --np-list LIST, --legs LIST, --dump DIR, --keep)" >&2
                 exit 2 ;;
   esac
done
legs_resolve poisson.sh
needs_fft=0 ; needs_mg=0 ; needs_bc=0
for leg in $LEGS; do
   case "$leg" in sv2|sv6) needs_fft=1 ;; sv9|sv13|bcmg|flail) needs_mg=1 ;; bc) needs_bc=1 ;; esac
done
if [[ -n "$DUMP" ]]; then needs_fft=1 ; fi
if [[ $needs_fft -eq 1 ]]; then
   if [[ ! -x "$EXE" ]]; then
      echo "poisson.sh: executable '$EXE' not found (build it or set POISSON_EXE)" >&2 ; exit 2
   fi
   EXE="$(realpath "$EXE")"
fi
if [[ $needs_mg -eq 1 ]]; then
   if [[ ! -x "$MG_EXE" ]]; then
      echo "poisson.sh: executable '$MG_EXE' not found (build it or set POISSON_MG_EXE)" >&2 ; exit 2
   fi
   MG_EXE="$(realpath "$MG_EXE")"
fi
if [[ $needs_bc -eq 1 ]]; then
   if [[ ! -x "$BC_EXE" ]]; then
      echo "poisson.sh: executable '$BC_EXE' not found (build test-flail-bc-* or set POISSON_BC_EXE)" >&2 ; exit 2
   fi
   BC_EXE="$(realpath "$BC_EXE")"
fi

TAG="$(basename "$EXE")"
work="$CASE_DIR/work-$TAG-poisson"
rm -rf "$work" ; mkdir -p "$work"
FAILED=0
suite() { # suite <label> <flags> [exe]: the test on every rank count of the list, its lines and verdict
   local label="$1" flags="$2" exe="${3:-$EXE}" np log
   for np in $NP_LIST; do
      log="$work/log-$label-np$np.txt"
      # shellcheck disable=SC2086
      if ! (cd "$work" && mpirun -np "$np" "$exe" $flags > "$log" 2>&1); then
         echo "   np $np: the run failed, see $log  FAIL"
         FAILED=1
         continue
      fi
      grep '^   ' "$log" | sed "s/^/   np $np:/"
      if ! grep -qE '^test_(poisson_(fft|mg)|flail_bc): PASSED' "$log"; then FAILED=1 ; fi
   done
}
refusal() { # refusal <label> <flag> <expected message> [exe]: the run must stop, naming the cause
   local label="$1" flag="$2" msg="$3" exe="${4:-$EXE}" log="$work/log-refusal-$1.txt"
   if (cd "$work" && mpirun -np 2 "$exe" "$flag" > "$log" 2>&1); then
      echo "   refusal, $label: accepted  FAIL"
      FAILED=1
   elif grep -q "$msg" "$log"; then
      echo "   refusal, $label: stops the run ($msg)  PASS"
   else
      echo "   refusal, $label: the run failed without the expected message, see $log  FAIL"
      FAILED=1
   fi
}
for leg in $LEGS; do
   case "$leg" in
      sv2)
         echo ">> SV-2: the FFT Poisson solver, periodic ($TAG)"
         suite periodic ""
         refusal "a non-periodic grid" --no-periodic 'direction 1 is not periodic' ;;
      sv6)
         echo ">> SV-6: the FFT Poisson solver, isolated ($TAG)"
         suite isolated "--isolated"
         refusal "isolated on a periodic grid" --refuse-isolated-periodic 'direction 1 is periodic'
         refusal "isolated with non-cubic cells" --refuse-isolated-anisotropic 'isolated boundaries need cubic cells'
         refusal "isolated in 2-D" --refuse-isolated-2d 'isolated boundaries are 3-D only' ;;
      sv9)
         echo ">> SV-9: the multigrid Poisson solver, AMR ($(basename "$MG_EXE"))"
         suite multigrid "" "$MG_EXE"
         refusal "multigrid, mixed boundaries without bc" --mixed 'mixed periodic and non-periodic directions' \
            "$MG_EXE"
         refusal "multigrid, a periodic axis in bc that the grid bounds" --refuse-bc-mismatch 'in the grid and' \
            "$MG_EXE"
         refusal "multigrid, a ghost-Dirichlet face against a face kind" --refuse-bc-mixed-ghost \
            'mixes a ghost-Dirichlet face' "$MG_EXE" ;;
      sv13)
         echo ">> SV-13: the multigrid Poisson solver, isolated on AMR ($(basename "$MG_EXE"))"
         suite multigrid-isolated "--isolated" "$MG_EXE"
         refusal "isolated, a refined leaf on the boundary" --refuse-boundary-refined 'touches the domain boundary' \
            "$MG_EXE" ;;
      bc)
         echo ">> FLAIL boundary kinds per face and the wide D(G), the FFT solver (issue #84, E1a, E2a; $(basename "$BC_EXE"))"
         suite flail-bc "" "$BC_EXE"
         refusal "a periodic face against a Dirichlet one" --refuse-mixed-periodic 'one periodic face and one' "$BC_EXE"
         refusal "a ghost-Dirichlet face against a face kind" --refuse-mixed-ghost 'mixes a ghost-Dirichlet face' \
            "$BC_EXE"
         refusal "the wide D(G) with isolated boundaries" --refuse-wide-isolated 'the compact operator only' \
            "$BC_EXE" ;;
      bcmg)
         echo ">> FLAIL boundary kinds per face, the multigrid solver (issue #84, E1b; $(basename "$MG_EXE"))"
         suite multigrid-bc "--bc" "$MG_EXE" ;;
      flail)
         echo ">> The FLAIL facade with boundary kinds (issue #84, E1c; $(basename "$MG_EXE"))"
         suite flail "--flail" "$MG_EXE"
         refusal "facade, boundary values on the FFT path" --refuse-flail-ghost-fft 'need the multigrid path' \
            "$MG_EXE"
         refusal "facade, isolated with bc" --refuse-flail-isolated-bc 'isolated boundaries take no bc' "$MG_EXE"
         refusal "facade, the wide D(G) on AMR" --refuse-flail-wide-amr 'the wide D(G) needs a uniform grid' "$MG_EXE" ;;
   esac
done
if [[ -n "$DUMP" ]]; then
   mkdir -p "$DUMP"
   DUMP="$(cd "$DUMP" && pwd)"
   for flags in "" "--isolated"; do
      # shellcheck disable=SC2086
      if ! (cd "$work" && mpirun -np 2 "$EXE" $flags --dump "$DUMP" > "$work/log-dump.txt" 2>&1); then
         echo "   dump ${flags:-periodic}: the run failed, see $work/log-dump.txt  FAIL"
         FAILED=1
      fi
   done
   echo "   dump: $(ls "$DUMP"/phi-*.bin | wc -l) potentials in $DUMP"
fi
if [[ $KEEP -eq 0 && $FAILED -eq 0 ]]; then rm -rf "$work" ; fi
if [[ $FAILED -ne 0 ]]; then
   echo "Poisson legs ($LEGS) FAILED ($TAG)" ; exit 1
fi
echo "Poisson legs ($LEGS) PASSED ($TAG)"
