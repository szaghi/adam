// Moved pages and sections of the docs restructure (2026-11): old site path -> new site path (with an anchor when the
// content now lives in a section). Paths are relative to the site base, without .html. Used by index.ts.
// Keep entries forever: external links and bookmarks still point at the old locations.

export const pageRedirects: Record<string, string> = {
  // developer section
  '/concepts/library': '/developer/library',
  '/concepts/library-fnl': '/developer/fnl',
  '/applications/flume/common': '/developer/flume-source#common',
  '/applications/flume/cpu': '/developer/flume-source#cpu-backend',
  '/applications/flume/fnl': '/developer/flume-source#fnl-backend',
  // merges
  '/concepts/features': '/concepts/architecture',
  '/start': '/',
  '/tests': '/reference/regression',
  '/tests/prism-regression': '/reference/regression#prism-suite',
  '/tests/prism-regression-tutorial': '/start/install#new-machine',
  '/applications/flume/verification/regression': '/reference/regression#flume-suite',
}

// 'old/page#old-anchor' -> 'new/page#new-anchor', for sections moved to another page or renamed.
export const anchorRedirects: Record<string, string> = {
  // developer section (old anchors arrive on the new page with the page redirect)
  '/developer/library#adam-tree-object-—-octree-quadtree-hierarchy': '/developer/library#adam-tree-object',
  '/developer/library#adam-tree-node-object-—-individual-tree-node': '/developer/library#adam-tree-node-object',
  '/developer/library#adam-tree-bucket-object-—-hash-map-bucket': '/developer/library#adam-tree-bucket-object',
  '/developer/library#adam-grid-object-—-block-structured-cartesian-grid': '/developer/library#adam-grid-object',
  '/developer/library#adam-field-object-—-block-geometry-ghost-exchange-refinement': '/developer/library#adam-field-object',
  '/developer/library#adam-maps-object-—-index-translation-and-mpi-communication-patterns': '/developer/library#adam-maps-object',
  '/developer/library#adam-weno-object-—-weno-reconstruction': '/developer/library#adam-weno-object',
  '/developer/library#adam-fdv-operators-library-—-finite-difference-volume-operators': '/developer/library#adam-fdv-operators-library',
  '/developer/library#adam-riemann-euler-library-—-euler-equation-riemann-solvers': '/developer/library#adam-riemann-euler-library',
  '/developer/library#adam-rk-object-—-explicit-runge-kutta-integration': '/developer/library#adam-rk-object',
  '/developer/library#adam-leapfrog-object-—-leapfrog-with-raw-filter': '/developer/library#unused-integrators',
  '/developer/library#adam-blanes-moan-object-—-blanes–moan-splitting-integrator': '/developer/library#unused-integrators',
  '/developer/library#adam-cfm-object-—-commutator-free-magnus-integrators': '/developer/library#unused-integrators',
  '/developer/library#adam-eos-ic-object-—-ideal-compressible-gas-eos': '/developer/library#adam-eos-ic-object',
  '/developer/library#adam-ib-object-—-immersed-boundary-method': '/developer/library#adam-ib-object',
  '/developer/library#adam-amr-object-—-amr-refinement-markers': '/developer/library#adam-amr-object',
  '/developer/library#adam-mpih-object-—-mpi-wrapper': '/developer/library#adam-mpih-object',
  '/developer/library#adam-io-object-—-parallel-i-o': '/developer/library#adam-io-object',
  '/developer/library#adam-slices-object-—-2d-diagnostic-slices': '/developer/library#adam-slices-object',
  '/developer/library#adam-adam-object-—-top-level-orchestrator': '/developer/library#adam-adam-object',
  '/developer/library#adam-parameters-—-global-constants': '/developer/library#adam-parameters',
  '/developer/library#adam-common-library-—-aggregate-entry-point': '/developer/library#module-summary',
  '/developer/library#adam-realm-object-—-the-per-realm-contract': '/developer/library#adam-realm-object',
  '/developer/library#adam-forest-object-—-the-orchestrator': '/developer/library#adam-forest-object',
  '/developer/library#intra-realm-amr-seams-a-second-independent-seam-family': '/developer/library#intra-realm-amr-seams',
  '/developer/library#seam-div-b-guard-rail-io-issue-29': '/developer/library#seam-div-b-guard-rail',
  '/developer/library#tracked-through': '/developer/library#forest-orchestration-multi-realm',
  '/developer/library#adam-flail-object-—-flail-the-elliptic-solver-library-issue-84': '/developer/flail#the-facade',
  '/developer/library#legacy-flail-smoothers-deleted-in-issue-84-e4': '/developer/flail',
  '/developer/library#adam-poisson-fft-object-—-fft-poisson-solver-uniform-grids-periodic-dirichlet-neumann-faces-or-isolated': '/developer/flail#fft-solver',
  '/developer/library#adam-poisson-mg-object-—-poisson-solver-on-amr-grids-multigrid-preconditioned-cg': '/developer/flail#multigrid-solver',
  '/developer/library#adam-lattice-green-library-—-lattice-green-s-function-of-the-7-point-laplacian-generated': '/developer/flail#lattice-green-s-function',
  '/developer/fnl#cpu-vs-gpu-array-layout': '/developer/fnl#cpu-and-gpu-array-layout',
  '/developer/fnl#gpu-direct-ghost-exchange-and-the-wsl2-rendezvous-abort-issue-12': '/developer/fnl#gpu-direct-ghost-exchange',
  '/developer/fnl#poisson-solver': '/developer/fnl#poisson-solvers',
  '/developer/flume-source#build-modes': '/developer/flume-source#cpu-build-modes',
  // merges
  '/tests/prism-regression#two-reference-files-per-case': '/reference/regression#goldens-and-tolerance',
  '/tests/prism-regression#case-inventory': '/reference/regression#prism-cases',
  '/applications/flume/verification/regression#cases': '/reference/regression#flume-cases',
  '/applications/flume/verification/regression#refreshing-a-golden': '/cookbook/regression#t3-add-a-case-and-capture-its-goldens',
  // PRISM hazards page
  '/applications/prism/input#hazards': '/applications/prism/hazards',
  '/applications/prism/input#divergence-correction': '/applications/prism/hazards#divergence-correction',
  '/applications/prism/input#pml-1': '/applications/prism/hazards#pml-constraints',
  '/applications/prism/input#elliptic-boundary-map': '/applications/prism/hazards#elliptic-boundary-map',
  '/applications/prism/input#refused-combinations': '/applications/prism/hazards#refused-combinations',
  '/applications/prism/numerics#the-2-1-seam-and-div-b': '/applications/prism/hazards#seam-div-b',
}
