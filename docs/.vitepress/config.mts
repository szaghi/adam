import { withMermaid } from 'vitepress-plugin-mermaid'
import apiSidebar from '../api/_sidebar.json'

export default withMermaid({
  title: 'ADAM Documentation',
  head: [
    ['style', {}, '.dark-only{display:none}.dark .dark-only{display:block}.dark .light-only{display:none}'],
  ],
  base: '/adam/',
  description: 'Accelerated fluid Dynamics on Adaptive Mesh refinement grids — fluid dynamics SDK and apps for High-Performance Computing — from laptop to exascale device-accelerated superpc',
  markdown: {
    math: true,
    languages: ['fortran-free-form', 'fortran-fixed-form'],
    languageAlias: {
      'fortran': 'fortran-free-form',
      'f90': 'fortran-free-form',
      'f95': 'fortran-free-form',
      'f03': 'fortran-free-form',
      'f08': 'fortran-free-form',
      'f77': 'fortran-fixed-form',
    },
    // The auto-generated API pages (formal/FORD) carry Fortran docstrings
    // with angle-bracket fragments in *plain text* — e.g. the literal
    // `solid <name>` inside a doc comment. markdown-it tokenises `<name>`
    // as an `html_inline` token and emits it verbatim; VitePress then runs
    // the page through the Vue SFC compiler, which sees an unclosed
    // `<name>` tag and aborts with "Element is missing end tag".
    // (Inline code in backticks is already escaped by markdown-it, so it
    // is not the problem — only raw-text html_inline / html_block tokens
    // are.) Escaping `<` / `>` on those tokens keeps them as literal text
    // the Vue compiler never parses as markup.
    config: (md) => {
      const escapeAngles = (s: string) =>
        s.replace(/</g, '&lt;').replace(/>/g, '&gt;')
      for (const rule of ['html_inline', 'html_block'] as const) {
        const fallback =
          md.renderer.rules[rule] ||
          ((tokens, idx) => tokens[idx].content)
        md.renderer.rules[rule] = (tokens, idx, options, env, self) =>
          escapeAngles(fallback(tokens, idx, options, env, self))
      }
    },
  },
  themeConfig: {
    nav: [
      { text: 'Home', link: '/' },
      {
        text: 'Guide',
        items: [
          { text: 'About',                  link: '/guide/' },
          { text: 'Features',               link: '/guide/features' },
          { text: 'Installation',           link: '/guide/installation' },
          { text: 'Architecture',           link: '/guide/architecture' },
          { text: 'Forest (multi-realm)',   link: '/guide/forest' },
          { text: 'Forest cookbook',        link: '/guide/forest-cookbook' },
          { text: 'Contributing',           link: '/guide/contributing' },
          { text: 'Changelog',              link: '/guide/changelog' },
        ],
      },
      { text: 'Library',      link: '/library/' },
      { text: 'Applications', link: '/applications/' },
      { text: 'Tests',        link: '/tests/' },
      { text: 'API',          link: '/api/' },
      { text: 'GitHub',       link: 'https://github.com/szaghi/adam' },
    ],
    sidebar: {
      '/guide/': [
        {
          text: 'Introduction',
          items: [
            { text: 'About',    link: '/guide/' },
            { text: 'Features', link: '/guide/features' },
          ],
        },
        {
          text: 'Getting Started',
          items: [
            { text: 'Installation',       link: '/guide/installation' },
            { text: 'Architecture',       link: '/guide/architecture' },
            { text: 'Forest (multi-realm)', link: '/guide/forest' },
            { text: 'Forest cookbook',      link: '/guide/forest-cookbook' },
          ],
        },
        {
          text: 'Project',
          items: [
            { text: 'Contributing',    link: '/guide/contributing' },
            { text: 'Changelog',       link: '/guide/changelog' },
            { text: 'Coverage Report', link: '/guide/coverage-analysis' },
          ],
        },
      ],
      '/library/': [
        {
          text: 'Library',
          items: [
            { text: 'Overview',          link: '/library/' },
            { text: 'Common',            link: '/library/common' },
            { text: 'FNL (OpenACC)',     link: '/library/fnl' },
            { text: 'NVF (CUDA Fortran)', link: '/library/nvf' },
            { text: 'GMP (OpenMP)',      link: '/library/gmp' },
          ],
        },
      ],
      '/applications/': [
        {
          text: 'Applications',
          items: [
            { text: 'Overview', link: '/applications/' },
            {
              text: 'NASTO',
              collapsed: false,
              items: [
                { text: 'Overview',      link: '/applications/nasto/' },
                { text: 'Common',        link: '/applications/nasto/common' },
                { text: 'CPU Backend',   link: '/applications/nasto/cpu' },
                { text: 'FNL Backend',   link: '/applications/nasto/fnl' },
                { text: 'NVF Backend',   link: '/applications/nasto/nvf' },
                { text: 'GMP Backend',   link: '/applications/nasto/gmp' },
              ],
            },
            {
              text: 'PRISM',
              collapsed: false,
              items: [
                { text: 'Overview',    link: '/applications/prism/' },
                { text: 'Common',      link: '/applications/prism/common' },
                { text: 'CPU Backend', link: '/applications/prism/cpu' },
                { text: 'FNL Backend', link: '/applications/prism/fnl' },
                { text: 'Verification', link: '/applications/prism/verification' },
              ],
            },
            {
              text: 'CHASE (deprecated)',
              collapsed: false,
              items: [
                { text: 'Overview',    link: '/applications/chase/' },
                { text: 'Common',      link: '/applications/chase/common' },
                { text: 'CPU Backend', link: '/applications/chase/cpu' },
              ],
            },
            {
              text: 'PATCH',
              collapsed: false,
              items: [
                { text: 'Overview',    link: '/applications/patch/' },
                { text: 'Common',      link: '/applications/patch/common' },
                { text: 'CPU Backend', link: '/applications/patch/cpu' },
              ],
            },
            {
              text: 'FLUME',
              collapsed: true,
              items: [
                { text: 'Overview',            link: '/applications/flume/' },
                { text: 'Physical Models',     link: '/applications/flume/models' },
                { text: 'Numerical Methods',   link: '/applications/flume/numerics' },
                { text: 'Initial Conditions',  link: '/applications/flume/initial-conditions' },
                { text: 'Boundary Conditions', link: '/applications/flume/boundary-conditions' },
                { text: 'Input Reference',     link: '/applications/flume/input' },
                {
                  text: 'Verification',
                  collapsed: true,
                  items: [
                    { text: 'Overview', link: '/applications/flume/verification/' },
                    {
                      text: 'Euler',
                      link: '/applications/flume/verification/euler/',
                      collapsed: true,
                      items: [
                        { text: 'Sod shock tube', link: '/applications/flume/verification/euler/sod' },
                        { text: 'Lax problem', link: '/applications/flume/verification/euler/lax' },
                        { text: 'Shu–Osher shock / entropy-wave interaction', link: '/applications/flume/verification/euler/shu-osher' },
                        { text: 'Isentropic vortex: order of accuracy', link: '/applications/flume/verification/euler/isentropic-vortex' },
                        { text: 'Shock over a cylinder', link: '/applications/flume/verification/euler/shock-cylinder' },
                        { text: 'Conservation across AMR seams', link: '/applications/flume/verification/euler/amr-conservation' },
                        { text: 'Multi-realm runs', link: '/applications/flume/verification/euler/multi-realm' },
                        { text: 'Forward-facing step: a three-realm forest on quadtrees', link: '/applications/flume/verification/euler/forward-facing-step' },
                        { text: 'Dual energy: the advected entropy of euler-de', link: '/applications/flume/verification/euler/dual-energy' },
                        { text: 'Passive scalars', link: '/applications/flume/verification/euler/passive-scalars' },
                      ],
                    },
                    {
                      text: 'Ideal MHD',
                      link: '/applications/flume/verification/mhd/',
                      collapsed: true,
                      items: [
                        { text: 'MHD on `weno-riemann`', link: '/applications/flume/verification/mhd/weno-riemann' },
                        { text: 'Ryu–Jones 2a: all seven waves', link: '/applications/flume/verification/mhd/ryu-jones-2a' },
                        { text: 'Brio–Wu and Ryu–Jones 4d', link: '/applications/flume/verification/mhd/brio-wu-rj4d' },
                        { text: 'MHD plumbing', link: '/applications/flume/verification/mhd/plumbing' },
                        { text: 'MHD with a zero field', link: '/applications/flume/verification/mhd/zero-field' },
                        { text: 'GLM divergence cleaning', link: '/applications/flume/verification/mhd/glm-cleaning' },
                        { text: 'Order of accuracy: linear waves, CPAW, magnetised vortex', link: '/applications/flume/verification/mhd/order-of-accuracy' },
                        { text: 'MHD conservation across AMR seams', link: '/applications/flume/verification/mhd/conservation' },
                        { text: 'Orszag–Tang vortex', link: '/applications/flume/verification/mhd/orszag-tang' },
                        { text: 'MHD rotor', link: '/applications/flume/verification/mhd/rotor' },
                        { text: 'Field loop advection', link: '/applications/flume/verification/mhd/field-loop' },
                        { text: 'Rotated shock tube', link: '/applications/flume/verification/mhd/rotated-shock-tube' },
                        { text: 'EGLM cleaning', link: '/applications/flume/verification/mhd/eglm' },
                        { text: 'Positivity limiter', link: '/applications/flume/verification/mhd/positivity-limiter' },
                        { text: 'Quadtree AMR', link: '/applications/flume/verification/mhd/quadtree-amr' },
                      ],
                    },
                    {
                      text: 'Infrastructure',
                      link: '/applications/flume/verification/infrastructure/',
                      collapsed: true,
                      items: [
                        { text: 'Ghost probe: every face and edge ghost on a linear field', link: '/applications/flume/verification/infrastructure/ghost-probe' },
                        { text: 'Scaling covariance and the reference layer', link: '/applications/flume/verification/infrastructure/scaling-covariance' },
                        { text: 'Ranks without blocks', link: '/applications/flume/verification/infrastructure/ranks-without-blocks' },
                        { text: 'Non-finite state guard', link: '/applications/flume/verification/infrastructure/nonfinite' },
                      ],
                    },
                    {
                      text: 'Dissipative terms (M4)',
                      link: '/applications/flume/verification/dissipative/',
                      collapsed: true,
                      items: [
                        { text: 'Input contract of the dissipative terms', link: '/applications/flume/verification/dissipative/input-contract' },
                        { text: 'Viscous waves', link: '/applications/flume/verification/dissipative/viscous-waves' },
                        { text: 'Compressible Couette flow', link: '/applications/flume/verification/dissipative/couette' },
                        { text: 'Becker\'s viscous shock', link: '/applications/flume/verification/dissipative/becker-shock' },
                        { text: 'Ohmic decay and the visco-resistive Alfvén wave', link: '/applications/flume/verification/dissipative/ohmic-alfven' },
                        { text: 'Ohmic heating and the divergence constraint', link: '/applications/flume/verification/dissipative/ohmic-heating' },
                        { text: 'Conservation across 2:1 seams and forest seams', link: '/applications/flume/verification/dissipative/seam-conservation' },
                        { text: 'Accuracy across 2:1 seams', link: '/applications/flume/verification/dissipative/seam-accuracy' },
                        { text: 'CPU and FNL agree run by run', link: '/applications/flume/verification/dissipative/cpu-fnl-agreement' },
                      ],
                    },
                    {
                      text: 'Runtime AMR (M5)',
                      link: '/applications/flume/verification/runtime-amr/',
                      collapsed: true,
                      items: [
                        { text: 'Regrid round trip and input contract', link: '/applications/flume/verification/runtime-amr/regrid-round-trip' },
                        { text: 'Runtime regrids in FLUME', link: '/applications/flume/verification/runtime-amr/runtime-regrids' },
                        { text: 'Accuracy of runtime AMR', link: '/applications/flume/verification/runtime-amr/accuracy' },
                        { text: 'MHD blast symmetry on the flux-splitting path', link: '/applications/flume/verification/runtime-amr/mhd-blast-symmetry' },
                      ],
                    },
                    {
                      text: 'Seam accuracy and divergence (M6)',
                      link: '/applications/flume/verification/seams/',
                      collapsed: true,
                      items: [
                        { text: 'Order across static 2:1 seams', link: '/applications/flume/verification/seams/order-across-seams' },
                        { text: 'Divergence across seams', link: '/applications/flume/verification/seams/divergence-across-seams' },
                        { text: 'Shock crossing a static 2:1 face', link: '/applications/flume/verification/seams/shock-crossing' },
                        { text: 'Integral drift and the WENO sensor near seams', link: '/applications/flume/verification/seams/integral-drift' },
                      ],
                    },
                    {
                      text: 'Seam studies (M6)',
                      link: '/applications/flume/verification/seams-studies/',
                      collapsed: true,
                      items: [
                        { text: 'Quintic coarse-to-fine fill', link: '/applications/flume/verification/seams-studies/quintic-fill' },
                        { text: 'Limited coarse-to-fine fills', link: '/applications/flume/verification/seams-studies/limited-fills' },
                        { text: 'Shock-crossing defect is a capturing error', link: '/applications/flume/verification/seams-studies/capturing-error' },
                        { text: 'Reflux policy: troubled against always', link: '/applications/flume/verification/seams-studies/reflux-policy' },
                        { text: '#78 is the reflux, not the restricted ghosts', link: '/applications/flume/verification/seams-studies/reflux-not-ghosts' },
                        { text: 'GLM damping length on a seam', link: '/applications/flume/verification/seams-studies/glm-damping' },
                        { text: 'EGLM reflux leaves B alone', link: '/applications/flume/verification/seams-studies/eglm-reflux' },
                        { text: '2:1 caps in a 2-D prototype', link: '/applications/flume/verification/seams-studies/caps-prototype' },
                        { text: 'Sensor-gated reflux prototype', link: '/applications/flume/verification/seams-studies/sensor-gated-reflux' },
                        { text: 'Footprint of a high-order restriction', link: '/applications/flume/verification/seams-studies/high-order-restriction' },
                      ],
                    },
                    {
                      text: 'External gravity (M7)',
                      link: '/applications/flume/verification/gravity/',
                      collapsed: true,
                      items: [
                        { text: 'Stratified atmospheres and the `hydrostatic` boundary', link: '/applications/flume/verification/gravity/atmospheres' },
                        { text: 'The gravity force', link: '/applications/flume/verification/gravity/force' },
                        { text: 'The order of the energy source', link: '/applications/flume/verification/gravity/energy-order' },
                        { text: 'The well-balanced mode', link: '/applications/flume/verification/gravity/well-balanced' },
                        { text: 'Rayleigh–Taylor', link: '/applications/flume/verification/gravity/rayleigh-taylor' },
                      ],
                    },
                    {
                      text: 'Self-gravity (M8)',
                      link: '/applications/flume/verification/self-gravity/',
                      collapsed: true,
                      items: [
                        { text: 'Gravity from a cell potential', link: '/applications/flume/verification/self-gravity/potential' },
                        { text: 'The FFT Poisson solver', link: '/applications/flume/verification/self-gravity/fft-poisson' },
                        { text: 'Jeans modes under self-gravity', link: '/applications/flume/verification/self-gravity/jeans' },
                        { text: 'The isolated Poisson solver', link: '/applications/flume/verification/self-gravity/isolated-solver' },
                        { text: 'The n = 1 polytrope: equilibrium and collapse', link: '/applications/flume/verification/self-gravity/polytrope' },
                        { text: 'The multigrid Poisson solver on AMR grids', link: '/applications/flume/verification/self-gravity/multigrid-amr' },
                        { text: 'Jeans modes on a static AMR grid', link: '/applications/flume/verification/self-gravity/amr-jeans' },
                        { text: 'Self-gravity under runtime AMR', link: '/applications/flume/verification/self-gravity/runtime-amr' },
                      ],
                    },
                    { text: 'Unit tests', link: '/applications/flume/verification/unit-tests' },
                    { text: 'Regression suite', link: '/applications/flume/verification/regression' },
                  ],
                },
                { text: 'Common',      link: '/applications/flume/common' },
                { text: 'CPU Backend', link: '/applications/flume/cpu' },
                { text: 'FNL Backend', link: '/applications/flume/fnl' },
              ],
            },
            { text: 'ASCOT', link: '/applications/ascot' },
          ],
        },
      ],
      '/api/': [
        {
          text: 'API Reference',
          items: [
            { text: 'Overview', link: '/api/' },
          ],
        },
        ...apiSidebar,
      ],
      '/tests/': [
        {
          text: 'Tests',
          items: [
            { text: 'Overview',            link: '/tests/' },
            { text: 'PRISM regression',    link: '/tests/prism-regression' },
            { text: 'Regression bring-up', link: '/tests/prism-regression-tutorial' },
          ],
        },
        {
          text: 'NASTO integration',
          items: [
            { text: 'Sod-X',         link: '/tests/sod-x' },
            { text: 'Sod-Y',         link: '/tests/sod-y' },
            { text: 'Sod-Z',         link: '/tests/sod-z' },
            { text: 'Shock-Sphere',  link: '/tests/shock-sphere' },
          ],
        },
      ],
    },
    search: {
      provider: 'local',
    },
  },
  mermaid: {},
  vite: {
    optimizeDeps: {
      include: ['mermaid'],
    },
  },
})
