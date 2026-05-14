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
          { text: 'About',         link: '/guide/' },
          { text: 'Features',      link: '/guide/features' },
          { text: 'Installation',  link: '/guide/installation' },
          { text: 'Architecture',  link: '/guide/architecture' },
          { text: 'Contributing',  link: '/guide/contributing' },
          { text: 'Changelog',     link: '/guide/changelog' },
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
            { text: 'Installation', link: '/guide/installation' },
            { text: 'Architecture', link: '/guide/architecture' },
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
              ],
            },
            {
              text: 'CHASE',
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
            { text: 'Overview',      link: '/tests/' },
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
