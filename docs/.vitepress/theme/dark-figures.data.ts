// docs/.vitepress/theme/dark-figures.data.ts
// The figures that have a silver-ink twin: every docs/public/<path>/<name>.png next to which <name>-dark.png exists
// (written by the figure scripts with --theme dark, docs/scripts/figure_theme.py). A VitePress data loader: the list is
// computed at build time (and refreshed by the dev server when a twin is added), so a new twin needs no edit here or in
// a page. index.ts swaps the images on the client.
import { existsSync } from 'node:fs'
import { relative, resolve } from 'node:path'

const PUBLIC = resolve(__dirname, '../../public')

export default {
  watch: ['../../public/**/*-dark.png'],
  load(files: string[]): string[] {
    return files
      .map((f) => relative(PUBLIC, f).split('\\').join('/').replace(/-dark\.png$/, '.png'))
      .filter((p) => existsSync(resolve(PUBLIC, p)))
      .sort()
  },
}
