// docs/.vitepress/theme/index.ts
// The sumi identity (custom.css, BRAND.md) on the default theme, plus the redirects of the pages and sections moved by the
// docs restructure: GitHub Pages has no server redirects, so an old URL lands on the 404 page (or on a page without the
// old anchor) and is replaced here, client side, by its new location (redirects.ts). And the dark figures: in the dark
// theme a page figure <img src=".../x.png"> shows x-dark.png, the silver-ink twin, when it exists (dark-figures.data.ts),
// and stays a washi card otherwise (custom.css).
import DefaultTheme from 'vitepress/theme'
import { h, nextTick, watch } from 'vue'
import { onContentUpdated, useData, useRoute, withBase } from 'vitepress'
import { pageRedirects, anchorRedirects } from './redirects'
import { data as darkFigures } from './dark-figures.data'
import './custom.css'

// '/adam/concepts/library.html#x' -> '/concepts/library'
function sitePath(pathname: string, base: string): string {
  let p = pathname.startsWith(base) ? '/' + pathname.slice(base.length) : pathname
  p = p.replace(/\.html$/, '').replace(/\/index$/, '/')
  if (p.length > 1 && p.endsWith('/')) p = p.slice(0, -1)
  return p || '/'
}

function redirect(isNotFound: boolean, base: string) {
  if (typeof window === 'undefined') return
  const path = sitePath(window.location.pathname, base)
  const hash = window.location.hash.slice(1)
  if (hash) {
    const target = anchorRedirects[`${path}#${decodeURIComponent(hash)}`]
    if (target && (isNotFound || !document.getElementById(decodeURIComponent(hash)))) {
      window.location.replace(withBase(target))
      return
    }
  }
  if (isNotFound) {
    const target = pageRedirects[path]
    if (target) window.location.replace(withBase(target) + (window.location.hash || ''))
  }
}

// The twins are swapped on the client only, after the content is mounted (the server-rendered page always carries the
// light figure, so hydration never mismatches): on every content update (first mount, route change) and theme change.
// The light source is kept in data-light-src; a twin that fails to load falls back to the washi card.
const twins = new Set<string>(darkFigures)
const TWIN = 'adam-dark-twin'

function swapFigures(dark: boolean, base: string) {
  if (typeof document === 'undefined') return
  document.querySelectorAll<HTMLImageElement>('.vp-doc img').forEach((img) => {
    const light = img.dataset.lightSrc ?? img.getAttribute('src') ?? ''
    if (!img.dataset.lightSrc) {
      const path = new URL(light, window.location.href).pathname
      if (!path.startsWith(base) || !twins.has(path.slice(base.length))) return
      img.dataset.lightSrc = light
    }
    if (dark && !img.classList.contains(TWIN)) {
      img.addEventListener('error', () => {
        img.classList.remove(TWIN)
        img.setAttribute('src', light)
      }, { once: true })
      img.classList.add(TWIN)
      img.setAttribute('src', light.replace(/\.png$/, '-dark.png'))
    } else if (!dark && img.classList.contains(TWIN)) {
      img.classList.remove(TWIN)
      img.setAttribute('src', light)
    }
  })
}

const Layout = {
  setup() {
    const { page, site, isDark } = useData()
    const route = useRoute()
    watch(() => [route.path, page.value.isNotFound],
          () => nextTick(() => redirect(!!page.value.isNotFound, site.value.base)),
          { immediate: true })
    onContentUpdated(() => swapFigures(isDark.value, site.value.base))
    watch(isDark, (dark) => nextTick(() => swapFigures(dark, site.value.base)))
    return () => h(DefaultTheme.Layout)
  },
}

export default { extends: DefaultTheme, Layout }
