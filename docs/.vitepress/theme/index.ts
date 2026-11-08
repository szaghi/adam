// docs/.vitepress/theme/index.ts
// The sumi identity (custom.css, BRAND.md) on the default theme, plus the redirects of the pages and sections moved by the
// docs restructure: GitHub Pages has no server redirects, so an old URL lands on the 404 page (or on a page without the
// old anchor) and is replaced here, client side, by its new location (redirects.ts).
import DefaultTheme from 'vitepress/theme'
import { h, nextTick, watch } from 'vue'
import { useData, useRoute, withBase } from 'vitepress'
import { pageRedirects, anchorRedirects } from './redirects'
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

const Layout = {
  setup() {
    const { page, site } = useData()
    const route = useRoute()
    watch(() => [route.path, page.value.isNotFound],
          () => nextTick(() => redirect(!!page.value.isNotFound, site.value.base)),
          { immediate: true })
    return () => h(DefaultTheme.Layout)
  },
}

export default { extends: DefaultTheme, Layout }
