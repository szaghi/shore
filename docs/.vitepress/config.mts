import { readdirSync, readFileSync } from 'node:fs'
import { fileURLToPath } from 'node:url'
import { defineConfig } from 'vitepress'

const docs = fileURLToPath(new URL('..', import.meta.url))
const apiSidebar = JSON.parse(readFileSync(`${docs}/reference/api/_sidebar.json`, 'utf8'))

// Every ADR, titled from its first heading: "0042 Blocking: ...".
function adrItems() {
  return readdirSync(`${docs}/adr`)
    .filter((f) => /^\d{4}-.*\.md$/.test(f))
    .sort()
    .map((f) => {
      const head = readFileSync(`${docs}/adr/${f}`, 'utf8').match(/^#\s+(.*)$/m)?.[1] ?? f
      const text = head.replace(/^(\d{4})\s*[—-]\s*/, '$1 ').replace(/`/g, '')
      return { text, link: `/adr/${f.replace(/\.md$/, '')}` }
    })
}

function adrGroups() {
  const all = adrItems()
  const phase = (lo: number, hi: number) =>
    all.filter((i) => {
      const n = Number(i.text.slice(0, 4))
      return n >= lo && n <= hi
    })
  return [
    { text: 'Design', items: [{ text: 'Index', link: '/adr/' }, { text: 'Roadmap', link: '/adr/roadmap' }] },
    { text: 'Foundations (0001–0006)', collapsed: true, items: phase(1, 6) },
    { text: 'Surface grids, collars (0007–0032)', collapsed: true, items: phase(7, 32) },
    { text: 'Projects, exports (0033–0036)', collapsed: true, items: phase(33, 36) },
    { text: 'GUI (0037–0041)', collapsed: true, items: phase(37, 41) },
    { text: 'Interactive blocking (0042–)', collapsed: true, items: phase(42, 9999) },
  ]
}

const cli = [
  {
    text: 'Grids',
    items: [
      { text: 'shore mesh', link: '/reference/mesh' },
      { text: 'shore primitive', link: '/reference/primitive' },
      { text: 'shore collar', link: '/reference/collar' },
      { text: 'shore cap-tip', link: '/reference/cap-tip' },
      { text: 'shore project', link: '/reference/project' },
      { text: 'shore run', link: '/reference/run' },
      { text: 'shore gui', link: '/reference/gui' },
    ],
  },
  {
    text: 'Bodies',
    items: [
      { text: 'shore profile', link: '/reference/profile' },
      { text: 'shore body', link: '/reference/body' },
    ],
  },
  {
    text: 'Xall and chimera',
    items: [
      { text: 'shore xall', link: '/reference/xall' },
      { text: 'shore grd', link: '/reference/grd' },
      { text: 'shore grd-merge', link: '/reference/grd-merge' },
      { text: 'shore cc-par', link: '/reference/cc-par' },
      { text: 'shore cc-par-merge', link: '/reference/cc-par-merge' },
      { text: 'shore proc-input', link: '/reference/proc-input' },
      { text: 'Marker boxes', link: '/reference/boxes' },
      { text: 'shore split', link: '/reference/split' },
      { text: 'shore balance', link: '/reference/balance' },
    ],
  },
  {
    text: 'Exports and inspection',
    items: [
      { text: 'shore cgns', link: '/reference/cgns' },
      { text: 'shore plot3d', link: '/reference/plot3d' },
      { text: 'shore export', link: '/reference/export' },
      { text: 'shore info', link: '/reference/info' },
      { text: 'shore check', link: '/reference/check' },
      { text: 'shore view', link: '/reference/view' },
    ],
  },
  {
    text: 'Files',
    items: [
      { text: 'Config files', link: '/reference/config' },
      { text: 'File formats', link: '/reference/formats' },
    ],
  },
]

export default defineConfig({
  title: 'SHORE',
  description:
    'Structured multi-block grids for overset (Chimera) CFD: hyperbolic near-body grids, collars, interactive blocking on STL and CAD, and the Xall solver bundle.',
  base: '/shore/',
  cleanUrls: true,
  rewrites: { 'adr/README.md': 'adr/index.md' },
  lastUpdated: true,
  head: [['link', { rel: 'icon', type: 'image/svg+xml', href: '/shore/logo.svg' }]],
  markdown: {
    math: true,
    config(md) {
      const fence = md.renderer.rules.fence!
      md.renderer.rules.fence = (tokens, idx, options, env, self) => {
        const token = tokens[idx]
        if (token.info.trim() === 'mermaid') {
          return `<Mermaid code="${encodeURIComponent(token.content)}" />`
        }
        return fence(tokens, idx, options, env, self)
      }
    },
  },
  themeConfig: {
    logo: '/logo.svg',
    nav: [
      { text: 'Get started', link: '/guide/', activeMatch: '/guide/' },
      { text: 'Tutorials', link: '/tutorials/', activeMatch: '/tutorials/' },
      { text: 'Cookbook', link: '/cookbook/', activeMatch: '/cookbook/' },
      { text: 'Concepts', link: '/concepts/', activeMatch: '/concepts/' },
      {
        text: 'Reference',
        activeMatch: '/reference/',
        items: [
          { text: 'CLI commands', link: '/reference/' },
          { text: 'Python API', link: '/reference/api/' },
          { text: 'Project scripts', link: '/reference/run' },
          { text: 'GUI', link: '/reference/gui' },
          { text: 'Config files', link: '/reference/config' },
          { text: 'File formats', link: '/reference/formats' },
        ],
      },
      { text: 'Design', link: '/adr/', activeMatch: '/adr/' },
    ],

    sidebar: {
      '/guide/': [
        {
          text: 'Get started',
          items: [
            { text: 'What is SHORE?', link: '/guide/' },
            { text: 'Installation', link: '/guide/installation' },
            { text: 'Quick start', link: '/guide/quickstart' },
            { text: 'Choosing a workflow', link: '/guide/workflows' },
          ],
        },
        {
          text: 'Project',
          items: [
            { text: 'Examples and tests', link: '/guide/running' },
            { text: 'Contributing', link: '/guide/contributing' },
            { text: 'Changelog', link: '/guide/changelog' },
          ],
        },
      ],

      '/tutorials/': [
        {
          text: 'Tutorials',
          items: [
            { text: 'Overview', link: '/tutorials/' },
            { text: '1. A sphere, end to end', link: '/tutorials/sphere' },
            { text: '2. A wing with a C-grid', link: '/tutorials/wing' },
            { text: '3. A Chimera assembly', link: '/tutorials/chimera-assembly' },
            { text: '4. The Xall bundle', link: '/tutorials/xall-bundle' },
            { text: '5. Collar at a junction', link: '/tutorials/collar' },
            { text: '6. Projects and replay', link: '/tutorials/projects' },
            { text: '7. Blocking a wing–body', link: '/tutorials/blocking' },
            { text: '8. The same in the GUI', link: '/tutorials/gui' },
          ],
        },
      ],

      '/cookbook/': [
        {
          text: 'Cookbook',
          items: [{ text: 'All recipes', link: '/cookbook/' }],
        },
        {
          text: 'Geometry',
          items: [
            { text: 'Make a test body', link: '/cookbook/test-bodies' },
            { text: 'Load STL and CAD', link: '/cookbook/load-geometry' },
            { text: 'Feature curves of an STL', link: '/cookbook/stl-curves' },
          ],
        },
        {
          text: 'Grids',
          items: [
            { text: 'Wall-normal spacing', link: '/cookbook/wall-spacing' },
            { text: 'Per-edge spacing', link: '/cookbook/per-edge-spacing' },
            { text: 'Redistribute an edge', link: '/cookbook/redistribute' },
            { text: 'Rebuild a face on a CAD surface', link: '/cookbook/cad-face' },
            { text: 'Associate with geometry', link: '/cookbook/associate' },
            { text: 'Half models and symmetry', link: '/cookbook/half-models' },
            { text: 'Check grid quality', link: '/cookbook/quality' },
          ],
        },
        {
          text: 'Blocking',
          items: [
            { text: 'Blocking recipes', link: '/cookbook/blocking' },
          ],
        },
        {
          text: 'Solver files',
          items: [
            { text: 'Export CGNS and Plot3D', link: '/cookbook/exports' },
            { text: 'Multi-component Xall run', link: '/cookbook/multi-component' },
            { text: 'Split and balance', link: '/cookbook/split-balance' },
            { text: 'After-overset splitting', link: '/cookbook/after-overset-split' },
          ],
        },
        {
          text: 'Projects',
          items: [
            { text: 'Save, load, replay', link: '/cookbook/projects' },
          ],
        },
      ],

      '/concepts/': [
        {
          text: 'Concepts',
          items: [
            { text: 'How SHORE works', link: '/concepts/' },
            { text: 'Block model', link: '/concepts/block-model' },
            { text: 'Geometry', link: '/concepts/geometry' },
            { text: 'Association', link: '/concepts/association' },
            { text: 'Projects and journals', link: '/concepts/projects' },
            { text: 'Grid quality', link: '/concepts/quality' },
          ],
        },
        {
          text: 'Near-body marching',
          collapsed: false,
          items: [
            { text: 'Pipeline', link: '/concepts/near-body-pipeline' },
            { text: 'Surface projection', link: '/concepts/projection' },
            { text: 'k=0 construction', link: '/concepts/k0-construction' },
            { text: 'Marching schemes', link: '/concepts/hyperbolic-marching' },
            { text: 'Seam-aware normals', link: '/concepts/seam-aware-normals' },
            { text: 'Spacing laws', link: '/concepts/spacing-laws' },
            { text: 'Limits', link: '/concepts/pipeline-limits' },
          ],
        },
        {
          text: 'Topologies',
          collapsed: true,
          items: [
            { text: 'Overview', link: '/concepts/topologies/' },
            { text: 'O-grid', link: '/concepts/topologies/o-grid' },
            { text: 'Cubed sphere', link: '/concepts/topologies/cubed-sphere' },
            { text: 'C-grid', link: '/concepts/topologies/c-grid' },
            { text: 'H-grid', link: '/concepts/topologies/h-grid' },
            { text: 'OH hybrid', link: '/concepts/topologies/oh-hybrid' },
            { text: 'CH hybrid', link: '/concepts/topologies/ch-hybrid' },
            { text: 'Primitives', link: '/concepts/topologies/primitives' },
            { text: 'C-grid construction', link: '/concepts/c-grid-construction' },
            { text: 'Airfoil profiles', link: '/concepts/airfoil-profiles' },
          ],
        },
        {
          text: 'Junctions',
          collapsed: true,
          items: [
            { text: 'Surface grids', link: '/concepts/surface-grids' },
            { text: 'Intersection curves', link: '/concepts/intersection-curves' },
            { text: 'Surface marching', link: '/concepts/surface-marching' },
            { text: 'Collar blocks', link: '/concepts/collars' },
          ],
        },
        {
          text: 'Interactive blocking',
          collapsed: true,
          items: [{ text: 'Blocking', link: '/concepts/blocking' }],
        },
        {
          text: 'Chimera and Xall',
          collapsed: true,
          items: [
            { text: 'Chimera grids', link: '/concepts/chimera' },
            { text: 'Marker boxes', link: '/concepts/marker-boxes' },
            { text: 'Load balancing', link: '/concepts/load-balancing' },
          ],
        },
      ],

      '/reference/api/': [
        { text: 'Python API', items: [{ text: 'Index', link: '/reference/api/' }] },
        ...apiSidebar,
      ],

      '/reference/': [{ text: 'CLI', items: [{ text: 'All commands', link: '/reference/' }] }, ...cli],

      '/adr/': adrGroups(),
    },

    socialLinks: [{ icon: 'github', link: 'https://github.com/szaghi/shore' }],
    editLink: {
      pattern: 'https://github.com/szaghi/shore/edit/main/docs/:path',
      text: 'Edit this page on GitHub',
    },
    outline: { level: [2, 3] },
    search: { provider: 'local' },
    footer: {
      message: 'Released under the <a href="https://opensource.org/licenses/MIT">MIT License</a>.',
      copyright: 'Copyright © Stefano Zaghi',
    },
  },
  vite: {
    build: { target: 'es2022', chunkSizeWarningLimit: 4000 },
  },
})
