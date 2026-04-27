import { defineConfig } from 'vitepress'

export default defineConfig({
  title: 'SHORE',
  description: 'Structured Hyperbolic Overlapping Remeshing Engine — generate structured body-fitted O-grid volumes from raw triangulated surfaces.',
  base: '/shore/',
  markdown: {
    math: true,
  },
  themeConfig: {
    nav: [
      { text: 'Home', link: '/' },
      {
        text: 'Guide',
        items: [
          { text: 'About',        link: '/guide/' },
          { text: 'Installation', link: '/guide/installation' },
          { text: 'Quick Start',  link: '/guide/quickstart' },
          { text: 'Running',      link: '/guide/running' },
          { text: 'Changelog',    link: '/guide/changelog' },
          { text: 'Contributing', link: '/guide/contributing' },
        ],
      },
      { text: 'Algorithm',  link: '/algorithm/' },
      { text: 'Topologies', link: '/topologies/' },
      { text: 'Reference',  link: '/reference/mesh' },
      { text: 'GitHub',     link: 'https://github.com/szaghi/shore' },
    ],

    sidebar: {
      '/guide/': [
        {
          text: 'Introduction',
          items: [
            { text: 'About SHORE',  link: '/guide/' },
            { text: 'Installation', link: '/guide/installation' },
            { text: 'Quick Start',  link: '/guide/quickstart' },
          ],
        },
        {
          text: 'Project',
          items: [
            { text: 'Changelog',    link: '/guide/changelog' },
            { text: 'Contributing', link: '/guide/contributing' },
            { text: 'Running',      link: '/guide/running' },
          ],
        },
      ],

      '/concepts/': [
        {
          text: 'Concepts',
          items: [
            { text: 'Pipeline overview', link: '/concepts/pipeline' },
            { text: 'Chimera context',   link: '/concepts/chimera' },
          ],
        },
        {
          text: 'See also',
          items: [
            { text: 'Algorithm details',  link: '/algorithm/' },
            { text: 'Topology catalogue', link: '/topologies/' },
          ],
        },
      ],

      '/algorithm/': [
        {
          text: 'Algorithm',
          items: [
            { text: 'Overview',             link: '/algorithm/' },
            { text: 'Surface projection',   link: '/algorithm/projection' },
            { text: 'k=0 mesh construction',link: '/algorithm/k0-construction' },
            { text: 'Hyperbolic marching',  link: '/algorithm/hyperbolic-marching' },
            { text: 'Seam-aware normals',   link: '/algorithm/seam-aware-normals' },
            { text: 'Spacing laws',         link: '/algorithm/spacing-laws' },
          ],
        },
      ],

      '/topologies/': [
        {
          text: 'Topologies',
          items: [
            { text: 'Overview & taxonomy',    link: '/topologies/' },
            { text: 'O-grid (single-block)',  link: '/topologies/o-grid' },
            { text: 'Cubed sphere (6-block)', link: '/topologies/cubed-sphere' },
          ],
        },
        {
          text: 'Primitives',
          items: [
            { text: 'Box / annulus / flat-caps', link: '/topologies/primitives' },
          ],
        },
      ],

      '/reference/': [
        {
          text: 'Commands',
          items: [
            { text: 'shore mesh',      link: '/reference/mesh' },
            { text: 'shore primitive', link: '/reference/primitive' },
            { text: 'shore project',   link: '/reference/project' },
            { text: 'shore info',      link: '/reference/info' },
            { text: 'shore check',     link: '/reference/check' },
            { text: 'shore export',    link: '/reference/export' },
            { text: 'shore view',      link: '/reference/view' },
          ],
        },
        {
          text: 'Configuration',
          items: [
            { text: 'Config files', link: '/reference/config' },
          ],
        },
        {
          text: 'Python API',
          items: [
            { text: 'surface.io',             link: '/reference/api-surface-io' },
            { text: 'surface.remesh',         link: '/reference/api-surface-remesh' },
            { text: 'volume.hyperbolic',      link: '/reference/api-hyperbolic' },
            { text: 'volume.jacobian',        link: '/reference/api-jacobian' },
            { text: 'topology.sphere_cap',    link: '/reference/api-sphere-cap' },
            { text: 'io.geo',                 link: '/reference/api-geo' },
            { text: 'io.vtk',                 link: '/reference/api-vtk' },
            { text: 'vis',                    link: '/reference/api-vis' },
          ],
        },
      ],
    },

    socialLinks: [
      { icon: 'github', link: 'https://github.com/szaghi/shore' },
    ],

    search: {
      provider: 'local',
    },

    footer: {
      message: 'Released under the <a href="https://opensource.org/licenses/MIT">MIT License</a>.',
      copyright: 'Copyright © Stefano Zaghi',
    },
  },
})
