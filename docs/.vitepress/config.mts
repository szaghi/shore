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
      { text: 'Concepts', link: '/concepts/pipeline' },
      { text: 'Reference', link: '/reference/mesh' },
      { text: 'GitHub', link: 'https://github.com/szaghi/shore' },
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
            { text: 'Pipeline overview',      link: '/concepts/pipeline' },
            { text: 'Hyperbolic marching',     link: '/concepts/hyperbolic-marching' },
            { text: 'Spacing laws',            link: '/concepts/spacing' },
            { text: 'O-grid topology',         link: '/concepts/o-grid' },
            { text: 'Sphere polar cap',        link: '/concepts/sphere-cap' },
            { text: 'Chimera mesh generation', link: '/concepts/chimera' },
          ],
        },
      ],

      '/reference/': [
        {
          text: 'Commands',
          items: [
            { text: 'shore mesh',    link: '/reference/mesh' },
            { text: 'shore project', link: '/reference/project' },
            { text: 'shore march',   link: '/reference/march' },
            { text: 'shore info',    link: '/reference/info' },
            { text: 'shore check',   link: '/reference/check' },
            { text: 'shore export',  link: '/reference/export' },
            { text: 'shore view',    link: '/reference/view' },
          ],
        },
        {
          text: 'Python API',
          items: [
            { text: 'surface.io',             link: '/reference/api-surface-io' },
            { text: 'surface.remesh',         link: '/reference/api-surface-remesh' },
            { text: 'volume.hyperbolic',      link: '/reference/api-hyperbolic' },
            { text: 'volume.jacobian',        link: '/reference/api-jacobian' },
            { text: 'volume.spacing',         link: '/reference/api-spacing' },
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
