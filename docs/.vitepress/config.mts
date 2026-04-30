import { defineConfig } from 'vitepress'

export default defineConfig({
  title: 'SHORE',
  description: 'Structured Hyperbolic Overlapping Remeshing Engine — hyperbolic mesh generation for Chimera CFD: body-fitted near-body grids and parametric backgrounds.',
  base: '/shore/',
  markdown: {
    math: true,
  },
  themeConfig: {
    nav: [
      { text: 'Home',       link: '/' },
      { text: 'Guide',      link: '/guide/' },
      { text: 'Algorithm',  link: '/algorithm/' },
      { text: 'Topologies', link: '/topologies/' },
      { text: 'Reference',  link: '/reference/mesh' },
      { text: 'GitHub',     link: 'https://github.com/szaghi/shore' },
    ],

    sidebar: {
      '/guide/': [
        {
          text: 'Getting started',
          items: [
            { text: 'About SHORE',  link: '/guide/' },
            { text: 'Installation', link: '/guide/installation' },
            { text: 'Quick start',  link: '/guide/quickstart' },
          ],
        },
        {
          text: 'Concepts',
          items: [
            { text: 'Pipeline overview', link: '/guide/pipeline' },
            { text: 'Chimera context',   link: '/guide/chimera' },
          ],
        },
        {
          text: 'Worked examples',
          items: [
            { text: 'Assembly: sphere in cylinder', link: '/guide/assembly' },
            { text: 'Per-edge spacing',             link: '/guide/per-edge-spacing' },
            { text: 'cc.par + boxes pipeline',      link: '/guide/cc-par-pipeline' },
          ],
        },
        {
          text: 'Running examples and tests',
          link: '/guide/running',
        },
        {
          text: 'Project',
          items: [
            { text: 'Changelog',    link: '/guide/changelog' },
            { text: 'Contributing', link: '/guide/contributing' },
          ],
        },
      ],

      '/algorithm/': [
        {
          text: 'Pipeline',
          items: [
            { text: 'Overview',              link: '/algorithm/' },
            { text: 'Surface projection',    link: '/algorithm/projection' },
            { text: 'k=0 mesh construction', link: '/algorithm/k0-construction' },
            { text: 'Hyperbolic marching',   link: '/algorithm/hyperbolic-marching' },
            { text: 'Seam-aware normals',    link: '/algorithm/seam-aware-normals' },
          ],
        },
        {
          text: 'Distributions',
          items: [
            { text: 'Spacing laws', link: '/algorithm/spacing-laws' },
          ],
        },
        {
          text: 'Auxiliary algorithms',
          items: [
            { text: 'Marker boxes (voxel-fill + greedy merge)', link: '/algorithm/boxes' },
            { text: 'C-grid construction',                       link: '/algorithm/c-grid' },
            { text: 'Airfoil profiles (NACA 4 + 5-digit)',       link: '/algorithm/airfoil-profiles' },
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
            { text: 'C-grid (sharp-TE)',      link: '/topologies/c-grid' },
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
          text: 'Mesh generation commands',
          items: [
            { text: 'shore mesh',      link: '/reference/mesh' },
            { text: 'shore primitive', link: '/reference/primitive' },
            { text: 'shore project',   link: '/reference/project' },
          ],
        },
        {
          text: 'Body generation commands',
          items: [
            { text: 'shore profile', link: '/reference/profile' },
            { text: 'shore body',    link: '/reference/body' },
          ],
        },
        {
          text: 'Adjacency and parallelism',
          items: [
            { text: 'shore split',   link: '/reference/split' },
            { text: 'shore balance', link: '/reference/balance' },
          ],
        },
        {
          text: 'Xall pipeline',
          items: [
            { text: 'shore grd',          link: '/reference/grd' },
            { text: 'shore grd-merge',    link: '/reference/grd-merge' },
            { text: 'shore cc-par',       link: '/reference/cc-par' },
            { text: 'shore cc-par-merge', link: '/reference/cc-par-merge' },
            { text: 'shore proc-input',   link: '/reference/proc-input' },
            { text: 'Marker boxes',       link: '/reference/boxes' },
          ],
        },
        {
          text: 'Inspection and visualisation',
          items: [
            { text: 'shore info',   link: '/reference/info' },
            { text: 'shore check',  link: '/reference/check' },
            { text: 'shore export', link: '/reference/export' },
            { text: 'shore view',   link: '/reference/view' },
          ],
        },
        {
          text: 'Configuration',
          items: [
            { text: 'Config files', link: '/reference/config' },
          ],
        },
        {
          text: 'Python API — core',
          items: [
            { text: 'shore.mesh.Mesh',                link: '/reference/api-mesh' },
            { text: 'shore.mesh.topology',            link: '/reference/api-topology' },
            { text: 'shore.mesh.spacing (Spacing1D)', link: '/reference/api-spacing' },
            { text: 'shore.primitive',                link: '/reference/api-primitives' },
            { text: 'shore.split',                    link: '/reference/api-split' },
            { text: 'shore.balance',                  link: '/reference/api-balance' },
          ],
        },
        {
          text: 'Python API — surface',
          items: [
            { text: 'shore.surface.io',                  link: '/reference/api-surface-io' },
            { text: 'shore.surface.projector',           link: '/reference/api-projector' },
            { text: 'shore.surface.cap & .equator',      link: '/reference/api-cap-equator' },
            { text: 'shore.surface.remesh (legacy)',     link: '/reference/api-surface-remesh' },
          ],
        },
        {
          text: 'Python API — volume',
          items: [
            { text: 'shore.volume.spacing',    link: '/reference/api-volume-spacing' },
            { text: 'shore.volume.jacobian',   link: '/reference/api-jacobian' },
            { text: 'shore.volume.hyperbolic', link: '/reference/api-hyperbolic' },
          ],
        },
        {
          text: 'Python API — I/O',
          items: [
            { text: 'shore.io.geo',       link: '/reference/api-geo' },
            { text: 'shore.io.grd',       link: '/reference/api-grd' },
            { text: 'shore.io.vtk',       link: '/reference/api-vtk' },
            { text: 'shore.io.adjacency',  link: '/reference/api-adjacency' },
            { text: 'shore.io.cc_par',     link: '/reference/api-cc-par' },
            { text: 'shore.io.proc_input', link: '/reference/api-proc-input' },
            { text: 'shore.io.boxes',      link: '/reference/api-boxes-io' },
            { text: 'shore.boxes',         link: '/reference/api-boxes' },
          ],
        },
        {
          text: 'Python API — visualisation',
          items: [
            { text: 'shore.vis', link: '/reference/api-vis' },
          ],
        },
        {
          text: 'Python API — internals',
          items: [
            { text: 'shore.topology.sphere_cap', link: '/reference/api-sphere-cap' },
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
