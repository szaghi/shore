import DefaultTheme from 'vitepress/theme'
import type { Theme } from 'vitepress'
import Figure from './components/Figure.vue'
import Mermaid from './components/Mermaid.vue'
import './custom.css'

export default {
  extends: DefaultTheme,
  enhanceApp({ app }) {
    app.component('Figure', Figure)
    app.component('Mermaid', Mermaid)
  },
} satisfies Theme
