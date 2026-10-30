<script setup lang="ts">
// Client-side Mermaid, drawn in the site palette and redrawn on theme change.
import { useData } from 'vitepress'
import { onMounted, ref, watch } from 'vue'

const props = defineProps<{ code: string }>()
const { isDark } = useData()
const svg = ref('')
let count = 0

const palette = (dark: boolean) => ({
  fontFamily: 'Inter, ui-sans-serif, system-ui, sans-serif',
  fontSize: '14px',
  background: 'transparent',
  primaryColor: dark ? '#173a47' : '#e3f1f5',
  primaryBorderColor: dark ? '#4aa3c2' : '#1f6f8b',
  primaryTextColor: dark ? '#e3e8ee' : '#1d2733',
  secondaryColor: dark ? '#1c3d38' : '#e2f3f0',
  secondaryBorderColor: dark ? '#43c2b3' : '#2a9d8f',
  tertiaryColor: dark ? '#3d3420' : '#fbf2db',
  tertiaryBorderColor: dark ? '#f0d283' : '#c9a227',
  lineColor: dark ? '#9aa8b6' : '#5b6b7b',
  textColor: dark ? '#e3e8ee' : '#1d2733',
  clusterBkg: dark ? '#20262d' : '#f5f8fa',
  clusterBorder: dark ? '#33404d' : '#d5dde5',
  edgeLabelBackground: dark ? '#1b1b1f' : '#ffffff',
  noteBkgColor: dark ? '#3d3420' : '#fbf2db',
  noteBorderColor: dark ? '#f0d283' : '#c9a227',
})

async function render() {
  const mermaid = (await import('mermaid')).default
  mermaid.initialize({
    startOnLoad: false,
    theme: 'base',
    themeVariables: palette(isDark.value),
  })
  const id = `mermaid-${Math.random().toString(36).slice(2)}-${count++}`
  svg.value = (await mermaid.render(id, decodeURIComponent(props.code))).svg
}

onMounted(render)
watch(isDark, render)
</script>

<template>
  <div class="shore-mermaid" v-html="svg" />
</template>
