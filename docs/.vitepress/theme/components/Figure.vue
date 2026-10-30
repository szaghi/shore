<script setup lang="ts">
// A figure rendered in a light and a dark variant (docs/scripts/figures):
// public/figures/{light,dark}/<name>.<ext>, switched with the site theme.
import { withBase } from 'vitepress'
import { computed } from 'vue'

const props = withDefaults(
  defineProps<{ name: string; alt?: string; caption?: string; ext?: string; width?: string }>(),
  { alt: '', caption: '', ext: 'webp', width: '' },
)
const light = computed(() => withBase(`/figures/light/${props.name}.${props.ext}`))
const dark = computed(() => withBase(`/figures/dark/${props.name}.${props.ext}`))
const style = computed(() => (props.width ? { maxWidth: props.width } : {}))
</script>

<template>
  <figure class="shore-figure" :style="style">
    <img class="only-light" :src="light" :alt="alt || caption" loading="lazy" />
    <img class="only-dark" :src="dark" :alt="alt || caption" loading="lazy" />
    <figcaption v-if="caption || $slots.default">
      <slot>{{ caption }}</slot>
    </figcaption>
  </figure>
</template>
