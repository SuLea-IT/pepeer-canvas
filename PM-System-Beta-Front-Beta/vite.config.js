import { defineConfig } from 'vite'
import vue from '@vitejs/plugin-vue'

// https://vitejs.dev/config/
export default defineConfig({
  transpileDependencies: true,
  base: './',
  plugins: [vue()],
  server: {
    port: 5174
  }
})
