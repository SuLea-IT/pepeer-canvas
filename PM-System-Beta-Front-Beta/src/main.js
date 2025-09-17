import { createApp } from 'vue'
import App from './App.vue'
import router from './router'
import ElementPlus from 'element-plus'
import 'element-plus/dist/index.css'
import * as ElementPlusIconsVue from '@element-plus/icons-vue'
import i18n from './i18n/i18n'
import './style.scss'
import './style.css'
import 'element-plus/theme-chalk/dark/css-vars.css'
import './styles/dark/css-vars.scss'
import { isDark } from './theme/composables/dark'

// 设置默认暗夜模式
isDark.value = true

// 创建 Vue 应用实例
const app = createApp(App)

// 注册 Element Plus
app.use(ElementPlus)

// 注册所有 Element Plus 图标组件
for (const [key, component] of Object.entries(ElementPlusIconsVue)) {
    app.component(key, component)
}

// 注册路由和 i18n
app.use(router)
app.use(i18n)

// 挂载应用
app.mount('#app')
