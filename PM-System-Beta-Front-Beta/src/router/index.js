import { createRouter, createWebHashHistory } from 'vue-router'
import Forum from '../components/Forum.vue'
import MainLayout from '../layouts/MainLayout.vue'
import GlobalLayout from '../layouts/GlobalLayout.vue'
import Portal from '../views/Portal.vue'
import IpDisplay from '../components/IpDisplay.vue'

const routes = [
    {
        path: '/',
        redirect: '/portal'
    },
    {
        path: '/portal',
        component: Portal
    },
    {
        path: '/about',
        redirect: '/portal'
    },
    {
        path: '/analyse',
        component: MainLayout,
        children: [
            {
                path: '',
                component: () => import('../views/Platform.vue')
            },
            {
                path: 'Platform',
                component: () => import('../views/Platform.vue')
            },
            {
                path: 'UploadData',
                component: () => import('../views/UploadData.vue')
            }
        ]
    },
    {
        path: '/xenium-data',
        component: GlobalLayout,
        children: [
            {
                path: '',
                component: () => import('../views/XeniumData.vue')
            }
        ]
    },
    {
        path: '/ip',
        component: GlobalLayout,
        children: [
            {
                path: '',
                component: IpDisplay
            }
        ]
    }
]

const router = createRouter({
    history: createWebHashHistory(),
    routes
})

export default router
