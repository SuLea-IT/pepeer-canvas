<template>
  <div class="portal-page">
    <!-- 顶部导航栏 -->
    <header class="header-bg sticky top-0 z-50 shadow-md">
      <div class="container mx-auto px-4">
        <div class="flex justify-between items-center py-3">
          <div class="flex items-center space-x-4">
            <div
              class="w-12 h-12 bg-white rounded-full flex items-center justify-center"
            >
              <span class="text-primary font-bold text-xl">C</span>
            </div>
            <div>
              <h1 class="text-xl font-bold text-white">CaSpatialDB</h1>
            </div>
          </div>

          <!-- 桌面导航 -->
          <nav class="hidden md:flex items-center space-x-1">
            <a href="#" class="nav-item relative">Home</a>
            <div class="relative group">
              <a href="#" class="nav-item flex items-center">
                <span>Data</span>
                <i class="fa fa-angle-down ml-1 text-xs"></i>
              </a>
              <div
                class="absolute left-0 mt-2 w-48 bg-white rounded-md shadow-lg opacity-0 invisible group-hover:opacity-100 group-hover:visible transition-all duration-300 z-50"
              >
                <a
                  href="https://ngdc.cncb.ac.cn/gsa/s/65g8iKs4"
                  target="_blank"
                  class="dropdown-item"
                  >Single-cell Data</a
                >
                <a
                  href="https://ngdc.cncb.ac.cn/gsa/s/qUIIG8Z9"
                  target="_blank"
                  class="dropdown-item"
                  >Spatial Transcriptomics Data</a
                >
              </div>
            </div>
            <a href="#" class="nav-item relative">Documentation</a>
          </nav>

          <!-- 移动导航 -->
          <div class="md:hidden">
            <button
              @click="toggleMobileMenu"
              class="text-white focus:outline-none"
            >
              <i class="fa fa-bars"></i>
            </button>
            <div
              v-show="showMobileMenu"
              class="absolute top-full left-0 right-0 bg-primary p-4"
            >
              <a
                href="#"
                class="block py-2 px-3 text-white hover:bg-primary-700 rounded"
                >Home</a
              >
              <div class="mobile-dropdown">
                <button
                  @click="toggleMobileSubmenu"
                  class="flex justify-between items-center w-full py-2 px-3 text-white hover:bg-primary-700 rounded"
                >
                  Data
                  <i class="fa fa-angle-down"></i>
                </button>
                <div v-show="showMobileSubmenu" class="pl-4 space-y-1 py-2">
                  <a
                    href="https://ngdc.cncb.ac.cn/gsa/s/65g8iKs4"
                    target="_blank"
                    class="block py-1 px-3 text-white hover:bg-primary-700 rounded"
                    >Single-cell Data</a
                  >
                  <a
                    href="https://ngdc.cncb.ac.cn/gsa/s/qUIIG8Z9"
                    target="_blank"
                    class="block py-1 px-3 text-white hover:bg-primary-700 rounded"
                    >Spatial Transcriptomics Data</a
                  >
                </div>
              </div>
              <a
                href="#"
                class="py-2 px-3 text-white hover:bg-primary-700 rounded"
                >Documentation</a
              >
            </div>
          </div>
        </div>
      </div>
    </header>

    <!-- 英雄区域 -->
    <section class="relative h-[400px] overflow-hidden">
      <img :src="bannerImage" alt="Chili pepper plant" class="hero-image" />
      <div class="absolute inset-0 bg-black/40 z-10"></div>
      <div class="absolute inset-0 z-20 flex items-center">
        <div class="container mx-auto px-4">
          <div class="max-w-4xl">
            <h2
              class="text-[clamp(1.5rem,3vw,2.5rem)] font-bold text-white leading-tight mb-4"
            >
              CaSpatialDB
            </h2>
            <p class="text-[clamp(1rem,2vw,1.25rem)] text-white/90 mb-8">
              A Spatiotemporal Single-Cell Atlas of Chili Pepper Development
            </p>
          </div>
        </div>
      </div>
    </section>

    <!-- 辣椒发育阶段与模式图 -->
    <section class="py-12 bg-white">
      <div class="container mx-auto px-4">
        <div class="flex flex-col md:flex-row items-center gap-8">
          <div class="w-full md:w-1/2">
            <img
              :src="figure1Image"
              alt="Chili pepper development stages"
              class="w-full h-auto rounded-lg shadow-lg"
            />
          </div>

          <div class="w-full md:w-1/2">
            <h2 class="text-2xl font-bold text-primary mb-4">Abstract</h2>
            <p class="text-gray-700 mb-4 leading-relaxed">
              Chili peppers (Capsicum annuum L.) are globally valued for their
              spiciness (capsaicin) and vivid red color (capsanthin). Although
              the biochemical pathways of these specialized metabolites have
              been partially elucidated, how their production is regulated at
              the single-cell level within the developing fruit remains poorly
              understood.
            </p>
            <p class="text-gray-700 mb-4 leading-relaxed">
              Here, we integrated comprehensive single-nucleus RNA sequencing
              with high-resolution spatial transcriptomics to construct the
              first spatiotemporal single-cell atlas of pepper development,
              encompassing 29 samples and nearly 200,000 cells from seedlings,
              roots, stems, leaves, flowers, and fruits. Our analyses revealed a
              distinct multilayered cellular architecture in multiple tissues.
            </p>
            <p class="text-gray-700 leading-relaxed">
              In particular, we identified a novel class of Laminar Patterning
              Transcription Factors (LPTFs) characterized by layer-specific
              expression and roles in determining cell fate and spatially
              restricting specialized metabolism. Among these, WRKY6 governs
              capsanthin biosynthesis in pericarp layers, while WRKY9 controls
              capsaicin production within distinct placental layers.
            </p>
          </div>
        </div>
      </div>
    </section>

    <!-- 统计信息方框 -->
    <section class="py-12 bg-gray-50">
      <div class="mx-auto">
        <div class="feature-card p-6 text-center rounded-none">
          <p class="text-xl md:text-2xl font-bold">
            <span class="text-primary">Single-nucleus Datasets: 14</span>
            &nbsp;&nbsp;
            <span class="text-primary"
              >Spatial Transcriptomics Datasets: 39</span
            >
            &nbsp;&nbsp;
            <span class="text-primary">Total Cells: 263,076</span>
          </p>
        </div>
      </div>
    </section>

    <!-- 基因表达和细胞类型模块 -->
    <section class="py-12 bg-gray-50">
      <div class="container mx-auto px-4">
        <div class="grid grid-cols-1 md:grid-cols-2 gap-6">
          <!-- 基因表达查询 -->
          <div class="feature-card flex flex-col h-full">
            <div class="md:flex flex-row-reverse flex-grow">
              <div class="md:w-1/2 feature-image-container h-[300px] md:h-auto">
                <img
                  :src="geneImage"
                  alt="Gene expression patterns"
                  class="feature-image"
                />
              </div>
              <div
                class="md:w-1/2 feature-content flex flex-col justify-center"
              >
                <h3 class="text-xl font-bold text-primary mb-4 text-center">
                  Gene Expression
                </h3>
                <p class="text-gray-700 mb-6 text-center">
                  Explore gene expression patterns across different
                  developmental stages and tissues in chili pepper to understand
                  gene function and regulatory networks.
                </p>
              </div>
            </div>
            <div class="p-4 text-center">
              <button
                @click="startGeneQuery"
                class="bg-secondary text-white px-6 py-3 rounded-lg text-lg hover:bg-orange-600 transition-colors"
                style="border:none"
              >
                <i class="fa fa-search"></i> Start Query
              </button>
            </div>
          </div>

          <!-- 细胞类型查询 -->
          <div class="feature-card flex flex-col h-full">
            <div class="md:flex flex-row-reverse flex-grow">
              <div class="md:w-1/2 feature-image-container h-[300px] md:h-auto">
                <img
                  :src="celltypeImage"
                  alt="Cell type distribution"
                  class="feature-image"
                />
              </div>
              <div
                class="md:w-1/2 feature-content flex flex-col justify-center"
              >
                <h3 class="text-xl font-bold text-primary mb-4 text-center">
                  Cell Types
                </h3>
                <p class="text-gray-700 mb-6 text-center">
                  Explore cell type distributions across different developmental
                  stages and tissues in chili pepper.
                </p>
              </div>
            </div>
            <div class="p-4 text-center">
              <button
                @click="startCellTypeQuery"
                class="bg-secondary text-white px-6 py-3 rounded-lg text-lg hover:bg-orange-600 transition-colors"
                style="border:none"
              >
                <i class="fa fa-search"></i> Start Query
              </button>
            </div>
          </div>
        </div>
      </div>
    </section>

    <!-- 引用部分 -->
    <section class="py-12 bg-gray-50">
      <div class="container mx-auto px-4">
        <div class="max-w-4xl mx-auto bg-white p-8 rounded-lg shadow-md">
          <div class="text-lg font-bold text-primary mb-4">Citation:</div>
          <div class="text-gray-700">
            {{ citationContent }}
          </div>
        </div>
      </div>
    </section>

    <!-- 页脚 -->
    <footer class="bg-primary text-white py-10">
      <div class="container mx-auto px-4">
        <div class="grid grid-cols-1 md:grid-cols-1 gap-8">
          <div class="text-center">
            <div class="flex justify-center mb-4">
              <div
                class="w-16 h-16 bg-white rounded-full flex items-center justify-center mx-auto"
              >
                <span class="text-primary font-bold text-2xl">C</span>
              </div>
            </div>
            <h4 class="text-xl font-bold mb-4">CaSpatialDB</h4>
            <p class="text-white/70 mb-6 max-w-2xl mx-auto">
              A Spatiotemporal Single-Cell Atlas of Chili Pepper Development
            </p>
            <p class="text-white/70 mb-2">Contact: example@example.com</p>
          </div>
        </div>

        <div class="text-center text-white/70 mt-8">
          Copyright &copy; 2025 All rights reserved | Best viewed with Chrome
          and Edge browsers
        </div>
      </div>
    </footer>
  </div>
</template>

<script setup>
import { ref } from "vue";
import { useRouter } from "vue-router";

const router = useRouter();

// 响应式数据
const showMobileMenu = ref(false);
const showMobileSubmenu = ref(false);
const citationContent = ref("");

// 图片路径（从assets目录导入）
import bannerImg from "../assets/bg.PNG";
import figure1Img from "../assets/img1.PNG";
import geneImg from "../assets/img2.PNG";
import celltypeImg from "../assets/umap.PNG";

const bannerImage = ref(bannerImg);
const figure1Image = ref(figure1Img);
const geneImage = ref(geneImg);
const celltypeImage = ref(celltypeImg);

// 方法
const toggleMobileMenu = () => {
  showMobileMenu.value = !showMobileMenu.value;
};

const toggleMobileSubmenu = () => {
  showMobileSubmenu.value = !showMobileSubmenu.value;
};

const startGeneQuery = () => {
  // 跳转到分析页面
  router.push('/analyse');
};

const startCellTypeQuery = () => {
  // 跳转到分析页面
  router.push('/analyse');
};
</script>

<style scoped>
/* 全局样式重置 */
.portal-page {
  background-color: #f9f9f9;
  min-height: 100vh;
}

/* 导航栏样式 */
.header-bg {
  background-image: linear-gradient(to right, #8b0000, #b22222);
}

.nav-item {
  @apply px-3 py-2 text-white hover:text-yellow-400 transition-colors;
  text-decoration: none !important;
  color: white !important;
}

.nav-item:hover {
  color: #fbbf24 !important;
  text-decoration: none !important;
}

.dropdown-item {
  @apply block px-4 py-2 text-sm text-gray-700 hover:bg-gray-100;
  text-decoration: none !important;
  color: #374151 !important;
}

.dropdown-item:hover {
  background-color: #f3f4f6 !important;
  color: #374151 !important;
  text-decoration: none !important;
}

/* 移除所有链接的默认样式 */
.portal-page a {
  text-decoration: none !important;
}

.portal-page a:hover {
  text-decoration: none !important;
}

/* 英雄区域样式 */
.hero-image {
  @apply w-full h-full object-contain;
  max-height: 400px;
  background-color: #000;
}

section.relative {
  background-color: #000;
  min-height: 400px;
  display: flex;
  align-items: center;
  justify-content: center;
}

/* 特色卡片样式 - 确保白底和投影 */
.feature-card {
  background-color: white !important;
  border-radius: 12px;
  box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06);
  overflow: hidden;
  transition: all 0.3s ease;
  height: 100%;
}

.feature-card:hover {
  box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05);
  transform: translateY(-2px);
}

.feature-content {
  @apply p-4 md:p-6;
}

.feature-image-container {
  @apply w-full overflow-hidden;
}

.feature-image {
  @apply w-full h-full object-cover transition-transform duration-500 hover:scale-105;
}

/* 颜色变量 */
.text-primary {
  color: #8b0000;
}

.bg-primary {
  background-color: #8b0000;
}

.bg-secondary {
  background-color: #ff8c00;
}

.text-accent {
  color: #ffff00;
}

/* 容器样式 */
.container {
  max-width: 1200px;
}

/* 背景色调整 */
section.bg-white {
  background-color: #fff;
}

section.bg-gray-50 {
  background-color: #f9f9f9;
}

/* 按钮样式增强 */
.action-button {
  transition: all 0.3s ease;
  border: none !important;
}

.action-button:hover {
  transform: translateY(-2px);
  box-shadow: 0 4px 12px rgba(255, 140, 0, 0.3);
}

/* 移动端适配 */
@media (max-width: 768px) {
  .hero-image {
    height: 300px;
    object-fit: contain;
  }
  
  .feature-card {
    margin-bottom: 1rem;
  }
  
  section.relative {
    min-height: 300px;
  }
}
</style>















