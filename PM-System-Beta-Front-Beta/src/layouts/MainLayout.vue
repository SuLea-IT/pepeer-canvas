<template>
  <div>
    <Navbar />
    <el-container class="main-layout">
      <el-aside width="200px">
        <!-- 功能栏内容 -->
        <el-menu
          :default-active="activeMenu"
          @select="handleSelect"
          class="el-menu-vertical-demo"
        >
          <el-menu-item index="Platform" class="menu-item">
            <div class="team">
              <div>
                <svg
                  t="1723453629090"
                  class="icon"
                  viewBox="0 0 1028 1024"
                  version="1.1"
                  xmlns="http://www.w3.org/2000/svg"
                  p-id="4359"
                  width="200"
                  height="200"
                >
                  <path
                    d="M955.57056 783.104a24.512 24.512 0 1 1 0-48.96 23.232 23.232 0 0 0 19.136-9.856 25.088 25.088 0 0 0 2.944-22.72 315.008 315.008 0 0 0-196.864-191.552 24.448 24.448 0 1 1 15.232-46.528 364.16 364.16 0 0 1 227.712 221.568 74.24 74.24 0 0 1-9.088 67.456 71.36 71.36 0 0 1-59.072 30.592z"
                    p-id="4360"
                  ></path>
                  <path
                    d="M773.68256 513.92a24.512 24.512 0 0 1-11.84-45.888 184.32 184.32 0 0 0 95.488-161.408c0-89.472-63.872-165.76-151.744-181.312a24.576 24.576 0 0 1 8.512-48.256 232.832 232.832 0 0 1 192.192 229.568 233.408 233.408 0 0 1-120.96 204.288 23.616 23.616 0 0 1-11.648 3.008zM750.83456 832.064a24.448 24.448 0 1 1 0-48.96 32 32 0 0 0 26.368-13.632 33.6 33.6 0 0 0 4.032-31.104 377.536 377.536 0 0 0-236.096-229.76 24.576 24.576 0 0 1-15.68-30.912 24.512 24.512 0 0 1 30.848-15.616 426.688 426.688 0 0 1 266.944 259.84c9.152 25.472 5.376 53.824-10.176 75.84a80 80 0 0 1-66.24 34.304zM81.58656 832.064c-26.752 0-50.944-12.48-66.368-34.368a83.2 83.2 0 0 1-10.176-75.776 426.752 426.752 0 0 1 267.008-259.84 24.448 24.448 0 0 1 15.232 46.528 377.344 377.344 0 0 0-236.096 229.76 33.6 33.6 0 0 0 4.032 31.104 32 32 0 0 0 26.368 13.632 24.448 24.448 0 0 1 0 48.96z"
                    p-id="4361"
                  ></path>
                  <path
                    d="M416.17856 543.552a272.064 272.064 0 0 1-271.744-271.808A272.064 272.064 0 0 1 416.17856 0a272.064 272.064 0 0 1 271.808 271.744A272.128 272.128 0 0 1 416.17856 543.552z m0-494.592c-122.88 0-222.784 99.968-222.784 222.784S293.29856 494.592 416.17856 494.592s222.784-99.968 222.784-222.848S539.05856 48.96 416.17856 48.96z"
                    p-id="4362"
                  ></path>
                </svg>
              </div>
              <span>{{ $t("platformData") }}</span>
            </div></el-menu-item
          >
          <el-menu-item index="UploadData" class="menu-item">
            <div class="team">
              <div>
                <svg
                  t="1723453689216"
                  class="icon"
                  viewBox="0 0 1024 1024"
                  version="1.1"
                  xmlns="http://www.w3.org/2000/svg"
                  p-id="6257"
                  width="200"
                  height="200"
                >
                  <path
                    d="M768.2048 839.68c-14.848 0-26.9312-10.9568-26.9312-24.4736s12.0832-24.4736 26.9312-24.4736c111.3088 0 201.8304-96.5632 201.8304-215.3472s-90.5216-215.3472-201.8304-215.3472c-10.752 0-21.9136 1.1264-34.2016 3.3792-10.1376 1.8432-20.6848-1.7408-26.8288-9.3184-47.2064-57.1392-114.9952-90.0096-185.9584-90.0096-107.2128 0-201.8304 72.9088-235.2128 181.4528-1.9456 6.4512-6.7584 11.8784-13.2096 14.9504-6.4512 3.1744-14.0288 3.7888-21.0944 1.6384-13.824-4.096-27.8528-6.2464-41.5744-6.2464-86.1184 0-156.0576 75.0592-156.0576 167.424 0 92.2624 70.0416 167.424 156.0576 167.424h311.1936c14.848 0 26.9312 10.9568 26.9312 24.4736S536.064 839.68 521.1136 839.68H210.0224C94.208 839.68 0 742.6048 0 623.3088s94.208-216.3712 210.0224-216.3712c10.6496 0 21.2992 0.8192 31.9488 2.56C287.8464 292.2496 397.312 215.04 521.1136 215.04c83.0464 0 161.8944 35.4304 218.8288 97.6896 9.728-1.2288 19.0464-1.7408 28.2624-1.7408C909.312 310.9888 1024 429.568 1024 575.3856S909.312 839.68 768.2048 839.68z"
                    p-id="6258"
                  ></path>
                  <path
                    d="M663.2448 655.9744c-5.0176 0-10.1376-1.2288-14.6432-3.9936l-109.8752-64.8192-109.8752 64.8192c-12.4928 7.3728-29.184 4.096-37.2736-7.168-8.0896-11.3664-4.5056-26.5216 7.8848-33.8944l124.5184-73.5232c8.9088-5.2224 20.3776-5.2224 29.2864 0l124.6208 73.5232c12.4928 7.3728 16.0768 22.528 7.8848 33.8944-5.0176 7.168-13.7216 11.1616-22.528 11.1616z"
                    p-id="6259"
                  ></path>
                  <path
                    d="M761.1392 839.68C623.616 839.68 511.6928 725.7088 511.6928 585.5232c0-13.5168 12.0832-24.4736 26.9312-24.4736s26.9312 10.9568 26.9312 24.4736c0 113.152 87.6544 205.1072 195.4816 205.1072 14.848 0 26.9312 10.9568 26.9312 24.4736S775.9872 839.68 761.1392 839.68z"
                    p-id="6260"
                  ></path>
                </svg>
              </div>
              <span>{{ $t("uploadData") }}</span>
            </div></el-menu-item
          >
        </el-menu>
      </el-aside>

      <el-container>
        <el-main>
          <router-view />
        </el-main>
      </el-container>
    </el-container>
  </div>
</template>

<script setup>
import { ref, onMounted } from "vue";
import { useRouter, useRoute } from "vue-router";
import Navbar from "../components/Navbar.vue";

const activeMenu = ref("");
const router = useRouter();
const route = useRoute();

const handleSelect = (key) => {
  activeMenu.value = key;
   router.push({ path: `/analyse/${key}` });
};

onMounted(() => {
  const currentPath = route.path.split("/").pop();

  // 如果访问的是 /analyse 则默认高亮 Platform
  if (route.path === "/analyse" || !currentPath) {
    activeMenu.value = "Platform";
    router.push({ path: "/analyse/Platform" });
  } else {
     activeMenu.value = currentPath;
  }
});
</script>


<style scoped>
.main-layout {
  height: calc(100vh - 76px); /* 76px 是导航条的高度，根据实际情况调整 */
}
.el-container {
  overflow: visible;
}
.el-menu-vertical-demo {
  width: 100%;
  height: 100%;
  border-right: none;
  background-color: var(--el-container-bg-color);
}
.el-menu-vertical-demo::after {
  content: "";
  position: absolute;
  bottom: 0px;
  right: 0;
  width: 1px;
  background-color: #cccccc;
  transform: scaleX(1);
  transition: all 0.3s ease;
  height: 100%;
}
.el-main {
  padding: 0 !important;
  background-color: var(--el-container-bg-color);
  color: var(--el-navbar-color) !important;
}
.el-menu-vertical-demo :deep(.el-menu-item.is-active) {
  background: #dfdfdf;
}
.el-menu-vertical-demo :deep(.el-menu-item.is-active svg) {
  fill: var(--el-navbar-after);
}
.el-menu-vertical-demo :deep(.el-menu-item.is-active span) {
  color: var(--el-navbar-after);
}
.menu-item {
  padding-left: 0 !important;
}
.el-menu-vertical-demo :deep(.el-menu-item.is-active)::after {
  content: "";
  position: absolute;
  bottom: 0px;
  right: 0;
  width: 2px;
  background-color: #334c67;
  transform: scaleX(1);
  transition: all 0.3s ease;
  height: 100%;
}
.team {
  display: flex;
  flex-direction: column;
  align-items: center;
  width: 100%;
  height: 100%;
  justify-content: flex-end;
}
.team span {
  height: 26px;
  line-height: 24px;
}
.team svg {
  width: 24px;
  height: 24px;
  fill: var(--el-navbar-color);
}
</style>
