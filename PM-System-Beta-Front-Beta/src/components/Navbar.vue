<template>
  <nav class="navbar">
    <div class="navbar-left">
      <span class="navbar-title">{{ $t("title") }}</span>
    </div>
    <div class="navbar-center">
      <ul class="navbar-menu">
        <li>
          <router-link to="/about">{{ $t("about") }}</router-link>
        </li>
        <li>
          <router-link to="/analyse">{{ $t("analyse") }}</router-link>
        </li>
        <li>
          <router-link to="/xenium-data">{{ $t("xeniumData") }}</router-link>
        </li>
      </ul>
    </div>
    <div class="navbar-right">
      <el-popover
        v-model:visible="popoverVisible"
        placement="bottom"
        width="60"
        trigger="click"
      >
        <template #reference>
          <div class="language-selector">
            <svg
              class="language-icon"
              xmlns="http://www.w3.org/2000/svg"
              xmlns:xlink="http://www.w3.org/1999/xlink"
              aria-hidden="true"
              viewBox="0 0 24 24"
              width="1.2em"
              height="1.2em"
              data-v-63d067da=""
            >
              <path
                fill="currentColor"
                d="m18.5 10l4.4 11h-2.155l-1.201-3h-4.09l-1.199 3h-2.154L16.5 10h2zM10 2v2h6v2h-1.968a18.222 18.222 0 0 1-3.62 6.301a14.864 14.864 0 0 0 2.336 1.707l-.751 1.878A17.015 17.015 0 0 1 9 13.725a16.676 16.676 0 0 1-6.201 3.548l-.536-1.929a14.7 14.7 0 0 0 5.327-3.042A18.078 18.078 0 0 1 4.767 8h2.24A16.032 16.032 0 0 0 9 10.877a16.165 16.165 0 0 0 2.91-4.876L2 6V4h6V2h2zm7.5 10.885L16.253 16h2.492L17.5 12.885z"
              ></path>
            </svg>
          </div>
        </template>
        <div class="language-options">
          <div
            v-for="lang in languages"
            :key="lang.value"
            :class="{ 'active-language': localeLabel == lang.label }"
            @click="changeLanguage(lang.value)"
          >
            {{ lang.label }}
          </div>
        </div>
      </el-popover>
      <span class="navbar-button" @click="toggleDarkMode()">
        <Toggle :is-dark-mode="isDarkMode"></Toggle>
      </span>
    </div>
  </nav>
</template>

<script setup>
import { useI18n } from "vue-i18n";
import { ref, watchEffect } from "vue";
import { useDark, useToggle } from "@vueuse/core";
import Toggle from "../components/ToggleDark.vue";

const isDark = useDark();
const toggleDark = useToggle(isDark);
const isDarkMode = ref(isDark.value);

const toggleDarkMode = () => {
  toggleDark();
  isDarkMode.value = !isDarkMode.value;
};
// 获取当前语言和切换语言的方法
const { locale } = useI18n();

// 定义语言选项
const languages = [
  { label: "English", value: "en" },
  { label: "中文", value: "zh" },
];

// 获取当前语言的显示文本
const localeLabel = ref(
  languages.find((lang) => lang.value === locale.value)?.label
);

// 控制弹出框的显示状态
const popoverVisible = ref(false);

// 监听 locale 的变化，以便在渲染时正确显示当前语言
watchEffect(() => {
  localeLabel.value = languages.find(
    (lang) => lang.value === locale.value
  )?.label;
});

// 切换语言的方法
const changeLanguage = (newLocale) => {
  locale.value = newLocale;
  popoverVisible.value = false; // 选择完语言后关闭弹出框
};
</script>

<style scoped>
.navbar {
  width: 100%;
  display: flex;
  justify-content: space-between;
  align-items: center;
  background-color: var(--el-navbar-bg-color);
  color: var(--el-navbar-color);
  padding: 10px 20px;
  box-sizing: border-box;
  font-family: Arial, sans-serif;
  position: relative;
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
  z-index: 1000;
  top: 0;
  position: sticky;
}

.navbar-left {
  width: 20%;
}

.navbar-title {
  font-weight: bold;
  font-size: 16px;
}

.navbar-center {
  flex: 3;
  display: flex;
  justify-content: center;
  width: 60%;
}
.navbar-center a {
  color: var(--el-navbar-color);
}

.navbar-menu {
  list-style: none;
  display: flex;
  gap: 20px;
}
.navbar-menu li {
  position: relative;
}

.navbar-menu li a {
  text-decoration: none;
  font-weight: 500;
  padding: 10px 0;
  position: relative;
}

.navbar-menu li .router-link-active::after {
  content: "";
  position: absolute;
  bottom: -20px;
  left: 0;
  right: 0;
  height: 2px;
  background-color: var(--el-navbar-after);
  transform: scaleX(1);
  transition: all 0.3s ease;
  width: 100%;
}

.navbar-right {
  flex: 1;
  display: flex;
  justify-content: flex-end;
  align-items: center;
  gap: 10px;
  width: 20%;
}

.language-selector {
  display: flex;
  align-items: center;
  cursor: pointer;
  margin: 0 10px;
  box-sizing: border-box;
  margin-right: -20px;
}

.language-icon {
  width: 28px;
  height: 28px;
  margin-right: 8px;
}

.language-options {
  display: flex;
  flex-direction: column;
  gap: 10px;
}

.language-options div {
  cursor: pointer;
  padding: 5px 10px;
  border-radius: 4px;
}

.language-options div:hover {
  background-color: var(--el-navbar-hover-color);
  color: var(--el-navbar-color);
}

.active-language {
  font-weight: bold;
  color: #545454;
  background-color: #b6e8ff; /* 突出显示当前语言 */
  border-radius: 4px;
}
</style>
