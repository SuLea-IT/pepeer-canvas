<template>
  <el-upload
    class="upload-demo"
    drag
    action="https://run.mocky.io/v3/9d059bf9-4660-45f2-925d-ce80ad6c4d15"
    multiple
    :accept="allowedExtensions"
    :before-upload="beforeUpload"
    :on-change="handleChange"
    :auto-upload="false"
    :show-file-list="showFileList"
    :upload-flag="uploadFlag"
  >
    <el-icon class="el-icon--upload"><upload-filled /></el-icon>
    <div class="el-upload__text" v-html="$t('dropFileText')"></div>
    <template #tip>
      <div class="el-upload__tip">
        {{ $t("fileSizeTip") }}
      </div>
    </template>
  </el-upload>
</template>

<script setup lang="ts">
import { UploadFilled } from "@element-plus/icons-vue";
import { ElMessage } from "element-plus";
import i18n from "../../i18n/i18n";
import { ref, watch } from "vue";

const props = defineProps({
  uploadFlag: Boolean,
  modelValue: Array, // v-model绑定
  selectedFileType: Number, // 父组件传递的selectedFileType
  selectedFun: Number, // 父组件传递的selectedFun
  fileTypeCluster: Object, // 父组件传递的fileTypeCluster对象
  fileTypeGene: Object, // 父组件传递的fileTypeGene对象
  showFileList: Boolean, // 是否显示文件列表
});

const emit = defineEmits(["update:modelValue"]); // v-model emit

const allowedExtensions = ref("");

// 监听 selectedFileType 和 selectedFun 的变化，更新允许的扩展名
watch(
  () => [props.selectedFileType, props.selectedFun],
  ([newFileType, newFun]) => {
    const restrictions =
      newFun === 1
        ? props.fileTypeCluster[newFileType]
        : props.fileTypeGene[newFileType];
    allowedExtensions.value = restrictions?.allowedExtensions.join(",") || "";
    console.log("Updated allowedExtensions:", allowedExtensions.value);
  },
  { immediate: true }
);

// 上传前的文件验证
const beforeUpload = (file) => {
  const isValid = validateFile(file, props.selectedFileType, props.selectedFun);
  if (!isValid) {
    ElMessage.error(i18n.global.t("fileUploadError", { type: file.name }));
  }
  return isValid;
};

// 处理文件变化，保留并累加文件
const handleChange = (file, fileList) => {
  if (props.uploadFlag) {
    console.log("子组件执行清空");

    fileList = [];
  }
  const isValid = validateFile(file, props.selectedFileType, props.selectedFun);
  if (isValid) {
    const newFiles = fileList.map((f) => f.raw); // 获取新选择的文件
    emit("update:modelValue", newFiles); // 覆盖现有文件列表，而不是累加
  }
};

// 验证文件类型和数量
const validateFile = (file, fileType, funType) => {
  const restrictions =
    funType === 1
      ? props.fileTypeCluster[fileType]
      : props.fileTypeGene[fileType];

  const fileName = file.name.toLowerCase();

  // 处理多重扩展名，例如 .mtx.gz
  const fileParts = fileName.split(".");
  const fileExtension =
    fileParts.length > 2
      ? `.${fileParts.slice(-2).join(".")}` // 如果有两个及以上部分，取最后两部分作为扩展名
      : `.${fileParts.pop()}`; // 否则取最后一部分

  const isValidType = restrictions.allowedExtensions.includes(fileExtension);
  const matchesRequiredName = restrictions.requiredFileNames.some(
    (name) => name === "*" || fileName.includes(name)
  );

  if (!isValidType || !matchesRequiredName) {
    console.error(
      i18n.global.t("fileUploadConsoleError", {
        type: fileName,
        allowedTypes: restrictions.allowedExtensions.join(", "),
        requiredNames: restrictions.requiredFileNames.join(", "),
      })
    );
    return false;
  }
  return true;
};

// 暴露给父组件
defineExpose({
  handleChange,
});
</script>

<style scoped>
.upload-demo {
  width: 60%;
  margin: 40px 0;
}
</style>
