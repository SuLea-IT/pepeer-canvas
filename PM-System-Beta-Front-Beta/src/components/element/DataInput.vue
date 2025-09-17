<template>
  <div class="platform-left-tag">
    <el-tag type="info" style="width: fit-content">{{ tagText }}</el-tag>
    <el-input
      :model-value="modelValue"
      :placeholder="placeholderText"
      style="width: 240px"
      @input="handleInput"
      clearable
    />
    <div v-if="errorMessage" class="error-message">{{ errorMessage }}</div>
  </div>
</template>

<script setup>
import { ref, watch } from "vue";

const props = defineProps({
  tagText: {
    type: String,
    default: "Data set",
  },

  modelValue: {
    type: [String, Number],
    default: "",
  },
  placeholderText: {
    type: String,
    default: "Select",
  },
  validateFn: {
    type: Function,
    default: null,
  },
});

const emit = defineEmits(["update:modelValue"]);

const errorMessage = ref("");

const handleInput = (value) => {
  if (props.validateFn) {
    const isValid = props.validateFn(value);
    errorMessage.value = isValid ? "" : "无效邮箱";
  }
  emit("update:modelValue", value);
};

watch(
  () => props.modelValue,
  (newVal) => {
    if (props.validateFn) {
      const isValid = props.validateFn(newVal);
      errorMessage.value = isValid ? "" : "无效邮箱";
    }
  }
);
</script>

<style scoped>
.platform-left-tag {
  width: 240px;
  display: flex;
  flex-direction: column;
  box-sizing: border-box;
  margin: 30px 0;
}

.platform-left-tag :deep(.el-input__wrapper)::after {
  content: "";
  position: absolute;
  bottom: -6px;
  left: 10%;
  height: 1px;
  background-color: #d5d5d5;
  transform: scaleX(1);
  transition: all 0.3s ease;
  width: 80%;
}

.error-message {
  color: red;
  font-size: 12px;
  margin-top: 5px;
}
</style>
