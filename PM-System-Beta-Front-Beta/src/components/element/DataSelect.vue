<template>
  <div class="platform-left-tag">
    <el-tag type="info" style="width: fit-content">{{ tagText }}</el-tag>
    <el-select
      :model-value="modelValue"
      filterable
      :placeholder="placeholderText"
      style="width: 240px"
      @update:model-value="updateValue"
    >
      <el-option
        v-for="item in options"
        :key="item.value"
        :label="item.label"
        :value="item.value"
      />
    </el-select>
  </div>
</template>

<script setup>
const props = defineProps({
  tagText: {
    type: String,
    default: "Data set",
  },
  placeholderText: {
    type: String,
    default: "Select",
  },
  options: {
    type: Array,
    default: () => [],
  },
  modelValue: {
    type: [String, Number],
    default: "",
  },
});

const emit = defineEmits(["update:modelValue"]);

const updateValue = (value) => {
  emit("update:modelValue", value);
};
</script>

<style scoped>
.platform-left-tag {
  width: 240px;
  display: flex;
  flex-direction: column;
  box-sizing: border-box;
  margin: 30px 0;
}
.platform-left-tag :deep(.el-select__wrapper)::after {
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
</style>
