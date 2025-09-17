<template>
  <div class="legend-container">
    <div v-for="(item, index) in items" :key="item.id" class="legend-item" @click="toggleClusterVisibility(item.id)">
      <el-color-picker
        v-model="item.color"
        @change="onColorChange(item, $event)"
        @click.stop
        size="small"
        style="margin-right: 10px;"
      />
      <input
        v-if="editingItem === item.id"
        v-model="editingName"
        @keyup.enter="saveEdit(item)"
        @blur="saveEdit(item)"
        @click.stop
        class="legend-input"
        v-focus
      />
      <span v-else class="legend-label" :style="{ opacity: isVisible(item.id) ? 1 : 0.5 }" @dblclick.stop="startEdit(item)">{{ item.name }}</span>
    </div>
  </div>
</template>

<script>
export default {
  name: 'Legend',
  props: {
    items: {
      type: Array,
      required: true,
      default: () => []
    },
    visibleClusters: {
        type: Array,
        required: true,
        default: () => []
    }
  },
  directives: {
    focus: {
      mounted(el) {
        el.focus();
      }
    }
  },
  data() {
    return {
      editingItem: null,
      editingName: '',
    };
  },
  methods: {
    isVisible(id) {
      return this.visibleClusters.includes(id);
    },
    toggleClusterVisibility(id) {
      if (this.editingItem) return;
      const newVisible = [...this.visibleClusters];
      const index = newVisible.indexOf(id);
      if (index > -1) {
        newVisible.splice(index, 1);
      } else {
        newVisible.push(id);
      }
      this.$emit('update:visible-clusters', newVisible);
    },
    startEdit(item) {
      this.editingItem = item.id;
      this.editingName = item.name;
    },
    saveEdit(item) {
      if (this.editingItem === item.id) {
        const updatedItem = { ...item, name: this.editingName };
        this.$emit('update:item', { item: updatedItem });
        this.editingItem = null;
      }
    },
    onColorChange(item, color) {
      const updatedItem = { ...item, color: color };
      this.$emit('update:item', { item: updatedItem });
    }
  }
}
</script>

<style scoped>
.legend-container {
  position: absolute;
  top: 10px;
  right: 10px;
  background-color: var(--el-legend-bg-color, rgba(255, 255, 255, 0.9));
  padding: 10px;
  border-radius: 8px;
  box-shadow: 0 2px 10px rgba(0,0,0,0.1);
  max-height: 400px;
  overflow-y: auto;
  z-index: 10;
}
.legend-item {
  display: flex;
  align-items: center;
  margin-bottom: 5px;
  padding: 5px;
  border-radius: 4px;
  transition: background-color 0.2s ease-in-out;
}
.legend-item:hover {
  background-color: var(--el-legend-hover-bg-color, #f5f5f5);
}
.legend-label {
  font-size: 14px;
  cursor: pointer;
}
.legend-input {
  border: 1px solid var(--el-legend-input-border-color, #ddd);
  padding: 4px 6px;
  font-size: 14px;
  width: 120px;
  border-radius: 4px;
}
</style>
