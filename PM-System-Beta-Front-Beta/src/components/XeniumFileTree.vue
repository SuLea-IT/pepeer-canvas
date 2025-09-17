<template>
  <el-card class="box-card">
    <template #header>
      <div class="card-header">
        <span>Xenium Data</span>
        <el-button class="button" text @click="loadDirectoryData">Refresh</el-button>
      </div>
    </template>
    <el-tabs v-model="activeTab">
      <el-tab-pane label="Tree View" name="tree">
        <el-tree
          :data="treeData"
          :props="defaultProps"
          @node-click="handleNodeClick"
          v-loading="loading"
        />
      </el-tab-pane>
      <el-tab-pane label="Table View" name="table">
        <el-table :data="tableData" v-loading="loading" style="width: 100%">
          <el-table-column prop="name" label="Folder Name" />
          <el-table-column prop="fileCount" label="File Count" />
          <el-table-column label="Actions">
            <template #default="scope">
              <el-button size="small" @click="handleDownloadFolder(scope.row.path)">Download</el-button>
            </template>
          </el-table-column>
        </el-table>
      </el-tab-pane>
    </el-tabs>
  </el-card>
</template>

<script setup>
import { ref, onMounted } from 'vue';
import { ElMessage } from 'element-plus';
import { fetchDirectory, downloadFolder } from '../services/xeniumApi';

const activeTab = ref('tree');
const treeData = ref([]);
const tableData = ref([]);
const loading = ref(true);
const defaultProps = {
  children: 'children',
  label: 'label',
};

const emit = defineEmits(['file-select']);

const loadDirectoryData = async () => {
  try {
    loading.value = true;
    const data = await fetchDirectory();
    
    treeData.value = data.map(dir => ({
      label: dir.name,
      children: dir.files.map(file => ({
        label: file.name,
        isLeaf: true,
        file: { ...file, path: `${dir.path}/${file.name}` },
      })),
    }));

    tableData.value = data.map(dir => ({
      name: dir.name,
      path: dir.path,
      fileCount: dir.files.length,
    }));

  } catch (error) {
    ElMessage.error('Failed to load directory: ' + error.message);
  } finally {
    loading.value = false;
  }
};

const handleNodeClick = (data) => {
  if (data.isLeaf) {
    emit('file-select', data.file);
  }
};

const handleDownloadFolder = async (path) => {
  try {
    await downloadFolder(path);
    ElMessage.success('Folder download started.');
  } catch (error) {
    ElMessage.error(`Download failed: ${error.message}`);
  }
};

onMounted(() => {
  loadDirectoryData();
});
</script>

<style scoped>
.card-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
}
</style>