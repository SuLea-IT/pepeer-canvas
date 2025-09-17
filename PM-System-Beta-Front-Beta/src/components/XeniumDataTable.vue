<template>
  <div>
    <!-- Header -->


    <!-- Main Content -->
    <el-card class="box-card">
      <template #header>
        <div class="card-header">
          <!-- Summary Info -->
          <div class="summary-info">
            <el-icon :size="24" style="margin-right: 8px;"><Folder /></el-icon>
            <span>Sample Folders</span>
            <el-tag type="success" style="margin-left: 10px;">{{ summary.folderCount }} Folders</el-tag>
            <el-tag type="warning" style="margin-left: 10px;">{{ summary.fileCount }} Files</el-tag>
            <el-tag type="info" style="margin-left: 10px;">{{ formatSize(summary.totalSize) }}</el-tag>
          </div>
          <!-- Action Buttons -->
          <div class="action-buttons">
            <el-button @click="loadDirectoryData">Refresh</el-button>
            <el-button type="primary" @click="handleDownloadAll">Download All</el-button>
          </div>
        </div>
      </template>

      <!-- Table -->
      <el-table
        :data="tableData"
        style="width: 100%"
        row-key="id"
        lazy
        :load="loadChildren"
        :tree-props="{ children: 'children', hasChildren: 'hasChildren' }"
        v-loading="loading"
      >
        <el-table-column prop="name" label="Folder Name / File Name" width="400">
          <template #default="scope">
            <div style="display: flex; align-items: center;">
              <el-icon v-if="scope.row.type === 'folder'" class="table-icon"><Folder /></el-icon>
              <el-icon v-else class="table-icon"><Document /></el-icon>
              <span>{{ scope.row.name }}</span>
              <el-tag v-if="scope.row.type === 'file'" size="small" style="margin-left: 10px;">
                {{ scope.row.name.split('.').pop().toUpperCase() }}
              </el-tag>
               <el-tag v-if="scope.row.type === 'folder'" type="info" size="small" style="margin-left: 10px;">
                {{ scope.row.children.length }} Files
              </el-tag>
            </div>
          </template>
        </el-table-column>
        <el-table-column prop="size" label="Total Size / Size" align="right">
           <template #default="scope">
            <span>{{ formatSize(scope.row.type === 'folder' ? scope.row.totalSize : scope.row.size) }}</span>
          </template>
        </el-table-column>
        <el-table-column label="Actions" align="right">
          <template #default="scope">
            <el-button v-if="scope.row.type === 'folder'" type="primary" size="small" @click="handleDownloadFolder(scope.row.path)">
              Download ZIP
            </el-button>
            <el-button v-else size="small" @click="handlePreview(scope.row)">
              Preview
            </el-button>
          </template>
        </el-table-column>
      </el-table>
    </el-card>
  </div>
</template>

<script setup>
import { ref, onMounted, computed } from 'vue';
import { ElMessage } from 'element-plus';
import { fetchDirectory, downloadFolder } from '../services/xeniumApi';
import { Folder, Document } from '@element-plus/icons-vue';

const loading = ref(false);
const rawData = ref([]);
const tableData = ref([]);

const loadDirectoryData = async () => {
  try {
    loading.value = true;
    const data = await fetchDirectory();
    rawData.value = data;
    tableData.value = processData(data);
  } catch (error) {
    ElMessage.error('Failed to load directory: ' + error.message);
  } finally {
    loading.value = false;
  }
};

const processData = (data) => {
  return data.map(item => {
    const folder = {
      ...item,
      id: item.path,
      type: 'folder',
      hasChildren: item.files && item.files.length > 0,
      children: item.files.map(file => ({
        ...file,
        id: file.path,
        type: 'file',
      })),
    };
    folder.totalSize = folder.children.reduce((acc, file) => acc + file.size, 0);
    return folder;
  });
};

const loadChildren = (row, treeNode, resolve) => {
  resolve(row.children);
};

const handleDownloadFolder = async (path) => {
  try {
    await downloadFolder(path);
    ElMessage.success('Folder download started.');
  } catch (error) {
    ElMessage.error(`Download failed: ${error.message}`);
  }
};

const handlePreview = (file) => {
  // Placeholder for preview functionality
  ElMessage.info(`Previewing ${file.name}`);
};

const summary = computed(() => {
  const folderCount = rawData.value.length;
  const fileCount = rawData.value.reduce((acc, dir) => acc + dir.files.length, 0);
  const totalSize = rawData.value.reduce((acc, dir) => {
    return acc + dir.files.reduce((fileAcc, file) => fileAcc + file.size, 0);
  }, 0);
  return { folderCount, fileCount, totalSize };
});

const formatSize = (bytes) => {
  if (bytes === 0) return '0 B';
  const k = 1024;
  const sizes = ['B', 'KB', 'MB', 'GB', 'TB'];
  const i = Math.floor(Math.log(bytes) / Math.log(k));
  return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
};

const handleDownloadAll = () => {
  // Placeholder for download all functionality
  ElMessage.info('Download All clicked');
};

onMounted(() => {
  loadDirectoryData();
});
</script>

<style scoped>
.header {
  margin-bottom: 20px;
}
.title {
  font-size: 24px;
  font-weight: bold;
}
.card-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
}
.summary-info span {
  font-size: 18px;
  font-weight: bold;
  margin-right: 10px;
}
.table-icon {
  margin-right: 8px;
}

.dark :deep(.el-card),
.dark :deep(.el-table),
.dark :deep(.el-table th),
.dark :deep(.el-table td) {
  --el-card-bg-color: #141414;
  --el-bg-color: #141414;
  --el-text-color-primary: #ffffff;
  --el-text-color-regular: #cfcfcf;
  --el-border-color-lighter: #424242;
  color: var(--el-text-color-primary);
}

.dark :deep(.el-table__row) {
  background-color: #1f1f1f;
}

.dark :deep(.el-table__row--level-1) {
  background-color: #2a2a2a;
}

.dark :deep(.el-table__row:hover > td) {
  background-color: #383838 !important;
}

:deep(.el-table__row--level-1) .cell {
  padding-left: 40px !important;
}

:deep(.el-tag) {
  border-radius: 12px;
  border: none;
}

:deep(.el-button) {
  border-radius: 6px;
}
</style>