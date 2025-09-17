<template>
  <div class="info-panel-container">
    <div v-if="mode === 'gene' && geneStats" class="info-section">
      <h3 class="section-header">{{ $t("geneExpression") }}</h3>
      <div class="section-content">
        <div class="stat-item">
          <span>{{ $t('expressingCells') }}:</span>
          <strong>{{ geneStats.expressingCells }} / {{ geneStats.totalCells }}</strong>
        </div>
        <div class="stat-item">
          <span>{{ $t('proportion') }}:</span>
          <strong>{{ geneStats.percentage.toFixed(2) }}%</strong>
        </div>
      </div>
    </div>
    <div v-if="mode === 'cluster'" class="info-section">
      <h3 class="section-header">{{ $t("selectedCluster") }}</h3>
      <div class="section-content cluster-list">
        <div v-if="selectedClusters.length">
          <div
            v-for="cluster in selectedClusters"
            :key="cluster.name"
            class="cluster-item"
          >
            <div class="cluster-color-box" :style="{ backgroundColor: cluster.color }"></div>
            <span>{{ cluster.name }}</span>
          </div>
        </div>
        <span v-else class="no-selection-text">{{ $t("none") || "None" }}</span>
      </div>
    </div>
    <div class="info-section points-section">
      <h3 class="section-header">{{ $t("selectedPoints") }} ({{ selectedPoints.length }})</h3>
      <div class="section-content points-list" v-html="formattedSelectedPoints"></div>
    </div>
  </div>
</template>

<script>
export default {
  name: "InfoPanel",
  props: {
    mode: String,
    selectedClusters: Array,
    selectedPoints: Array,
    geneStats: Object,
  },
  computed: {
    formattedSelectedPoints() {
      if (!this.selectedPoints || this.selectedPoints.length === 0) {
        return `<div class="no-selection-text">${this.$t("none") || "None"}</div>`;
      }
      return this.selectedPoints.map(p => `<div class="point-item">${p.cellName}</div>`).join("");
    },
  },
};
</script>

<style scoped>
.info-panel-container {
  width: 300px;
  background-color: var(--el-info-panel-bg-color, #ffffff);
  border-left: 1px solid var(--el-info-panel-border-color, #e4e7ed);
  padding: 15px;
  height: calc(100vh - 96px);
  display: flex;
  flex-direction: column;
  color: var(--el-text-color-primary, #303133);
  box-sizing: border-box;
  gap: 15px;
}

.info-section {
  display: flex;
  flex-direction: column;
  border: 1px solid var(--el-info-panel-border-color, #e0e0e0);
  border-radius: 4px;
  background-color: var(--el-info-panel-section-bg-color, #f9f9f9);
}

.section-header {
  font-size: 16px;
  font-weight: 600;
  padding: 10px;
  margin: 0;
  border-bottom: 1px solid var(--el-info-panel-border-color, #e0e0e0);
  background-color: var(--el-info-panel-header-bg-color, #f5f5f5);
  border-top-left-radius: 4px;
  border-top-right-radius: 4px;
}

.section-content {
  padding: 10px;
}

.stat-item {
  display: flex;
  justify-content: space-between;
  padding: 4px 0;
}

.cluster-list {
  display: flex;
  flex-direction: column;
  gap: 8px;
}

.cluster-item {
  display: flex;
  align-items: center;
}

.cluster-color-box {
  width: 20px;
  height: 20px;
  margin-right: 10px;
  border: 1px solid #ccc;
}

.points-section {
  flex: 1;
  min-height: 0; 
}

.points-list {
  height: 100%;
  overflow-y: auto;
  padding-right: 5px; 
}

.point-item {
  padding: 2px 0;
  font-size: 12px;
}

.no-selection-text {
  color: #909399;
  font-style: italic;
}
</style>
