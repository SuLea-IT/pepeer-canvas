<template>
  <div class="platform-container">
    <ControlPanel
      :dataTypes="dataTypes"
      :clusterOptions="clusterOptions"
      :show-proportion="showProportionChart"
      @update:dataSource="dataSource = $event"
      @update:dataType="dataType = $event"
      @update:mode="mode = $event"
      @update:pointSize="pointSize = $event"
      @gene-input="geneName = $event"
      @gene-set-input="geneSet = $event"
      @update:visible-clusters="visibleClusters = $event"
      @update:show-cluster-bg="showClusterBg = $event"
      @update:gene-opacity="geneOpacity = $event"
      @update:cluster-bg-opacity="clusterBgOpacity = $event"
      @update:show-cluster-labels="showClusterLabels = $event"
      @update:show-proportion="showProportionChart = $event"
      @update:lod-enabled="lodEnabled = $event"
      @update:step-reveal-enabled="stepRevealEnabled = $event"
      @download-request="downloadWithLegend"
    />
    <div class="canvas-wrapper">
      <div v-if="showProportionChart" class="proportion-chart-container">
        <div class="chart-header">
          <el-button-group>
            <el-button size="small" @click="chartType = 'pie'" :type="chartType === 'pie' ? 'primary' : 'default'">Pie</el-button>
            <el-button size="small" @click="chartType = 'bar'" :type="chartType === 'bar' ? 'primary' : 'default'">Bar</el-button>
          </el-button-group>
          <el-button @click="closeChart" class="close-button" size="small" circle>
            <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" width="16" height="16"><path d="M19 6.41L17.59 5 12 10.59 6.41 5 5 6.41 10.59 12 5 17.59 6.41 19 12 13.41 17.59 19 19 17.59 13.41 12z"/></svg>
          </el-button>
        </div>
        <ClusterProportionChart
          :chartData="chartReadyData"
          :visibleClusters="visibleClusters"
          :chartType="chartType"
          @update:visible-clusters="visibleClusters = $event"
        />
      </div>
      <CanvasDisplay
        ref="canvasDisplay"
        :dataSource="dataSource"
        :dataType="dataType"
        :mode="mode"
        :pointSize="pointSize"
        :geneName="geneName"
        :geneSet="geneSet"
        :selectedClustersProp="selectedClusters"
        :selectedPointsProp="selectedPoints"
        :visibleClusters="visibleClusters"
        :showClusterBg="showClusterBg"
        :showClusterLabels="showClusterLabels"
        :geneOpacity="geneOpacity"
        :clusterBgOpacity="clusterBgOpacity"
        :lodEnabled="lodEnabled"
        :stepRevealEnabled="stepRevealEnabled"
        :stepRevealCount="stepRevealCount"
        @update:selected-clusters="onSelectionUpdate"
        @update:selected-points="selectedPoints = $event"
        @update:cluster-meta="onClusterMetaUpdate"
        @update:all-points="allPointsData = $event"
        @update:gene-stats="geneStats = $event"
      />
      <Legend
        v-if="clusterMeta.length > 0"
        :items="clusterMeta"
        :visibleClusters="visibleClusters"
        @update:visible-clusters="visibleClusters = $event"
        @update:item="handleItemUpdate"
      />
    </div>
    <InfoPanel
      :mode="mode"
      :selectedClusters="selectedClusters"
      :selectedPoints="selectedPoints"
      :geneStats="geneStats"
    />
  </div>
</template>

<script>
import ControlPanel from "../components/ControlPanel.vue";
import CanvasDisplay from "../components/CanvasDisplay.vue";
import InfoPanel from "../components/InfoPanel.vue";
import Legend from "../components/Legend.vue";
import ClusterProportionChart from "../components/ClusterProportionChart.vue";
import axios from "axios";
import { apiConfig } from "../config/apiConfig.js";

export default {
  name: "Platform",
  components: {
    ControlPanel,
    CanvasDisplay,
    InfoPanel,
    Legend,
    ClusterProportionChart,
  },
  data() {
    return {
      lodEnabled: true,
      stepRevealEnabled: true,
      stepRevealCount: 0,
      stepRevealTimer: null,
      stepRevealInterval: 400,
      dataSource: "data",
      dataType: "",
      mode: "cluster",
      pointSize: 1,
      geneName: "",
      geneSet: [],
      selectedClusters: [],
      selectedPoints: [],
      dataTypes: [],
      clusterOptions: [],
      visibleClusters: [],
      clusterMeta: [],
      showClusterBg: false,
      showClusterLabels: false,
      geneOpacity: 1,
      clusterBgOpacity: 1,
      showProportionChart: false,
      allPointsData: [],
      chartType: 'pie',
      geneStats: null,
    };
  },
  mounted() {
    this.fetchDataTypes();
  },
  beforeUnmount() {
    this.stopStepRevealAutoPlay();
  },
  methods: {
    async fetchDataTypes() {
      try {
        const response = await axios.get(`${apiConfig.baseURL}/data/types`);
        if (response.data.success && response.data.data.length > 0) {
          this.dataTypes = response.data.data;
          this.dataType = this.dataTypes[0];
        }
      } catch (error) {
        console.error("获取数据类型失败:", error);
      }
    },
    closeChart() {
      this.showProportionChart = false;
      this.chartType = 'pie'; // Reset to default
    },
    onClusterMetaUpdate(meta) {
      this.clusterMeta = meta;
      this.clusterOptions = meta.map((g) => ({ value: g.id, label: g.name }));
      this.visibleClusters = meta.map(g => g.id); // Default all to visible
      // Reset step reveal on new data to show grey shapes first
      this.stepRevealEnabled = true;
      this.stepRevealCount = 0;
      this.startStepRevealAutoPlay();
    },
    startStepRevealAutoPlay() {
      if (!this.stepRevealEnabled || this.mode !== 'cluster') return;
      const max = this.clusterMeta ? this.clusterMeta.length : 0;
      if (!max) return;
      this.stopStepRevealAutoPlay();
      this.stepRevealTimer = setInterval(() => {
        if (this.stepRevealCount >= max) {
          this.stopStepRevealAutoPlay();
          return;
        }
        this.stepRevealCount += 1;
      }, this.stepRevealInterval);
    },
    stopStepRevealAutoPlay() {
      if (this.stepRevealTimer) {
        clearInterval(this.stepRevealTimer);
        this.stepRevealTimer = null;
      }
    },
    onSelectionUpdate(selectedClusterIds) {
        // Filter full cluster objects based on ids from the event
        this.selectedClusters = this.clusterMeta.filter(metaItem => 
            selectedClusterIds.includes(metaItem.id)
        );
    },
    downloadWithLegend() {
      const umapCanvas = this.$refs.canvasDisplay.getCanvas();
      if (!umapCanvas) return;

      const legendItems = this.clusterMeta.filter((item) =>
        this.visibleClusters.includes(item.name)
      );
      const legendWidth = 200;
      const legendPadding = 20;
      const itemHeight = 30;
      const totalHeight = Math.max(
        umapCanvas.height,
        legendItems.length * itemHeight + legendPadding * 2
      );

      const offscreenCanvas = document.createElement("canvas");
      offscreenCanvas.width = umapCanvas.width + legendWidth;
      offscreenCanvas.height = totalHeight;
      const ctx = offscreenCanvas.getContext("2d");

      ctx.fillStyle = "white";
      ctx.fillRect(0, 0, offscreenCanvas.width, offscreenCanvas.height);
      ctx.drawImage(umapCanvas, 0, 0);

      let y = legendPadding;
      ctx.font = "16px Arial";
      ctx.fillStyle = "black";

      for (const item of legendItems) {
        ctx.fillStyle = item.color;
        ctx.fillRect(umapCanvas.width + legendPadding, y, 20, 20);
        ctx.fillStyle = "black";
        ctx.fillText(item.name, umapCanvas.width + legendPadding + 30, y + 15);
        y += itemHeight;
      }

      const link = document.createElement("a");
      link.download = `umap-with-legend-${new Date().toISOString()}.png`;
      link.href = offscreenCanvas.toDataURL("image/png");
      link.click();
    },
    handleItemUpdate({ index, item }) {
      const originalItemIndex = this.clusterMeta.findIndex(m => m.id === item.id);
      if (originalItemIndex !== -1) {
        this.clusterMeta.splice(originalItemIndex, 1, item);
      }
    },
  },
  watch: {
    stepRevealEnabled(newVal) {
      if (newVal) this.startStepRevealAutoPlay(); else this.stopStepRevealAutoPlay();
    },
    mode(newVal) {
      if (newVal === 'cluster') this.startStepRevealAutoPlay(); else this.stopStepRevealAutoPlay();
    }
  },
  computed: {
    chartReadyData() {
      if (!this.allPointsData.length || !this.clusterMeta.length) {
        return { data: [], total: 0 };
      }
      
      const metaMapById = new Map(this.clusterMeta.map(m => [m.id, m]));
      const counts = {};

      for (const point of this.allPointsData) {
        const clusterId = point.clusterId;
        if (counts[clusterId] === undefined) {
          counts[clusterId] = 0;
        }
        counts[clusterId]++;
      }

      const chartData = Object.entries(counts).map(([id, value]) => {
        const meta = metaMapById.get(Number(id));
        return {
          id: Number(id),
          name: meta ? meta.name : `Cluster ${id}`,
          value: value,
          itemStyle: {
            color: meta ? meta.color : '#ccc'
          }
        };
      });

      return {
        data: chartData,
        total: this.allPointsData.length,
      };
    }
  }
};
</script>

<style scoped>
.platform-container {
  display: flex;
  height: calc(100vh - 96px);
}
.canvas-wrapper {
  flex: 1;
  position: relative;
  display: flex;
}
.canvas-wrapper > :deep(.canvas-div) {
  flex: 1;
  min-height: 0;
}
.proportion-chart-container {
  position: absolute;
  top: 20px;
  right: 20px;
  width: 480px;
  height: 450px;
  background-color: rgba(255, 255, 255, 0.95);
  border: 1px solid #ccc;
  border-radius: 8px;
  box-shadow: 0 4px 12px rgba(0,0,0,0.15);
  z-index: 100;
  padding: 15px;
  display: flex;
  flex-direction: column;
}
.chart-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  margin-bottom: 10px;
}
.close-button {
  border: none;
  background: transparent;
  color: #000; /* Default color for light mode */
}

.dark .proportion-chart-container {
  background-color: rgba(40, 40, 40, 0.95);
  border-color: #555;
  box-shadow: 0 4px 12px rgba(0,0,0,0.5);
}

.dark .close-button {
  color: #fff; /* White color for dark mode */
}
</style>
