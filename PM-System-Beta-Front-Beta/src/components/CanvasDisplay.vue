<template>
  <div class="canvas-div" ref="canvasDiv">
    <div class="toolbar">
      <button @click="switchTool('move')" class="tool-button" :class="{ active: !isSelecting }" title="Move">
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" width="24" height="24"><path d="M12 2L8 6h3v3H6v4h3v3H8l4 4 4-4h-3v-3h4V9h-3V6h3L12 2z"/></svg>
      </button>
      <button @click="switchTool('select')" class="tool-button" :class="{ active: isSelecting }" title="Select">
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" width="24" height="24"><path d="M16 4H8C5.79 4 4 5.79 4 8v8c0 2.21 1.79 4 4 4h8c2.21 0 4-1.79 4-4V8c0-2.21-1.79-4-4-4zm-2 10h-4v-4h4v4z"/></svg>
      </button>
      <button @click="resetView" class="tool-button" title="Reset">
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" width="24" height="24"><path d="M12 5V1L7 6l5 5V7c3.31 0 6 2.69 6 6s-2.69 6-6 6-6-2.69-6-6H4c0 4.42 3.58 8 8 8s8-3.58 8-8-3.58-8-8-8z"/></svg>
      </button>
    </div>
    <div ref="deckContainer" class="deck-container"></div>
    <div ref="selectionBox" class="selection" id="selection-box" style="display: none"></div>
    <div v-if="mode !== 'cluster' && allPoints.length > 0" class="legend-bar">
      <div class="legend-gradient" :style="legendGradientStyle"></div>
      <div class="legend-labels">
        <span>{{ geneValueMin.toFixed(2) }}</span>
        <span>{{ ((geneValueMin + geneValueMax) / 2).toFixed(2) }}</span>
        <span>{{ geneValueMax.toFixed(2) }}</span>
      </div>
    </div>
  </div>
</template>

<script>
import { ElLoading } from "element-plus";
import { apiConfig } from "../config/apiConfig";
import { quadtree as d3_quadtree } from "d3-quadtree";
import { Deck, OrthographicView } from '@deck.gl/core';
import { ScatterplotLayer, TextLayer } from '@deck.gl/layers';

export default {
  name: "CanvasDisplay",
  props: {
    dataSource: String,
    dataType: String,
    mode: String,
    pointSize: Number,
    geneName: String,
    geneSet: Array,
    selectedClustersProp: Array,
    selectedPointsProp: Array,
    visibleClusters: Array,
    showClusterBg: Boolean,
    showClusterLabels: Boolean,
    geneOpacity: {
      type: Number,
      default: 1.0
    },
    clusterBgOpacity: {
      type: Number,
      default: 1.0
    },
    minColor: {
      type: String,
      default: '#0000ff'
    },
    midColor: {
      type: String,
      default: '#ffff00'
    },
    maxColor: {
      type: String,
      default: '#ff0000'
    },
    // Cluster step-reveal settings
    stepRevealEnabled: {
      type: Boolean,
      default: true,
    },
    stepRevealCount: {
      // number of clusters to reveal in color (0 => only grey background)
      type: Number,
      default: 0,
    },
    // LOD and aggregation configs (cluster mode)
    lodEnabled: {
      type: Boolean,
      default: true,
    },
    lodScaleThreshold: {
      // when scaleFactor < threshold, enable LOD
      type: Number,
      default: 25,
    },
    gridSizePx: {
      type: Number,
      default: 6,
    },
    clusterAggFactor: {
      // e.g. 3 => 3 points render as 1 representative
      type: Number,
      default: 3,
    },
    maxAggPerCell: {
      // safety cap to avoid too many reps per grid cell
      type: Number,
      default: 6,
    },
  },
  data() {
    return {
      allPoints: [],
      pointsByCluster: new Map(),
      clusterPoints: [],
      clusterLabels: [],
      clusterOrder: [],
      // deck.gl
      deck: null,
      pointLayerId: 'points-layer',
      bgLayerId: 'bg-layer',
      labelLayerId: 'label-layer',
      viewState: null,
      isDragging: false,
      isSelecting: false,
      startX: null,
      startY: null,
      quadtree: null,
      geneValueMin: 0,
      geneValueMax: 0,
      // rAF-based redraw control
      _rafId: 0,
      _pendingRedraw: false,
    };
  },
  computed: {
    legendGradientStyle() {
      return {
        background: `linear-gradient(to right, ${this.minColor}, ${this.midColor}, ${this.maxColor})`
      };
    }
  },
  methods: {
    // small deterministic hash for stable jitter
    hashMix(key, i = 0) {
      let h = 2166136261;
      const s = key + '#' + i;
      for (let k = 0; k < s.length; k++) {
        h ^= s.charCodeAt(k);
        h = Math.imul(h, 16777619);
      }
      // final avalanching
      h ^= h >>> 15; h = Math.imul(h, 2246822519);
      h ^= h >>> 13; h = Math.imul(h, 3266489917);
      h ^= h >>> 16;
      return (h >>> 0);
    },
    async loadAndDrawUMAPData() {
      if (!this.dataSource || !this.dataType || !this.mode) return;

      const loading = ElLoading.service({ lock: true, text: "Loading..." });

      if (this.mode === 'gene' && this.showClusterBg) {
        try {
          const clusterUrl = apiConfig.endpoints.getClusterData(this.dataType);
          const response = await fetch(clusterUrl);
          if (!response.ok) throw new Error(`HTTP error! status: ${response.status}`);
          const arrayBuffer = await response.arrayBuffer();
          this.clusterPoints = this.parseClusterDataForBackground(arrayBuffer);
        } catch (error) {
          console.error("Failed to load cluster background data:", error);
          this.$message.error(`Failed to load cluster background. ${error.message}`);
          this.clusterPoints = []; // Ensure it's reset on failure
        }
      } else {
        this.clusterPoints = [];
      }
      
      let url;
      let options = {};

      // Determine the correct API endpoint based on the mode
      if (this.mode === 'cluster') {
          url = apiConfig.endpoints.getClusterData(this.dataType);
      } else if (this.mode === 'gene' && this.geneName) {
          url = apiConfig.endpoints.getGeneData(this.dataSource, this.geneName, this.dataType);
      } else if (this.mode === 'gene' && this.geneSet && this.geneSet.length > 0) {
          const geneSetApiConfig = apiConfig.endpoints.getGeneSetData(this.dataSource, this.geneSet, this.dataType);
          url = geneSetApiConfig.url;
          options = { method: geneSetApiConfig.method, headers: geneSetApiConfig.headers, body: geneSetApiConfig.body };
      } else {
          loading.close();
          return;
      }

      try {
        const response = await fetch(url, options);
        if (!response.ok) throw new Error(`HTTP error! status: ${response.status}`);
        const arrayBuffer = await response.arrayBuffer();
        
        this.parseAndProcessSimpleBinary(arrayBuffer);
        
        this.$nextTick(() => {
          this.autoAdjustView();
          this.initCanvasEvents();
        });

      } catch (error) {
        console.error("Failed to load data:", error);
        this.$message.error(`Failed to load or parse data. ${error.message}`);
      } finally {
        loading.close();
      }
    },

    parseAndProcessSimpleBinary(arrayBuffer) {
        const dataView = new DataView(arrayBuffer);
        const numPoints = arrayBuffer.byteLength / 12; // Each point is 3 floats * 4 bytes/float = 12 bytes
        
        if (arrayBuffer.byteLength % 12 !== 0) {
            throw new Error("Binary file size is not a multiple of 12, format is incorrect.");
        }

        const points = [];
        let minV = Infinity, maxV = -Infinity;

        for (let i = 0; i < numPoints; i++) {
            const offset = i * 12;
            const x = dataView.getFloat32(offset, true);
            const y = dataView.getFloat32(offset + 4, true);
            const value = dataView.getFloat32(offset + 8, true);
            
            points.push({ x, y, value });
            if (value < minV) minV = value;
            if (value > maxV) maxV = value;
        }

        this.allPoints = points.map((p, i) => ({
            x: p.x,
            y: p.y,
            value: p.value,
            cellName: `cell_${i}`, // Generic cell name
            color: 'gray' // Default color
        }));

        if (this.mode === 'cluster') {
            const uniqueClusters = [...new Set(this.allPoints.map(p => p.value))].sort((a, b) => a - b);
            const colors = this.generateDistinctColors(uniqueClusters.length);
            const clusterMap = new Map(uniqueClusters.map((id, index) => [id, colors[index]]));
            this.clusterOrder = uniqueClusters;
            
            const clusterMeta = uniqueClusters.map((id, index) => ({
                id: id,
                name: `Cluster ${id}`,
                color: colors[index]
            }));

            // Build per-cluster buckets once for faster drawing
            this.pointsByCluster = new Map();
            this.allPoints.forEach(p => {
                p.color = clusterMap.get(p.value);
                p.cellName = `${p.value}_${p.cellName}`;
                p.clusterId = p.value; // Keep original id
                if (!this.pointsByCluster.has(p.clusterId)) this.pointsByCluster.set(p.clusterId, []);
                this.pointsByCluster.get(p.clusterId).push(p);
            });
            
            this.$emit('update:cluster-meta', clusterMeta);

        } else { // Gene or Gene-set mode
            this.geneValueMin = minV;
            this.geneValueMax = maxV;
            // Colors for gene mode are calculated on-the-fly in redrawUMAP
            this.$emit('update:cluster-meta', []); // No legend for gene mode
        }

        this.quadtree = d3_quadtree()
            .x(d => d.x)
            .y(d => d.y)
            .addAll(this.allPoints);
        
        this.$emit('update:all-points', this.allPoints);

        if (this.mode === 'gene') {
            const expressingCells = this.allPoints.filter(p => p.value > 0).length;
            const totalCells = this.allPoints.length;
            const percentage = totalCells > 0 ? (expressingCells / totalCells) * 100 : 0;
            this.$emit('update:gene-stats', {
                expressingCells,
                totalCells,
                percentage,
            });
        } else {
            this.$emit('update:gene-stats', null);
        }
    },

    parseClusterDataForBackground(arrayBuffer) {
        const dataView = new DataView(arrayBuffer);
        const numPoints = arrayBuffer.byteLength / 12;
        const points = [];
        const clusterSums = {};

        for (let i = 0; i < numPoints; i++) {
            const offset = i * 12;
            const x = dataView.getFloat32(offset, true);
            const y = dataView.getFloat32(offset + 4, true);
            const value = dataView.getFloat32(offset + 8, true);
            points.push({ x, y, value });

            if (!clusterSums[value]) {
                clusterSums[value] = { x: 0, y: 0, count: 0 };
            }
            clusterSums[value].x += x;
            clusterSums[value].y += y;
            clusterSums[value].count++;
        }

        this.clusterLabels = Object.entries(clusterSums).map(([id, data]) => ({
            id: `C${id}`,
            x: data.x / data.count,
            y: data.y / data.count,
        }));

        const uniqueClusters = [...new Set(points.map(p => p.value))].sort((a, b) => a - b);
        const colors = this.generateDistinctColors(uniqueClusters.length); // Do not apply opacity here
        const clusterMap = new Map(uniqueClusters.map((id, index) => [id, colors[index]]));

        points.forEach(p => {
            p.color = clusterMap.get(p.value);
        });
        return points;
    },

    generateDistinctColors(count, opacity = 1) {
        const colors = [];
        const hueStep = 360 / count;
        for (let i = 0; i < count; i++) {
            colors.push(`hsla(${i * hueStep}, 90%, 55%, ${opacity})`);
        }
        return colors;
    },

    getGeneColor(intensity, withAlpha = false) {
        const min = this.hexToRgb(this.minColor);
        const mid = this.hexToRgb(this.midColor);
        const max = this.hexToRgb(this.maxColor);

        let r, g, b;
        if (intensity < 0.5) {
            const t = intensity * 2;
            r = Math.round(min.r + (mid.r - min.r) * t);
            g = Math.round(min.g + (mid.g - min.g) * t);
            b = Math.round(min.b + (mid.b - min.b) * t);
        } else {
            const t = (intensity - 0.5) * 2;
            r = Math.round(mid.r + (max.r - mid.r) * t);
            g = Math.round(mid.g + (max.g - mid.g) * t);
            b = Math.round(mid.b + (max.b - mid.b) * t);
        }
        
        let alpha;
        if (this.showClusterBg) {
            alpha = intensity > 0.01 ? this.geneOpacity : 0.0;
        } else {
            alpha = this.geneOpacity;
        }

        return `rgba(${r}, ${g}, ${b}, ${alpha})`;
    },

    // For deck.gl: return [r,g,b,a]
    getGeneColorArray(intensity, opacity = 1.0) {
        const min = this.hexToRgb(this.minColor);
        const mid = this.hexToRgb(this.midColor);
        const max = this.hexToRgb(this.maxColor);
        let r, g, b;
        if (intensity < 0.5) {
            const t = intensity * 2;
            r = Math.round(min.r + (mid.r - min.r) * t);
            g = Math.round(min.g + (mid.g - min.g) * t);
            b = Math.round(min.b + (mid.b - min.b) * t);
        } else {
            const t = (intensity - 0.5) * 2;
            r = Math.round(mid.r + (max.r - mid.r) * t);
            g = Math.round(mid.g + (max.g - mid.g) * t);
            b = Math.round(mid.b + (max.b - mid.b) * t);
        }
        let alpha;
        if (this.showClusterBg) {
            alpha = intensity > 0.01 ? opacity : 0.0;
        } else {
            alpha = opacity;
        }
        return [r, g, b, Math.round(alpha * 255)];
    },

    cssColorToRgbaArray(colorStr, alphaOverride = null) {
        if (!colorStr) return [128,128,128,255];
        if (colorStr.startsWith('#')) {
            const rgb = this.hexToRgb(colorStr);
            const a = alphaOverride == null ? 255 : alphaOverride;
            return [rgb.r, rgb.g, rgb.b, a];
        }
        const hslaMatch = /hsla?\(([^)]+)\)/i.exec(colorStr);
        if (hslaMatch) {
            const parts = hslaMatch[1].split(',').map(s => s.trim());
            const h = parseFloat(parts[0]);
            const s = parseFloat(parts[1]) / 100;
            const l = parseFloat(parts[2]) / 100;
            const a = parts[3] !== undefined ? parseFloat(parts[3]) : 1;
            const [r, g, b] = this.hslToRgb(h / 360, s, l);
            const A = alphaOverride == null ? Math.round(a * 255) : alphaOverride;
            return [r, g, b, A];
        }
        const rgbaMatch = /rgba?\(([^)]+)\)/i.exec(colorStr);
        if (rgbaMatch) {
            const parts = rgbaMatch[1].split(',').map(s => s.trim());
            const r = parseInt(parts[0]);
            const g = parseInt(parts[1]);
            const b = parseInt(parts[2]);
            const a = parts[3] !== undefined ? parseFloat(parts[3]) : 1;
            const A = alphaOverride == null ? Math.round(a * 255) : alphaOverride;
            return [r, g, b, A];
        }
        return [128,128,128, alphaOverride == null ? 255 : alphaOverride];
    },

    hslToRgb(h, s, l) {
        let r, g, b;
        if (s === 0) {
            r = g = b = l; // achromatic
        } else {
            const hue2rgb = (p, q, t) => {
                if (t < 0) t += 1;
                if (t > 1) t -= 1;
                if (t < 1/6) return p + (q - p) * 6 * t;
                if (t < 1/2) return q;
                if (t < 2/3) return p + (q - p) * (2/3 - t) * 6;
                return p;
            };
            const q = l < 0.5 ? l * (1 + s) : l + s - l * s;
            const p = 2 * l - q;
            r = hue2rgb(p, q, h + 1/3);
            g = hue2rgb(p, q, h);
            b = hue2rgb(p, q, h - 1/3);
        }
        return [Math.round(r * 255), Math.round(g * 255), Math.round(b * 255)];
    },

    hexToRgb(hex) {
        const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
        return result ? {
            r: parseInt(result[1], 16),
            g: parseInt(result[2], 16),
            b: parseInt(result[3], 16)
        } : null;
    },

    autoAdjustView() {
      if (!this.allPoints || this.allPoints.length === 0 || !this.deck) return;
      const container = this.$refs.canvasDiv;
      if (!container) return;
      let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
      for (const p of this.allPoints) {
        if (p.x < minX) minX = p.x;
        if (p.x > maxX) maxX = p.x;
        if (p.y < minY) minY = p.y;
        if (p.y > maxY) maxY = p.y;
      }
      const dataWidth = Math.max(1e-6, maxX - minX);
      const dataHeight = Math.max(1e-6, maxY - minY);
      const padding = 50;
      const width = container.clientWidth || 1;
      const height = container.clientHeight || 1;
      const scaleX = (width - padding * 2) / dataWidth;
      const scaleY = (height - padding * 2) / dataHeight;
      const scale = Math.max(1e-6, Math.min(scaleX, scaleY));
      const zoom = Math.log2(scale);
      const centerX = (minX + maxX) / 2;
      const centerY = (minY + maxY) / 2;
      this.viewState = { target: [centerX, centerY, 0], zoom, minZoom: -10, maxZoom: 20 };
      // Only set initialViewState; let deck manage viewState for interactivity
      this.deck.setProps({ initialViewState: this.viewState });
      this.updateDeckLayers();
    },
    
    requestRedraw() {
      this.updateDeckLayers();
    },

    updateDeckLayers() {
      if (!this.deck) return;
      const selectedSet = new Set((this.selectedPointsProp || []).map(p => p.cellName));
      const hasSelection = selectedSet.size > 0;
      const visibleSet = new Set(this.visibleClusters || []);
      const mode = this.mode;
      const minV = this.geneValueMin;
      const maxV = this.geneValueMax;
      const self = this;

      // default coloring layer (used when not in step-reveal mode or when selection active)
      const pointsLayer = new ScatterplotLayer({
        id: this.pointLayerId,
        data: this.allPoints,
        pickable: true,
        radiusUnits: 'pixels',
        getRadius: this.pointSize,
        getPosition: d => [d.x, d.y, 0],
        getFillColor: d => {
          if (mode === 'cluster') {
            if (visibleSet.size && !visibleSet.has(d.clusterId)) return [0,0,0,0];
            if (hasSelection && !selectedSet.has(d.cellName)) return [200,200,200,128];
            return self.cssColorToRgbaArray(d.color, 255);
          } else {
            const range = maxV - minV;
            const intensity = range > 0 ? (d.value - minV) / range : 0;
            return self.getGeneColorArray(intensity, self.geneOpacity);
          }
        },
        updateTriggers: {
          getFillColor: [mode, hasSelection, [...selectedSet].join(','), [...visibleSet].join(','), minV, maxV, this.geneOpacity, this.minColor, this.midColor, this.maxColor],
          getRadius: [this.pointSize]
        }
      });

      const layers = [];

      // Step reveal mode for cluster visualization (no active selection)
      if (mode === 'cluster' && this.stepRevealEnabled && !hasSelection) {
        const revealCount = Math.max(0, Math.min(this.stepRevealCount, this.clusterOrder.length));
        if (revealCount === 0) {
          // 仅灰色形状，不上色任何簇
          layers.push(new ScatterplotLayer({
            id: this.pointLayerId + '-grey',
            data: this.allPoints,
            pickable: true,
            radiusUnits: 'pixels',
            getRadius: this.pointSize,
            getPosition: d => [d.x, d.y, 0],
            getFillColor: d => {
              if (visibleSet.size && !visibleSet.has(d.clusterId)) return [0,0,0,0];
              return [200,200,200,180];
            },
            updateTriggers: {
              getFillColor: [[...visibleSet].join(',')],
              getRadius: [this.pointSize]
            }
          }));
        } else {
          const revealSet = new Set(this.clusterOrder.slice(0, revealCount));

        // Base grey shape for visible clusters
        layers.push(new ScatterplotLayer({
          id: this.pointLayerId + '-grey',
          data: this.allPoints,
          pickable: true,
          radiusUnits: 'pixels',
          getRadius: this.pointSize,
          getPosition: d => [d.x, d.y, 0],
          getFillColor: d => {
            if (visibleSet.size && !visibleSet.has(d.clusterId)) return [0,0,0,0];
            return [200,200,200,180];
          },
          updateTriggers: {
            getFillColor: [[...visibleSet].join(',')],
            getRadius: [this.pointSize]
          }
        }));

        // Colored overlay for the first N clusters
        layers.push(new ScatterplotLayer({
          id: this.pointLayerId + '-reveal',
          data: this.allPoints,
          pickable: true,
          radiusUnits: 'pixels',
          getRadius: this.pointSize,
          getPosition: d => [d.x, d.y, 0],
          getFillColor: d => {
            if (visibleSet.size && !visibleSet.has(d.clusterId)) return [0,0,0,0];
            if (!revealSet.has(d.clusterId)) return [0,0,0,0];
            return self.cssColorToRgbaArray(d.color, 255);
          },
          updateTriggers: {
            getFillColor: [[...visibleSet].join(','), revealCount, this.minColor, this.midColor, this.maxColor],
            getRadius: [this.pointSize]
          }
        }));
        }
      } else {
        // Normal coloring or with active selection
        layers.push(pointsLayer);
      }

      if (this.mode === 'gene' && this.showClusterBg && this.clusterPoints.length > 0) {
        layers.unshift(new ScatterplotLayer({
          id: this.bgLayerId,
          data: this.clusterPoints,
          pickable: false,
          radiusUnits: 'pixels',
          getRadius: this.pointSize,
          getPosition: d => [d.x, d.y, 0],
          getFillColor: d => self.cssColorToRgbaArray(d.color, Math.round(self.clusterBgOpacity * 255)),
          updateTriggers: { getFillColor: [this.clusterBgOpacity] }
        }));
      }

      if (this.mode === 'gene' && this.showClusterBg && this.showClusterLabels && this.clusterLabels.length > 0) {
        layers.push(new TextLayer({
          id: this.labelLayerId,
          data: this.clusterLabels,
          getPosition: d => [d.x, d.y, 0],
          getText: d => d.id,
          getColor: [0,0,0,255],
          getSize: 12,
          sizeUnits: 'pixels'
        }));
      }

      this.deck.setProps({ layers });
    },

    onWindowResize() {
      this.requestRedraw();
    },

    initCanvasEvents() {
        const canvas = this.$refs.deckContainer;
        if (!canvas) return;
        canvas.addEventListener("mousedown", this.onMouseDown);
        canvas.addEventListener("mousemove", this.onMouseMove);
        canvas.addEventListener("mouseup", this.onMouseUp);
        canvas.addEventListener("mouseleave", () => (this.isDragging = false));
        // Do not intercept wheel; let deck.gl handle zoom
    },
    
    getMousePos(e) {
        const rect = (this.deck && this.deck.canvas) ? this.deck.canvas.getBoundingClientRect() : this.$refs.deckContainer.getBoundingClientRect();
        return { 
            mouseX: e.clientX - rect.left, 
            mouseY: e.clientY - rect.top 
        };
    },

    onMouseDown(e) {
        const { mouseX, mouseY } = this.getMousePos(e);
        if (this.isSelecting) {
            this.isDragging = true;
            this.startX = mouseX;
            this.startY = mouseY;
            const selectionBox = this.$refs.selectionBox;
            selectionBox.style.left = `${mouseX}px`;
            selectionBox.style.top = `${mouseY}px`;
            selectionBox.style.width = '0px';
            selectionBox.style.height = '0px';
            selectionBox.style.display = 'block';
        } else {
            this.isDragging = true;
            this.startX = mouseX;
            this.startY = mouseY;
        }
    },

    onMouseMove(e) {
        if (!this.isDragging) return;
        const { mouseX, mouseY } = this.getMousePos(e);
        if (this.isSelecting) {
            const width = mouseX - this.startX;
            const height = mouseY - this.startY;
            const selectionBox = this.$refs.selectionBox;
            selectionBox.style.width = `${Math.abs(width)}px`;
            selectionBox.style.height = `${Math.abs(height)}px`;
            selectionBox.style.left = `${width > 0 ? this.startX : mouseX}px`;
            selectionBox.style.top = `${height > 0 ? this.startY : mouseY}px`;
        } else {
            // panning is handled by deck.gl controller when not selecting
        }
    },

    onMouseUp(e) {
        if (this.isSelecting && this.isDragging) {
            this.updateSelectedPoints();
            this.$refs.selectionBox.style.width = '0px';
            this.$refs.selectionBox.style.height = '0px';
            this.$refs.selectionBox.style.display = 'none';
        }
        this.isDragging = false;
    },
    
    onWheel(e) { /* no-op, deck handles */ },

    switchTool(tool) {
        this.isSelecting = tool === 'select';
        if (this.deck) {
          const enable = !this.isSelecting;
          this.deck.setProps({ controller: { dragPan: enable, scrollZoom: enable, doubleClickZoom: enable, touchZoom: enable, touchRotate: false, dragRotate: false } });
        }
    },

    resetView() {
        this.isSelecting = false;
        this.$emit('update:selected-clusters', []);
        this.$emit('update:selected-points', []);
        this.autoAdjustView();
    },
    
    updateSelectedPoints() {
      const selectionBox = this.$refs.selectionBox;
      const rect = selectionBox.getBoundingClientRect();
      const canvasRect = this.deck && this.deck.canvas ? this.deck.canvas.getBoundingClientRect() : this.$refs.deckContainer.getBoundingClientRect();
      const x = rect.left - canvasRect.left;
      const y = rect.top - canvasRect.top;
      const width = rect.right - rect.left;
      const height = rect.bottom - rect.top;
      let selectedPoints = [];
      if (this.deck && width > 0 && height > 0) {
        const layerIds = [this.pointLayerId, this.pointLayerId + '-grey', this.pointLayerId + '-reveal'];
        const infos = this.deck.pickObjects({x, y, width, height, layerIds}) || [];
        selectedPoints = infos.map(info => info.object);
      }

      const selectedClusterIds = [...new Set(selectedPoints.map(p => p.clusterId))];
      this.$emit('update:selected-points', selectedPoints);
      this.$emit('update:selected-clusters', selectedClusterIds);
    },

    getCanvas() {
        return this.deck ? this.deck.canvas : null;
    },
  },
  watch: {
    dataSource: 'loadAndDrawUMAPData',
    dataType: 'loadAndDrawUMAPData',
    mode: 'loadAndDrawUMAPData',
    geneName: 'loadAndDrawUMAPData',
    geneSet: 'loadAndDrawUMAPData',
    showClusterBg: 'loadAndDrawUMAPData',
    pointSize() { this.requestRedraw(); },
    selectedPointsProp() { this.requestRedraw(); },
    visibleClusters() { this.requestRedraw(); },
    geneOpacity() { this.requestRedraw(); },
    clusterBgOpacity() { this.requestRedraw(); },
    showClusterLabels() { this.requestRedraw(); },
    minColor() { this.requestRedraw(); },
    midColor() { this.requestRedraw(); },
    maxColor() { this.requestRedraw(); },
    stepRevealEnabled() { this.requestRedraw(); },
    stepRevealCount() { this.requestRedraw(); },
  },
  mounted() {
    const parent = this.$refs.deckContainer;
    this.deck = new Deck({
      parent,
      views: [new OrthographicView({id: 'ortho'})],
      controller: { dragPan: true, scrollZoom: true, doubleClickZoom: true, touchZoom: true, touchRotate: false, dragRotate: false },
      onViewStateChange: ({viewState}) => { this.viewState = viewState; }
    });
    window.addEventListener("resize", this.onWindowResize);
    this.loadAndDrawUMAPData();
  },
  beforeUnmount() {
    window.removeEventListener("resize", this.onWindowResize);
    const container = this.$refs.deckContainer;
    if (container) {
        container.removeEventListener("mousedown", this.onMouseDown);
        container.removeEventListener("mousemove", this.onMouseMove);
        container.removeEventListener("mouseup", this.onMouseUp);
        container.removeEventListener("mouseleave", () => (this.isDragging = false));
    }
    if (this.deck) { this.deck.finalize(); this.deck = null; }
  }
};
</script>

<style scoped>
.canvas-div {
  position: relative;
  width: 100%;
  height: 100%;
}
.toolbar {
  position: absolute;
  top: 10px;
  left: 10px;
  z-index: 10;
  display: flex;
  gap: 5px;
  background-color: var(--el-toolbar-bg-color, rgba(255, 255, 255, 0.8));
  padding: 5px;
  border-radius: 5px;
}
.tool-button {
  background: var(--el-toolbar-button-bg-color, #fff);
  border: 1px solid var(--el-toolbar-button-border-color, #ccc);
  padding: 5px;
  cursor: pointer;
  border-radius: 3px;
  display: flex;
  align-items: center;
  justify-content: center;
}
.tool-button.active {
  background: var(--el-toolbar-button-active-bg-color, #e0e0e0);
  border-color: var(--el-toolbar-button-active-border-color, #999);
}
.tool-button svg {
  width: 20px;
  height: 20px;
}
.selection {
  border: 1px dashed red;
  position: absolute;
  background-color: rgba(255, 0, 0, 0.1);
  pointer-events: none;
}
.deck-container {
  position: absolute;
  inset: 0;
}
.legend-bar {
  position: absolute;
  bottom: 20px;
  left: 20px;
  width: 200px;
  background-color: rgba(255, 255, 255, 0.8);
  padding: 8px;
  border-radius: 5px;
  box-shadow: 0 2px 4px rgba(0,0,0,0.1);
  z-index: 10;
}
.legend-gradient {
  height: 20px;
  border-radius: 3px;
  background: linear-gradient(to right, rgb(0, 0, 255), rgb(255, 255, 0), rgb(255, 0, 0));
}
.legend-labels {
  display: flex;
  justify-content: space-between;
  font-size: 12px;
  color: #333;
  margin-top: 4px;
}
</style>
