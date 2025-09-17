<template>
  <div ref="chart" class="proportion-chart"></div>
</template>

<script>
import * as echarts from 'echarts';
import { isDark } from '../theme/composables/dark';
import { watch } from 'vue';

export default {
  name: 'ClusterProportionChart',
  props: {
    chartData: {
      type: Object,
      required: true,
    },
    visibleClusters: {
      type: Array,
      default: () => [],
    },
    chartType: {
      type: String,
      default: 'pie',
    }
  },
  data() {
    return {
      chartInstance: null,
    };
  },
  mounted() {
    this.chartInstance = echarts.init(this.$refs.chart);
    this.chartInstance.on('legendselectchanged', this.handleLegendChange);
    this.updateChart();
    watch(isDark, () => this.updateChart());
  },
  watch: {
    chartData: {
      handler() {
        this.updateChart();
      },
      deep: true,
    },
    visibleClusters: {
      handler() {
        this.updateChart();
      },
      deep: true,
    },
    chartType() {
      this.updateChart();
    }
  },
  methods: {
    handleLegendChange(params) {
      const visibleIds = [];
      const nameToIdMap = new Map(this.chartData.data.map(d => [d.name, d.id]));
      for (const name in params.selected) {
        if (params.selected[name]) {
          const id = nameToIdMap.get(name);
          if (id !== undefined) {
            visibleIds.push(id);
          }
        }
      }
      this.$emit('update:visible-clusters', visibleIds);
    },
    updateChart() {
      if (!this.chartInstance || !this.chartData || !this.chartData.data || !this.chartData.data.length) {
        return;
      }

      const { data, total } = this.chartData;

      const legendSelected = {};
      data.forEach(item => {
        legendSelected[item.name] = this.visibleClusters.includes(item.id);
      });

      const textColor = isDark.value ? '#fff' : '#333';

      let option;
      if (this.chartType === 'pie') {
        option = {
          title: {
            text: 'Cluster Proportion',
            subtext: ``,
            left: 'center',
            textStyle: {
              color: textColor,
            },
          },
          tooltip: {
            trigger: 'item',
            formatter: '{a} <br/>{b} : {c} ({d}%)',
          },
          legend: {
            type: 'scroll',
            orient: 'vertical',
            left: 'left',
            top: 20,
            bottom: 20,
            selected: legendSelected,
            textStyle: {
              color: textColor,
            },
          },
          series: [
            {
              name: 'Proportion',
              type: 'pie',
              radius: '50%',
              center: ['65%', '50%'],
              data: data,
              emphasis: {
                itemStyle: {
                  shadowBlur: 10,
                  shadowOffsetX: 0,
                  shadowColor: 'rgba(0, 0, 0, 0.5)',
                },
              },
              label: {
                show: false,
              },
              labelLine: {
                show: false,
              }
            },
          ],
        };
      } else { // Bar chart
        option = {
          title: {
            text: 'Cluster Counts',
            subtext: ``,
            left: 'center',
            textStyle: {
              color: textColor,
            },
          },
          tooltip: {
            trigger: 'axis',
            axisPointer: {
              type: 'shadow'
            }
          },
          legend: {
            show: false
          },
          xAxis: {
            type: 'category',
            data: data.map(d => d.name),
            axisLabel: {
              rotate: 45,
              color: textColor,
            }
          },
          yAxis: {
            type: 'value',
            axisLabel: {
              color: textColor,
            },
          },
          series: [
            {
              name: 'Count',
              type: 'bar',
              data: data,
            },
          ],
        };
      }

      this.chartInstance.setOption(option, true);
    },
  },
  beforeDestroy() {
    if (this.chartInstance) {
      this.chartInstance.off('legendselectchanged', this.handleLegendChange);
      this.chartInstance.dispose();
    }
  },
};
</script>

<style scoped>
.proportion-chart {
  width: 100%;
  height: 400px;
}
</style>
