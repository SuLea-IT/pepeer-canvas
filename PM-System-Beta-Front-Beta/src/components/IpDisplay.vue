<template>
  <div class="min-h-screen bg-gray-100 dark:bg-gray-900 py-10">
    <div class="max-w-7xl mx-auto px-4 space-y-8">
      <header class="flex flex-col md:flex-row md:items-center md:justify-between gap-4">
        <div>
          <h1 class="text-3xl font-bold text-gray-900 dark:text-white">IP Access Management</h1>
          <p class="text-sm text-gray-500 dark:text-gray-400">Track visitors, filter logs, and keep your portal secure.</p>
        </div>
        <div class="flex flex-wrap items-center gap-3">
          <label class="flex items-center text-sm text-gray-600 dark:text-gray-300">
            <input v-model="autoRefresh" type="checkbox" class="mr-2 rounded border-gray-300" />
            Auto refresh (60s)
          </label>
          <button
            @click="loadLogs"
            :disabled="loading"
            class="inline-flex items-center px-4 py-2 rounded-md bg-white dark:bg-gray-800 border border-gray-300 dark:border-gray-700 text-gray-700 dark:text-gray-200 shadow-sm hover:bg-gray-50 disabled:opacity-60"
          >
            <i class="fa fa-sync mr-2" aria-hidden="true"></i>
            Refresh
          </button>
          <button
            @click="downloadCsv"
            :disabled="loading || !hasData"
            class="inline-flex items-center px-4 py-2 rounded-md bg-primary text-white shadow-sm hover:bg-primary/90 disabled:opacity-60"
          >
            <i class="fa fa-download mr-2" aria-hidden="true"></i>
            Export CSV
          </button>
          <button
            @click="confirmClearAll"
            :disabled="loading || !hasData"
            class="inline-flex items-center px-4 py-2 rounded-md bg-red-600 text-white shadow-sm hover:bg-red-500 disabled:opacity-60"
          >
            <i class="fa fa-trash mr-2" aria-hidden="true"></i>
            Clear All Logs
          </button>
        </div>
      </header>

      <section class="bg-white dark:bg-gray-800 shadow rounded-lg p-6 space-y-6">
        <form @submit.prevent="loadLogs" class="grid grid-cols-1 md:grid-cols-[2fr,repeat(2,1fr),auto] gap-4">
          <div>
            <label class="block text-xs font-semibold text-gray-500 dark:text-gray-400 uppercase tracking-wide mb-1">Search</label>
            <input
              v-model="filters.search"
              type="text"
              placeholder="IP substring (e.g. 192.168)"
              class="w-full rounded-md border border-gray-300 dark:border-gray-700 bg-white dark:bg-gray-900 text-gray-800 dark:text-gray-100 px-3 py-2 focus:outline-none focus:ring-2 focus:ring-primary"
            />
          </div>
          <div>
            <label class="block text-xs font-semibold text-gray-500 dark:text-gray-400 uppercase tracking-wide mb-1">Start</label>
            <input
              v-model="filters.start"
              type="datetime-local"
              class="w-full rounded-md border border-gray-300 dark:border-gray-700 bg-white dark:bg-gray-900 text-gray-800 dark:text-gray-100 px-3 py-2 focus:outline-none focus:ring-2 focus:ring-primary"
            />
          </div>
          <div>
            <label class="block text-xs font-semibold text-gray-500 dark:text-gray-400 uppercase tracking-wide mb-1">End</label>
            <input
              v-model="filters.end"
              type="datetime-local"
              class="w-full rounded-md border border-gray-300 dark:border-gray-700 bg-white dark:bg-gray-900 text-gray-800 dark:text-gray-100 px-3 py-2 focus:outline-none focus:ring-2 focus:ring-primary"
            />
          </div>
          <div class="flex items-end gap-2">
            <button
              type="submit"
              :disabled="loading"
              class="flex-1 inline-flex items-center justify-center px-4 py-2 rounded-md bg-primary text-white shadow hover:bg-primary/90 disabled:opacity-60"
            >
              <i class="fa fa-search mr-2" aria-hidden="true"></i>
              Apply
            </button>
            <button
              type="button"
              @click="resetFilters"
              :disabled="loading"
              class="inline-flex items-center justify-center px-4 py-2 rounded-md border border-gray-300 dark:border-gray-700 text-gray-700 dark:text-gray-200 hover:bg-gray-50 dark:hover:bg-gray-700"
            >
              Reset
            </button>
          </div>
        </form>

        <div v-if="loading" class="flex items-center text-sm text-gray-600 dark:text-gray-300">
          <i class="fa fa-spinner fa-spin mr-2" aria-hidden="true"></i>
          Loading IP logs...
        </div>
        <div v-if="error" class="flex items-center text-sm text-red-500">
          <i class="fa fa-exclamation-triangle mr-2" aria-hidden="true"></i>
          {{ error }}
        </div>
        <div v-if="hasData" class="grid grid-cols-1 md:grid-cols-4 gap-4">
          <div class="p-4 rounded-lg border border-gray-200 dark:border-gray-700 bg-gray-50 dark:bg-gray-900">
            <p class="text-xs uppercase tracking-wide text-gray-500 dark:text-gray-400 mb-1">Unique IPs</p>
            <p class="text-2xl font-semibold text-gray-900 dark:text-white">{{ summary.uniqueIps }}</p>
          </div>
          <div class="p-4 rounded-lg border border-gray-200 dark:border-gray-700 bg-gray-50 dark:bg-gray-900">
            <p class="text-xs uppercase tracking-wide text-gray-500 dark:text-gray-400 mb-1">Total Visits</p>
            <p class="text-2xl font-semibold text-gray-900 dark:text-white">{{ summary.totalVisits }}</p>
          </div>
          <div class="p-4 rounded-lg border border-gray-200 dark:border-gray-700 bg-gray-50 dark:bg-gray-900">
            <p class="text-xs uppercase tracking-wide text-gray-500 dark:text-gray-400 mb-1">First Visit</p>
            <p class="text-sm text-gray-700 dark:text-gray-200">{{ summary.firstVisit ? formatDate(summary.firstVisit) : '?' }}</p>
          </div>
          <div class="p-4 rounded-lg border border-gray-200 dark:border-gray-700 bg-gray-50 dark:bg-gray-900">
            <p class="text-xs uppercase tracking-wide text-gray-500 dark:text-gray-400 mb-1">Last Visit</p>
            <p class="text-sm text-gray-700 dark:text-gray-200">{{ summary.lastVisit ? formatDate(summary.lastVisit) : '?' }}</p>
          </div>
        </div>
        <div v-if="hasData && summary.topIps && summary.topIps.length" class="bg-gray-50 dark:bg-gray-900 rounded-lg p-4">
          <p class="text-xs uppercase tracking-wide text-gray-500 dark:text-gray-400 mb-3">Top Visitors</p>
          <ul class="space-y-2 text-sm text-gray-700 dark:text-gray-200">
            <li v-for="item in summary.topIps" :key="item.ip" class="flex justify-between items-start">
              <div class="flex flex-col">
                <span class="font-mono">{{ item.ip }}</span>
                <span class="text-xs text-gray-500 dark:text-gray-400">{{ formatCountry(item.countryName, item.countryCode) }}</span>
              </div>
              <span>{{ item.count }} visits</span>
            </li>
          </ul>
        </div>
        <div v-if="lastUpdated" class="text-xs text-gray-500 dark:text-gray-400">
          Last updated: {{ formatDate(lastUpdated) }}
        </div>

        <div v-if="!loading && !hasData" class="text-center py-12 text-gray-500 dark:text-gray-400">
          <i class="fa fa-info-circle text-xl mb-2" aria-hidden="true"></i>
          <p>No IP logs match the current filters.</p>
        </div>

        <div v-if="hasData" class="overflow-x-auto">
          <table class="min-w-full divide-y divide-gray-200 dark:divide-gray-700">
            <thead class="bg-gray-100 dark:bg-gray-700/50">
              <tr>
                <th class="px-4 py-3 text-left text-xs font-semibold text-gray-500 dark:text-gray-300 uppercase tracking-wide">IP Address</th>
                <th class="px-4 py-3 text-left text-xs font-semibold text-gray-500 dark:text-gray-300 uppercase tracking-wide">Country</th>
                <th class="px-4 py-3 text-left text-xs font-semibold text-gray-500 dark:text-gray-300 uppercase tracking-wide">Visits</th>
                <th class="px-4 py-3 text-left text-xs font-semibold text-gray-500 dark:text-gray-300 uppercase tracking-wide">First Visit</th>
                <th class="px-4 py-3 text-left text-xs font-semibold text-gray-500 dark:text-gray-300 uppercase tracking-wide">Last Visit</th>
                <th class="px-4 py-3 text-right text-xs font-semibold text-gray-500 dark:text-gray-300 uppercase tracking-wide">Actions</th>
              </tr>
            </thead>
            <tbody class="divide-y divide-gray-200 dark:divide-gray-700">
              <template v-for="row in rows" :key="row.ip">
                <tr class="hover:bg-gray-50 dark:hover:bg-gray-700/30">
                  <td class="px-4 py-3 font-mono text-sm text-gray-900 dark:text-white">
                    <button @click="toggle(row.ip)" class="flex items-center gap-2 focus:outline-none">
                      <i :class="['fa', expanded.includes(row.ip) ? 'fa-chevron-up' : 'fa-chevron-down', 'text-xs text-gray-400']" aria-hidden="true"></i>
                      {{ row.ip }}
                    </button>
                  </td>
                  <td class="px-4 py-3 text-sm text-gray-600 dark:text-gray-300">{{ formatCountry(row.countryName, row.countryCode) }}</td>
                  <td class="px-4 py-3 text-sm text-gray-600 dark:text-gray-300">{{ row.count }}</td>
                  <td class="px-4 py-3 text-sm text-gray-600 dark:text-gray-300">{{ row.firstVisit ? formatDate(row.firstVisit) : '?' }}</td>
                  <td class="px-4 py-3 text-sm text-gray-600 dark:text-gray-300">{{ row.lastVisit ? formatDate(row.lastVisit) : '?' }}</td>
                  <td class="px-4 py-3 text-right">
                    <button
                      @click.stop="confirmRemove(row.ip)"
                      class="inline-flex items-center px-3 py-1.5 rounded-md bg-red-50 text-red-600 border border-red-200 hover:bg-red-100"
                    >
                      <i class="fa fa-times mr-1" aria-hidden="true"></i>
                      Remove
                    </button>
                  </td>
                </tr>
                <tr v-if="expanded.includes(row.ip)">
                  <td colspan="6" class="px-6 pb-4 pt-0">
                    <div class="bg-gray-50 dark:bg-gray-900 rounded-lg p-4">
                      <p class="text-xs uppercase tracking-wide text-gray-500 dark:text-gray-400 mb-2">
                        Visit Timeline ({{ row.timestamps.length }} entries)
                      </p>
                      <div class="max-h-48 overflow-y-auto pr-2">
                        <ul class="space-y-1 text-sm text-gray-700 dark:text-gray-200 font-mono">
                          <li v-for="timestamp in row.timestamps" :key="timestamp">{{ formatDate(timestamp) }}</li>
                        </ul>
                      </div>
                    </div>
                  </td>
                </tr>
              </template>
            </tbody>
          </table>
        </div>
      </section>
    </div>
  </div>
</template>

<script>
import {
  clearAllLogs,
  exportIpLogs,
  fetchIpLogs,
  removeIpAddress
} from '../services/ipApi';

const defaultSummary = () => ({
  uniqueIps: 0,
  totalVisits: 0,
  firstVisit: null,
  lastVisit: null,
  rangeDays: 0,
  topIps: []
});

export default {
  data() {
    return {
      filters: {
        search: '',
        start: '',
        end: ''
      },
      logs: {},
      summary: defaultSummary(),
      expanded: [],
      loading: false,
      error: null,
      autoRefresh: false,
      refreshTimer: null,
      lastUpdated: null
    };
  },
  computed: {
    hasData() {
      return Object.keys(this.logs).length > 0;
    },
    rows() {
      return Object.entries(this.logs)
        .map(([ip, entry]) => ({
          ip,
          count: entry.count || entry.timestamps?.length || 0,
          firstVisit: entry.firstVisit || entry.timestamps?.[0] || null,
          lastVisit: entry.lastVisit || entry.timestamps?.[entry.timestamps.length - 1] || null,
          timestamps: entry.timestamps || [],
          countryName: entry.countryName || null,
          countryCode: entry.countryCode || null
        }))
        .sort((a, b) => {
          if (b.count !== a.count) {
            return b.count - a.count;
          }
          const lastA = a.lastVisit ? new Date(a.lastVisit).getTime() : 0;
          const lastB = b.lastVisit ? new Date(b.lastVisit).getTime() : 0;
          return lastB - lastA;
        });
    }
  },
  watch: {
    autoRefresh(value) {
      if (value) {
        this.startAutoRefresh();
      } else {
        this.stopAutoRefresh();
      }
    }
  },
  created() {
    this.loadLogs();
  },
  beforeUnmount() {
    this.stopAutoRefresh();
  },
  methods: {
    async loadLogs() {
      this.loading = true;
      this.error = null;
      try {
        const payload = await fetchIpLogs(this.normalisedFilters());
        this.logs = payload.logs || {};
        this.summary = { ...defaultSummary(), ...(payload.summary || {}) };
        this.lastUpdated = new Date().toISOString();
        // Remove expansion state for IPs no longer present
        this.expanded = this.expanded.filter((ip) => Boolean(this.logs[ip]));
      } catch (err) {
        console.error('Failed to load IP logs:', err);
        this.error = 'Unable to load IP logs. Please try again.';
      } finally {
        this.loading = false;
      }
    },
    normalisedFilters() {
      const params = {};
      if (this.filters.search) {
        params.search = this.filters.search;
      }
      if (this.filters.start) {
        params.start = this.filters.start;
      }
      if (this.filters.end) {
        params.end = this.filters.end;
      }
      return params;
    },
    resetFilters() {
      this.filters = {
        search: '',
        start: '',
        end: ''
      };
      this.loadLogs();
    },
    toggle(ip) {
      const index = this.expanded.indexOf(ip);
      if (index > -1) {
        this.expanded.splice(index, 1);
      } else {
        this.expanded.push(ip);
      }
    },
    startAutoRefresh() {
      this.stopAutoRefresh();
      this.refreshTimer = window.setInterval(() => {
        this.loadLogs();
      }, 60000);
    },
    stopAutoRefresh() {
      if (this.refreshTimer) {
        window.clearInterval(this.refreshTimer);
        this.refreshTimer = null;
      }
    },
    formatCountry(name, code) {
      const rawName = typeof name === "string" ? name.trim() : "";
      const rawCode = typeof code === "string" ? code.trim() : "";
      if (rawName && rawCode && rawName.toUpperCase() !== rawCode.toUpperCase()) {
        return `${rawName} (${rawCode})`;
      }
      if (rawName) {
        return rawName;
      }
      if (rawCode) {
        return rawCode;
      }
      return 'Unknown';
    },
    formatDate(value) {
      if (!value) {
        return '';
      }
      const date = new Date(value);
      if (Number.isNaN(date.getTime())) {
        return value;
      }
      return date.toLocaleString();
    },
    async confirmRemove(ip) {
      const confirmed = window.confirm(`Remove ${ip} from the log?`);
      if (!confirmed) {
        return;
      }
      try {
        await removeIpAddress(ip);
        await this.loadLogs();
      } catch (err) {
        console.error(`Failed to remove IP ${ip}:`, err);
        this.error = `Failed to remove IP ${ip}.`;
      }
    },
    async confirmClearAll() {
      const confirmed = window.confirm('This will delete all stored IP logs. Continue?');
      if (!confirmed) {
        return;
      }
      try {
        await clearAllLogs();
        await this.loadLogs();
      } catch (err) {
        console.error('Failed to clear IP logs:', err);
        this.error = 'Failed to clear IP logs.';
      }
    },
    async downloadCsv() {
      try {
        const blobData = await exportIpLogs(this.normalisedFilters());
        const blob = new Blob([blobData], { type: 'text/csv;charset=utf-8;' });
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.setAttribute('download', 'ip-access-log.csv');
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        URL.revokeObjectURL(url);
      } catch (err) {
        console.error('Failed to export IP logs:', err);
        this.error = 'Failed to export IP logs.';
      }
    }
  }
};
</script>
