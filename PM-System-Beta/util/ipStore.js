const fs = require("fs/promises");
const path = require("path");

const LOG_FILE_PATH = path.join(__dirname, "../data/ip_log.json");
const DATE_ONLY_REGEX = /^\d{4}-\d{2}-\d{2}$/;

async function ensureStorage() {
  await fs.mkdir(path.dirname(LOG_FILE_PATH), { recursive: true });
  try {
    await fs.access(LOG_FILE_PATH);
  } catch (err) {
    if (err && err.code !== "ENOENT") {
      throw err;
    }
    await fs.writeFile(LOG_FILE_PATH, "{}", "utf8");
  }
}

function convertLegacyPayload(payload) {
  if (Array.isArray(payload)) {
    return payload.reduce((acc, entry) => {
      if (!entry || !entry.ip) {
        return acc;
      }
      const timestamps = Array.isArray(entry.timestamps)
        ? entry.timestamps
        : entry.timestamp
          ? [entry.timestamp]
          : [];
      if (!acc[entry.ip]) {
        acc[entry.ip] = { timestamps: [] };
      }
      acc[entry.ip].timestamps.push(...timestamps);
      return acc;
    }, {});
  }
  if (payload && typeof payload === "object") {
    return { ...payload };
  }
  return {};
}

function normaliseEntry(entry) {
  if (!entry || typeof entry !== "object") {
    return { count: 0, timestamps: [], countryCode: null, countryName: null };
  }
  const timestamps = Array.isArray(entry.timestamps)
    ? entry.timestamps.filter(Boolean)
    : [];
  const sortedUnique = Array.from(new Set(
    timestamps
      .map((value) => {
        if (typeof value !== "string") {
          return null;
        }
        const trimmed = value.trim();
        return trimmed ? trimmed : null;
      })
      .filter(Boolean)
  )).sort((a, b) => new Date(a).getTime() - new Date(b).getTime());

  let countryCode = null;
  if (typeof entry.countryCode === "string") {
    const trimmedCode = entry.countryCode.trim().toUpperCase();
    if (trimmedCode) {
      countryCode = trimmedCode;
    }
  }

  let countryName = null;
  if (typeof entry.countryName === "string") {
    const trimmedName = entry.countryName.trim();
    if (trimmedName) {
      countryName = trimmedName;
    }
  }

  const rawCountry = typeof entry.country === "string" ? entry.country.trim() : "";
  if (!countryCode && rawCountry && rawCountry.length <= 3) {
    countryCode = rawCountry.toUpperCase();
  }
  if (!countryName && rawCountry && rawCountry.length > 3) {
    countryName = rawCountry;
  }

  return {
    count: sortedUnique.length,
    timestamps: sortedUnique,
    countryCode: countryCode || null,
    countryName: countryName || null
  };
}

async function readRawLogs() {
  await ensureStorage();
  try {
    const data = await fs.readFile(LOG_FILE_PATH, "utf8");
    if (!data || !data.trim()) {
      return {};
    }
    const parsed = JSON.parse(data);
    const legacySafe = convertLegacyPayload(parsed);
    return Object.entries(legacySafe).reduce((acc, [ip, entry]) => {
      const normalised = normaliseEntry(entry);
      if (normalised.timestamps.length === 0) {
        return acc;
    if (normalised.timestamps.length > 0) {
      const payload = {
        count: normalised.count,
        timestamps: normalised.timestamps
      };
      if (normalised.countryCode) {
        payload.countryCode = normalised.countryCode;
      }
      if (normalised.countryName) {
        payload.countryName = normalised.countryName;
      }
      acc[ip] = payload;
    }, {});
  } catch (err) {
    console.error("Failed to read ip_log.json:", err);
    return {};
  }
}

async function writeRawLogs(logs) {
  await ensureStorage();
  const safePayload = Object.entries(logs || {}).reduce((acc, [ip, entry]) => {
    const normalised = normaliseEntry(entry);
    if (normalised.timestamps.length > 0) {
      const payload = {
        count: normalised.count,
        timestamps: normalised.timestamps
      };
      if (normalised.countryCode) {
        payload.countryCode = normalised.countryCode;
      }
      if (normalised.countryName) {
        payload.countryName = normalised.countryName;
      }
      acc[ip] = payload;
    }
    return acc;
  }, {});
  await fs.writeFile(LOG_FILE_PATH, JSON.stringify(safePayload, null, 2), "utf8");
}

function parseBoundary(value, type) {
  if (!value || typeof value !== "string") {
    return null;
  }
  const trimmed = value.trim();
  if (!trimmed) {
    return null;
  }
  const date = new Date(trimmed);
  if (Number.isNaN(date.getTime())) {
    return null;
  }
  const time = date.getTime();
  if (type === "start") {
    return time;
  }
  if (type === "end" && DATE_ONLY_REGEX.test(trimmed)) {
    return time + 24 * 60 * 60 * 1000 - 1;
  }
  return time;
}

function applyFilters(logs, filters = {}) {
  const { search, start, end } = filters;
  const searchValue = typeof search === "string" ? search.trim().toLowerCase() : "";
  const startBoundary = parseBoundary(start, "start");
  const endBoundary = parseBoundary(end, "end");
  const filtered = {};

  Object.entries(logs || {}).forEach(([ip, entry]) => {
    if (searchValue && !ip.toLowerCase().includes(searchValue)) {
      return;
    }
    const normalised = normaliseEntry(entry);
    const filteredTimestamps = normalised.timestamps.filter((timestamp) => {
      const date = new Date(timestamp);
      const time = date.getTime();
      if (Number.isNaN(time)) {
        return false;
      }
      if (startBoundary !== null && time < startBoundary) {
        return false;
      }
      if (endBoundary !== null && time > endBoundary) {
        return false;
      }
      return true;
    });

    if (filteredTimestamps.length === 0) {
      return;
    }

    const sorted = filteredTimestamps.sort((a, b) => new Date(a).getTime() - new Date(b).getTime());

    filtered[ip] = {
      count: sorted.length,
      timestamps: sorted,
      firstVisit: sorted[0],
      lastVisit: sorted[sorted.length - 1],
      countryCode: normalised.countryCode || null,
      countryName: normalised.countryName || null
    };
  });

  return filtered;
}

function buildSummary(logs) {
  const entries = Object.entries(logs || {});
  if (entries.length === 0) {
    return {
      uniqueIps: 0,
      totalVisits: 0,
      firstVisit: null,
      lastVisit: null,
      rangeDays: 0,
      topIps: []
    };
  }

  let earliest = null;
  let latest = null;
  let totalVisits = 0;

  const topIps = entries
    .map(([ip, entry]) => {
      totalVisits += entry.count;
      if (entry.firstVisit && (!earliest || entry.firstVisit < earliest)) {
        earliest = entry.firstVisit;
      }
      if (entry.lastVisit && (!latest || entry.lastVisit > latest)) {
        latest = entry.lastVisit;
      }
      return {
        ip,
        count: entry.count,
        lastVisit: entry.lastVisit || null,
        countryCode: entry.countryCode || null,
        countryName: entry.countryName || null
      };
    })
    .sort((a, b) => {
      if (b.count !== a.count) {
        return b.count - a.count;
      }
      if (a.lastVisit && b.lastVisit) {
        return new Date(b.lastVisit).getTime() - new Date(a.lastVisit).getTime();
      }
      if (a.lastVisit) {
        return -1;
      }
      if (b.lastVisit) {
        return 1;
      }
      return 0;
    })
    .slice(0, 5);

  let rangeDays = 0;
  if (earliest && latest) {
    const diff = new Date(latest).getTime() - new Date(earliest).getTime();
    rangeDays = diff >= 0 ? Math.max(1, Math.round(diff / (24 * 60 * 60 * 1000)) + 1) : 0;
  }

  return {
    uniqueIps: entries.length,
    totalVisits,
    firstVisit: earliest,
    lastVisit: latest,
    rangeDays,
    topIps
  };
}

async function getFilteredLogs(filters = {}) {
  const rawLogs = await readRawLogs();
  const logs = applyFilters(rawLogs, filters);
  return {
    logs,
    summary: buildSummary(logs)
  };
}

async function recordVisit(ip, metadata = {}) {
  if (!ip || typeof ip !== "string") {
    return;
  }

  let timestamp = new Date().toISOString();
  let countryCode = null;
  let countryName = null;

  if (typeof metadata === "string") {
    timestamp = metadata;
    metadata = {};
  }

  if (metadata && typeof metadata === "object") {
    if (typeof metadata.timestamp === "string" && metadata.timestamp.trim()) {
      timestamp = metadata.timestamp.trim();
    }
    if (typeof metadata.countryCode === "string") {
      const trimmedCode = metadata.countryCode.trim().toUpperCase();
      if (trimmedCode) {
        countryCode = trimmedCode;
      }
    }
    if (typeof metadata.countryName === "string") {
      const trimmedName = metadata.countryName.trim();
      if (trimmedName) {
        countryName = trimmedName;
      }
    }
  }

  const rawLogs = await readRawLogs();
  const existing = normaliseEntry(rawLogs[ip]);
  const timestampSet = new Set(existing.timestamps);
  if (typeof timestamp === "string" && timestamp.trim()) {
    timestampSet.add(timestamp.trim());
  }

  const nextTimestamps = Array.from(timestampSet).sort(
    (a, b) => new Date(a).getTime() - new Date(b).getTime()
  );

  const payload = {
    count: nextTimestamps.length,
    timestamps: nextTimestamps
  };

  const finalCountryCode = countryCode || existing.countryCode || null;
  const finalCountryName = countryName || existing.countryName || null;
  if (finalCountryCode) {
    payload.countryCode = finalCountryCode;
  }
  if (finalCountryName) {
    payload.countryName = finalCountryName;
  }

  rawLogs[ip] = payload;
  await writeRawLogs(rawLogs);
}

async function removeIp(ip) {
  if (!ip || typeof ip !== "string") {
    return false;
  }
  const rawLogs = await readRawLogs();
  if (!Object.prototype.hasOwnProperty.call(rawLogs, ip)) {
    return false;
  }
  delete rawLogs[ip];
  await writeRawLogs(rawLogs);
  return true;
}

async function clearLogs() {
  await writeRawLogs({});
}

function escapeForCsv(value) {
  if (value === null || value === undefined) {
    return "";
  }
  const stringValue = String(value);
  if (stringValue === "") {
    return "";
  }
  if (/[",\n]/.test(stringValue)) {
    return '"' + stringValue.replace(/"/g, '""') + '"';
  }
  return stringValue;
}

function convertToCsv(logs) {
  const header = ["IP Address", "Country", "Visit Count", "First Visit", "Last Visit", "Timestamps"];
  const rows = [header.join(",")];

  Object.entries(logs || {})
    .sort((a, b) => {
      const countDiff = (b[1].count || 0) - (a[1].count || 0);
      if (countDiff !== 0) {
        return countDiff;
      }
      const lastA = a[1].lastVisit ? new Date(a[1].lastVisit).getTime() : 0;
      const lastB = b[1].lastVisit ? new Date(b[1].lastVisit).getTime() : 0;
      return lastB - lastA;
    })
    .forEach(([ip, entry]) => {
      const values = [
        escapeForCsv(ip),
        escapeForCsv(entry.countryName || entry.countryCode || ""),
        escapeForCsv(entry.count ?? entry.timestamps?.length ?? 0),
        escapeForCsv(entry.firstVisit || ""),
        escapeForCsv(entry.lastVisit || ""),
        escapeForCsv((entry.timestamps || []).join(" | "))
      ];
      rows.push(values.join(","));
    });

  return rows.join("\n");
}

module.exports = {
  applyFilters,
  buildSummary,
  clearLogs,
  convertToCsv,
  getFilteredLogs,
  readRawLogs,
  recordVisit,
  removeIp,
  writeRawLogs
};
