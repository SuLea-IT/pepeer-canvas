const express = require("express");
const {
  clearLogs,
  convertToCsv,
  getFilteredLogs,
  removeIp
} = require("../util/ipStore");

const router = express.Router();

function parseFilters(query = {}) {
  const pick = (value) => {
    if (Array.isArray(value)) {
      return value[0];
    }
    return value;
  };

  const filters = {};

  const search = pick(query.search);
  if (typeof search === "string" && search.trim()) {
    filters.search = search.trim();
  }

  const start = pick(query.start);
  if (typeof start === "string" && start.trim()) {
    filters.start = start.trim();
  }

  const end = pick(query.end);
  if (typeof end === "string" && end.trim()) {
    filters.end = end.trim();
  }

  return filters;
}

router.get("/", async (req, res) => {
  try {
    const filters = parseFilters(req.query);
    const result = await getFilteredLogs(filters);
    res.json(result);
  } catch (error) {
    console.error("Failed to fetch IP logs:", error);
    res.status(500).json({ message: "Failed to read IP logs." });
  }
});

router.get("/export", async (req, res) => {
  try {
    const filters = parseFilters(req.query);
    const { logs } = await getFilteredLogs(filters);
    const csv = convertToCsv(logs);
    res.setHeader("Content-Type", "text/csv; charset=utf-8");
    res.setHeader("Content-Disposition", "attachment; filename=ip-access-log.csv");
    res.send(csv);
  } catch (error) {
    console.error("Failed to export IP logs:", error);
    res.status(500).json({ message: "Failed to export IP logs." });
  }
});

router.delete("/", async (_req, res) => {
  try {
    await clearLogs();
    res.json({ success: true });
  } catch (error) {
    console.error("Failed to clear IP logs:", error);
    res.status(500).json({ message: "Failed to clear IP logs." });
  }
});

router.delete("/:ip", async (req, res) => {
  const targetIp = decodeURIComponent(req.params.ip || "");

  if (!targetIp) {
    return res.status(400).json({ message: "IP address is required." });
  }

  try {
    const removed = await removeIp(targetIp);
    if (!removed) {
      return res.status(404).json({ message: "IP address not found." });
    }
    res.json({ success: true });
  } catch (error) {
    console.error(`Failed to remove IP ${targetIp}:`, error);
    res.status(500).json({ message: "Failed to remove IP address." });
  }
});

module.exports = router;
