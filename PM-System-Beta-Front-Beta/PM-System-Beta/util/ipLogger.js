const { recordVisit } = require("./ipStore");
const { lookupCountry } = require("./ipGeo");

const LOCAL_IDENTIFIERS = new Set([
  "127.0.0.1",
  "0.0.0.0",
  "::1",
  "0:0:0:0:0:0:0:1",
  "localhost"
]);

function normaliseIp(rawValue) {
  if (!rawValue) {
    return { ip: null, isLocal: false };
  }

  let value = String(rawValue).split(",")[0].trim();
  if (!value) {
    return { ip: null, isLocal: false };
  }

  if (value.includes(" ")) {
    value = value.split(" ")[0].trim();
  }

  if (value.startsWith("[") && value.endsWith("]")) {
    value = value.slice(1, -1);
  }

  if (value.startsWith("::ffff:")) {
    value = value.substring(7);
  }

  if (value.includes("%")) {
    value = value.split("%", 1)[0];
  }

  if (/^\d+\.\d+\.\d+\.\d+:\d+$/.test(value)) {
    value = value.split(":", 1)[0];
  }

  const lower = value.toLowerCase();
  if (LOCAL_IDENTIFIERS.has(lower) || LOCAL_IDENTIFIERS.has(value)) {
    return { ip: lower === "::1" ? "::1" : value, isLocal: true };
  }

  return { ip: value, isLocal: false };
}

function extractClientIp(req) {
  const forwardedHeader = req.headers["x-forwarded-for"];
  const forwarded = Array.isArray(forwardedHeader)
    ? forwardedHeader
    : typeof forwardedHeader === "string"
      ? forwardedHeader.split(",")
      : [];

  const candidates = [
    ...forwarded,
    req.headers["x-real-ip"],
    req.socket?.remoteAddress,
    req.connection?.remoteAddress,
    req.ip
  ].filter(Boolean);

  for (const candidate of candidates) {
    const result = normaliseIp(candidate);
    if (result.ip) {
      return result;
    }
  }

  return { ip: null, isLocal: false };
}

const ipLogger = async (req, res, next) => {
  if (req.method === "OPTIONS") {
    return next();
  }

  const { ip, isLocal } = extractClientIp(req);

  if (!ip || isLocal) {
    return next();
  }

  const location = lookupCountry(ip);
  const visitMetadata = {};
  if (location && location.countryCode) {
    visitMetadata.countryCode = location.countryCode;
  }
  if (location && location.countryName) {
    visitMetadata.countryName = location.countryName;
  }

  try {
    await recordVisit(ip, visitMetadata);
  } catch (error) {
    console.error(`Failed to persist IP visit for ${ip}:`, error);
  }

  return next();
};

module.exports = ipLogger;
