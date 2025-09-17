const geoip = require('geoip-lite');

let regionFormatter;

function getRegionFormatter() {
  if (regionFormatter !== undefined) {
    return regionFormatter;
  }
  if (typeof Intl !== 'undefined' && typeof Intl.DisplayNames === 'function') {
    try {
      regionFormatter = new Intl.DisplayNames(['en'], { type: 'region' });
      return regionFormatter;
    } catch (error) {
      console.warn('Failed to initialise region formatter:', error);
    }
  }
  regionFormatter = null;
  return regionFormatter;
}

function resolveCountryNameFromCode(code) {
  if (!code || typeof code !== 'string') {
    return null;
  }
  const formatter = getRegionFormatter();
  if (formatter) {
    try {
      const displayName = formatter.of(code);
      if (displayName && displayName.toUpperCase() !== code.toUpperCase()) {
        return displayName;
      }
    } catch (_error) {
      // Ignore - fall back to the ISO country code
    }
  }
  return null;
}

function lookupCountry(ip) {
  if (!ip || typeof ip !== 'string') {
    return null;
  }
  try {
    const result = geoip.lookup(ip);
    if (!result || !result.country) {
      return null;
    }
    const countryCode = String(result.country).toUpperCase();
    const countryName = resolveCountryNameFromCode(countryCode);
    return {
      countryCode,
      countryName: countryName || null
    };
  } catch (error) {
    console.error('Failed to lookup country for IP ' + ip + ':', error);
    return null;
  }
}

module.exports = {
  lookupCountry
};
