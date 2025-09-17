import axios from 'axios';

const API_BASE_URL =
  import.meta.env.VITE_IP_API_URL ||
  import.meta.env.VITE_API_BASE_URL ||
  (typeof window !== 'undefined' ? `${window.location.origin}/api` : 'http://localhost:3002/api');

const client = axios.create({
  baseURL: API_BASE_URL,
  timeout: 15000
});

export async function fetchIpLogs(params = {}) {
  const response = await client.get('/ip', { params });
  return response.data;
}

export async function exportIpLogs(params = {}) {
  const response = await client.get('/ip/export', {
    params,
    responseType: 'blob'
  });
  return response.data;
}

export async function removeIpAddress(ip) {
  if (!ip) {
    throw new Error('IP address is required');
  }
  return client.delete(`/ip/${encodeURIComponent(ip)}`);
}

export async function clearAllLogs() {
  return client.delete('/ip');
}

export default {
  fetchIpLogs,
  exportIpLogs,
  removeIpAddress,
  clearAllLogs
};
