// Lightweight API client with timeout and unified error handling
(function(global) {
  const DEFAULT_TIMEOUT_MS = 20000;

  async function fetchJson(url, opts = {}, timeoutMs = DEFAULT_TIMEOUT_MS) {
    const controller = new AbortController();
    const id = setTimeout(() => controller.abort(), timeoutMs);
    try {
      const res = await fetch(url, { ...opts, signal: controller.signal });
      const data = await res.json().catch(() => ({}));
      if (!res.ok) {
        const msg = data && data.error ? data.error : `Request failed (${res.status})`;
        throw new Error(msg);
      }
      return data;
    } finally {
      clearTimeout(id);
    }
  }

  async function getHealth() {
    return await fetchJson('/health', { method: 'GET', headers: { 'Accept': 'application/json' } }, 5000);
  }

  async function postPredict(body) {
    return await fetchJson('/predict', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body)
    });
  }

  async function postSmilesSearch(body) {
    return await fetchJson('/smiles-search', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body)
    });
  }

  async function postAnalyze(body) {
    return await fetchJson('/analyze', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body)
    }, 30000);
  }

  async function postSecondaryRetrieval(body) {
    return await fetchJson('/secondary-retrieval', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body)
    });
  }

  global.ApiClient = {
    fetchJson,
    getHealth,
    postPredict,
    postSmilesSearch,
    postAnalyze,
    postSecondaryRetrieval
  };
})(window);


