// Health feature module - backend status updates and polling
(function(global) {
  async function updateBackendStatus() {
    const statusElement = document.getElementById('backend-status');
    if (!statusElement) return;
    
    const indicator = statusElement.querySelector('.status-indicator');
    const statusText = statusElement.querySelector('.status-text');
    if (!indicator || !statusText) return;
    
    indicator.className = 'fas fa-circle status-indicator loading';
    statusText.textContent = 'Checking backend status...';
    
    try {
      const data = await ApiClient.getHealth();
      if (data && data.model_loaded) {
        indicator.className = 'fas fa-circle status-indicator ready';
        statusText.textContent = `Backend ready - Model loaded (uptime: ${data.uptime_seconds || 'unknown'}s)`;
      } else if (data.status === 'loading') {
        indicator.className = 'fas fa-circle status-indicator loading';
        statusText.textContent = `Model loading in background... (uptime: ${data.uptime_seconds || 'unknown'}s)`;
      } else if (data.status === 'initializing') {
        indicator.className = 'fas fa-circle status-indicator initializing';
        statusText.textContent = `Model initializing... (uptime: ${data.uptime_seconds || 'unknown'}s)`;
      } else if (data.error) {
        indicator.className = 'fas fa-circle status-indicator error';
        statusText.textContent = 'Model initialization failed';
      } else {
        indicator.className = 'fas fa-circle status-indicator error';
        statusText.textContent = 'Backend unavailable';
      }
    } catch (error) {
      console.warn('Backend status check failed:', error);
      indicator.className = 'fas fa-circle status-indicator error';
      if (error.name === 'AbortError') {
        statusText.textContent = 'Backend connection timeout - still starting up';
      } else {
        statusText.textContent = 'Cannot connect to backend';
      }
    }
  }

  function startPolling(intervalMs = 10000) {
    setTimeout(() => { updateBackendStatus(); }, 100);
    const statusCheckInterval = setInterval(async () => {
      try { await updateBackendStatus(); } catch (e) { updateBackendStatus(); }
    }, intervalMs);
    return statusCheckInterval;
  }

  global.Health = { updateBackendStatus, startPolling };
})(window);


