// Error UI module: renders error messages safely
(function(global) {
  function showMessage(message, type = 'info') {
    const messageDiv = document.createElement('div');
    messageDiv.className = `${type}-message`;
    messageDiv.textContent = message;
    const container = document.querySelector('.container');
    if (container) {
      container.insertBefore(messageDiv, container.firstChild);
      setTimeout(() => { messageDiv.remove(); }, 5000);
    }
  }

  function showDetailedError(userMessage, error) {
    const resultsGrid = document.getElementById('results-grid');
    if (!resultsGrid) return;

    const errorDiv = document.createElement('div');
    errorDiv.className = 'error-details';

    const header = document.createElement('div');
    header.className = 'error-header';
    const icon = document.createElement('i');
    icon.className = 'fas fa-exclamation-triangle';
    const h3 = document.createElement('h3');
    h3.textContent = 'Prediction Failed';
    header.appendChild(icon);
    header.appendChild(h3);
    errorDiv.appendChild(header);

    const msg = document.createElement('div');
    msg.className = 'error-message';
    const p = document.createElement('p');
    const strong = document.createElement('strong');
    strong.textContent = 'Error: ';
    p.appendChild(strong);
    const span = document.createElement('span');
    span.textContent = userMessage || 'Unexpected error';
    p.appendChild(span);
    msg.appendChild(p);
    errorDiv.appendChild(msg);

    const troubleshooting = document.createElement('div');
    troubleshooting.className = 'error-troubleshooting';
    const h4 = document.createElement('h4');
    h4.textContent = 'Troubleshooting Steps:';
    troubleshooting.appendChild(h4);
    const ul = document.createElement('ul');
    ['Check if the Docker container is running: docker ps',
     'Verify the data folder is mounted correctly',
     'Check server logs: docker-compose logs -f',
     'Restart the container: ./restart.sh',
     'Ensure all required files are in the data/ directory']
      .forEach(text => { const li = document.createElement('li'); li.textContent = text; ul.appendChild(li); });
    troubleshooting.appendChild(ul);
    errorDiv.appendChild(troubleshooting);

    const actions = document.createElement('div');
    actions.className = 'error-actions';
    const retryBtn = document.createElement('button');
    retryBtn.className = 'btn btn-primary';
    retryBtn.innerHTML = '<i class="fas fa-redo"></i> Retry Prediction';
    retryBtn.addEventListener('click', function() { if (typeof retryPrediction === 'function') retryPrediction(); });
    const healthBtn = document.createElement('button');
    healthBtn.className = 'btn btn-secondary';
    healthBtn.innerHTML = '<i class="fas fa-heartbeat"></i> Check Backend Health';
    healthBtn.addEventListener('click', function() { if (typeof updateBackendStatus === 'function') updateBackendStatus(); });
    actions.appendChild(retryBtn);
    actions.appendChild(healthBtn);
    errorDiv.appendChild(actions);

    const details = document.createElement('details');
    details.className = 'error-technical';
    const summary = document.createElement('summary');
    summary.textContent = 'Technical Details (Click to expand)';
    details.appendChild(summary);
    const pre = document.createElement('pre');
    const stack = (error && (error.stack || error.message)) ? String(error.stack || error.message) : '';
    pre.textContent = stack;
    details.appendChild(pre);
    errorDiv.appendChild(details);

    resultsGrid.innerHTML = '';
    resultsGrid.appendChild(errorDiv);
  }

  global.ErrorUI = { showMessage, showDetailedError };
})(window);


