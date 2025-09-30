// Global state
let currentTab = 'hsqc';
let helpVisible = false;

// Update input summary
function updateInputSummary() {
    const summaryContent = document.getElementById('input-summary-content');
    const inputs = [];
    
    // Check HSQC NMR
    const hsqcRows = document.querySelectorAll('#hsqc-table tbody tr');
    const hsqcData = Array.from(hsqcRows).filter(row => {
        const inputs = row.querySelectorAll('input');
        return inputs[0].value || inputs[1].value || inputs[2].value;
    });
    if (hsqcData.length > 0) {
        inputs.push({ type: 'HSQC NMR', count: hsqcData.length, icon: 'fas fa-atom' });
    }
    
    // Check H NMR
    const hNmrRows = document.querySelectorAll('#h_nmr-table tbody tr');
    const hNmrData = Array.from(hNmrRows).filter(row => {
        const input = row.querySelector('input');
        return input.value;
    });
    if (hNmrData.length > 0) {
        inputs.push({ type: 'H NMR', count: hNmrData.length, icon: 'fas fa-wave-square' });
    }
    
    // Check C NMR
    const cNmrRows = document.querySelectorAll('#c_nmr-table tbody tr');
    const cNmrData = Array.from(cNmrRows).filter(row => {
        const input = row.querySelector('input');
        return input.value;
    });
    if (cNmrData.length > 0) {
        inputs.push({ type: 'C NMR', count: cNmrData.length, icon: 'fas fa-wave-square' });
    }
    
    // Check Mass Spec
    const massSpecRows = document.querySelectorAll('#mass_spec-table tbody tr');
    const massSpecData = Array.from(massSpecRows).filter(row => {
        const inputs = row.querySelectorAll('input');
        return inputs[0].value || inputs[1].value;
    });
    if (massSpecData.length > 0) {
        inputs.push({ type: 'Mass Spec', count: massSpecData.length, icon: 'fas fa-chart-bar' });
    }
    
    // Check Molecular Weight
    const mwInput = document.getElementById('mw-input');
    if (mwInput && mwInput.value) {
        inputs.push({ type: 'Molecular Weight', count: 1, icon: 'fas fa-weight', value: mwInput.value + ' g/mol' });
    }
    
    if (inputs.length === 0) {
        summaryContent.innerHTML = '<p class="no-inputs">No spectral data entered yet</p>';
    } else {
        summaryContent.innerHTML = inputs.map(input => 
            `<div class="input-item">
                <i class="${input.icon}"></i>
                <span>${input.type}: ${input.count} ${input.count === 1 ? 'entry' : 'entries'}${input.value ? ` (${input.value})` : ''}</span>
            </div>`
        ).join('');
    }
}

// Render molecular structure
function renderMolecule(smiles, containerId) {
    fetch('/render-molecule', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
        },
        body: JSON.stringify({ smiles: smiles })
    })
    .then(response => response.json())
    .then(data => {
        if (data.error) {
            console.error('Error rendering molecule:', data.error);
            return;
        }
        
        const container = document.getElementById(containerId);
        container.innerHTML = `
            <div class="molecule-structure">
                <div class="molecule-svg">${data.svg}</div>
                <div class="molecule-info">
                    <div class="formula">${data.formula}</div>
                    <div class="mw">MW: ${data.molecular_weight.toFixed(2)} g/mol</div>
                </div>
            </div>
        `;
    })
    .catch(error => {
        console.error('Error:', error);
    });
}

// Help section functions
function toggleHelp() {
    const helpSection = document.getElementById('help-section');
    helpVisible = !helpVisible;
    
    if (helpVisible) {
        helpSection.style.display = 'block';
        helpSection.scrollIntoView({ behavior: 'smooth' });
    } else {
        helpSection.style.display = 'none';
    }
}

function switchHelpTab(tabName) {
    // Hide all help tab contents
    document.querySelectorAll('.help-tab-content').forEach(tab => {
        tab.classList.remove('active');
    });
    
    // Remove active class from all help tab buttons
    document.querySelectorAll('.help-tab').forEach(button => {
        button.classList.remove('active');
    });
    
    // Show selected help tab content
    document.getElementById(`${tabName}-help`).classList.add('active');
    
    // Add active class to clicked help tab button
    event.target.classList.add('active');
}

// Tab switching functionality
function switchTab(tabName) {
    // Hide all tab contents
    document.querySelectorAll('.tab-content').forEach(tab => {
        tab.classList.remove('active');
    });
    
    // Remove active class from all tab buttons
    document.querySelectorAll('.tab-button').forEach(button => {
        button.classList.remove('active');
    });
    
    // Show selected tab content
    document.getElementById(`${tabName}-tab`).classList.add('active');
    
    // Add active class to clicked tab button
    event.target.classList.add('active');
    
    currentTab = tabName;
}

// Add row to table
function addRow(tableId) {
    const table = document.getElementById(`${tableId}-table`);
    const tbody = table.querySelector('tbody');
    
    let newRow;
    
    switch (tableId) {
        case 'hsqc':
            newRow = createRow([
                { type: 'number', step: '0.01', placeholder: 'e.g., 7.2' },
                { type: 'number', step: '0.01', placeholder: 'e.g., 120.5' },
                { type: 'number', step: '0.01', placeholder: 'e.g., 1.0' },
                { type: 'button', onclick: 'removeRow(this)', icon: 'fas fa-times', title: 'Remove row' }
            ]);
            break;
        case 'h_nmr':
            newRow = createRow([
                { type: 'number', step: '0.01', placeholder: 'e.g., 7.2' },
                { type: 'button', onclick: 'removeRow(this)', icon: 'fas fa-times', title: 'Remove row' }
            ]);
            break;
        case 'c_nmr':
            newRow = createRow([
                { type: 'number', step: '0.01', placeholder: 'e.g., 120.5' },
                { type: 'button', onclick: 'removeRow(this)', icon: 'fas fa-times', title: 'Remove row' }
            ]);
            break;
        case 'mass_spec':
            newRow = createRow([
                { type: 'number', step: '0.01', placeholder: 'e.g., 180.5' },
                { type: 'number', step: '0.01', placeholder: 'e.g., 1000' },
                { type: 'button', onclick: 'removeRow(this)', icon: 'fas fa-times', title: 'Remove row' }
            ]);
            break;
    }
    
    tbody.appendChild(newRow);
}

// Create table row with specified cells
function createRow(cellConfigs) {
    const row = document.createElement('tr');
    
    cellConfigs.forEach(config => {
        const cell = document.createElement('td');
        
        if (config.type === 'number') {
            const input = document.createElement('input');
            input.type = 'number';
            if (config.step) input.step = config.step;
            if (config.placeholder) input.placeholder = config.placeholder;
            cell.appendChild(input);
        } else if (config.type === 'button') {
            const button = document.createElement('button');
            button.className = 'btn-icon';
            button.onclick = function() { eval(config.onclick); };
            if (config.title) button.title = config.title;
            
            const icon = document.createElement('i');
            icon.className = config.icon;
            button.appendChild(icon);
            
            cell.appendChild(button);
        }
        
        row.appendChild(cell);
    });
    
    return row;
}

// Remove row from table
function removeRow(button) {
    const row = button.closest('tr');
    const tbody = row.closest('tbody');
    
    // Don't remove if it's the last row
    if (tbody.children.length > 1) {
        row.remove();
    }
}

// Clear table
function clearTable(tableId) {
    const table = document.getElementById(`${tableId}-table`);
    const tbody = table.querySelector('tbody');
    
    // Keep only the first row
    while (tbody.children.length > 1) {
        tbody.removeChild(tbody.lastChild);
    }
    
    // Clear all inputs in the remaining row
    tbody.querySelectorAll('input').forEach(input => {
        input.value = '';
    });
}


// Parse clipboard data and fill table
function parseAndFillTable(tableId, text) {
    const lines = text.trim().split('\n');
    const table = document.getElementById(`${tableId}-table`);
    const tbody = table.querySelector('tbody');
    
    // Clear existing data
    clearTable(tableId);
    
    let dataRows = [];
    lines.forEach(line => {
        const values = line.split(/[\t,;]/).map(v => v.trim()).filter(v => v);
        if (values.length > 0) {
            dataRows.push(values);
        }
    });
    
    // Ensure we have at least one row
    if (dataRows.length === 0) return;
    
    // Add rows as needed
    while (tbody.children.length < dataRows.length) {
        addRow(tableId);
    }
    
    // Fill data
    dataRows.forEach((values, rowIndex) => {
        const row = tbody.children[rowIndex];
        const inputs = row.querySelectorAll('input');
        
        values.forEach((value, colIndex) => {
            if (inputs[colIndex]) {
                inputs[colIndex].value = value;
            }
        });
    });
}

// HSQC column swapping
function swapHSQCColumns() {
    const table = document.getElementById('hsqc-table');
    const hCol = table.querySelector('.hsqc-h-col');
    const cCol = table.querySelector('.hsqc-c-col');
    
    // Swap headers
    const hText = hCol.textContent;
    const cText = cCol.textContent;
    hCol.textContent = cText;
    cCol.textContent = hText;
    
    // Swap column classes
    hCol.classList.toggle('swapped');
    cCol.classList.toggle('swapped');
    
    // Swap data in all rows
    const rows = table.querySelectorAll('tbody tr');
    rows.forEach(row => {
        const cells = row.querySelectorAll('td');
        if (cells.length >= 3) {
            const hInput = cells[0].querySelector('input');
            const cInput = cells[1].querySelector('input');
            
            const hValue = hInput.value;
            const cValue = cInput.value;
            
            hInput.value = cValue;
            cInput.value = hValue;
        }
    });
    
    // Update input summary
    updateInputSummary();
}

// Show message to user
function showMessage(message, type = 'info') {
    const messageDiv = document.createElement('div');
    messageDiv.className = `${type}-message`;
    messageDiv.textContent = message;
    
    // Insert at the top of the container
    const container = document.querySelector('.container');
    container.insertBefore(messageDiv, container.firstChild);
    
    // Remove after 5 seconds
    setTimeout(() => {
        messageDiv.remove();
    }, 5000);
}

// Show detailed error message with troubleshooting
function showDetailedError(userMessage, error) {
    const resultsGrid = document.getElementById('results-grid');
    
    const errorDiv = document.createElement('div');
    errorDiv.className = 'error-details';
    errorDiv.innerHTML = `
        <div class="error-header">
            <i class="fas fa-exclamation-triangle"></i>
            <h3>Prediction Failed</h3>
        </div>
        <div class="error-message">
            <p><strong>Error:</strong> ${userMessage}</p>
        </div>
        <div class="error-troubleshooting">
            <h4>Troubleshooting Steps:</h4>
            <ul>
                <li>Check if the Docker container is running: <code>docker ps</code></li>
                <li>Verify the data folder is mounted correctly</li>
                <li>Check server logs: <code>docker-compose logs -f</code></li>
                <li>Restart the container: <code>./restart.sh</code></li>
                <li>Ensure all required files are in the data/ directory</li>
            </ul>
        </div>
        <div class="error-actions">
            <button class="btn btn-primary" onclick="retryPrediction()">
                <i class="fas fa-redo"></i> Retry Prediction
            </button>
            <button class="btn btn-secondary" onclick="checkBackendHealth()">
                <i class="fas fa-heartbeat"></i> Check Backend Health
            </button>
        </div>
        <details class="error-technical">
            <summary>Technical Details (Click to expand)</summary>
            <pre>${error.stack || error.message}</pre>
        </details>
    `;
    
    resultsGrid.innerHTML = '';
    resultsGrid.appendChild(errorDiv);
}

// Retry prediction function
function retryPrediction() {
    runPrediction();
}

// Check backend health function
async function checkBackendHealth() {
    try {
        const response = await fetch('/health');
        const data = await response.json();
        
        if (response.ok && data.model_loaded) {
            showMessage('Backend is healthy and model is ready!', 'success');
        } else if (data.status === 'initializing') {
            showMessage('Backend is running but model is still initializing...', 'warning');
        } else if (data.error) {
            showMessage(`Model initialization failed: ${data.error}`, 'error');
    } else {
            showMessage('Backend health check failed', 'error');
        }
    } catch (error) {
        showMessage('Cannot connect to backend health endpoint', 'error');
    }
}

// Collect data from all inputs
function collectInputData() {
    const data = {};
    
    // HSQC data
    const hsqcTable = document.getElementById('hsqc-table');
    const hsqcRows = hsqcTable.querySelectorAll('tbody tr');
    const hsqcData = [];
    
    hsqcRows.forEach((row, index) => {
        const inputs = row.querySelectorAll('input');
        
        if (inputs.length >= 3) {
            const hShift = parseFloat(inputs[0].value);
            const cShift = parseFloat(inputs[1].value);
            const intensity = parseFloat(inputs[2].value);
            
            if (!isNaN(hShift) && !isNaN(cShift) && !isNaN(intensity) && inputs[0].value.trim() !== '' && inputs[1].value.trim() !== '' && inputs[2].value.trim() !== '') {
                // Check if columns are swapped by looking at header classes
                const table = document.getElementById('hsqc-table');
                const hCol = table.querySelector('.hsqc-h-col');
                const isSwapped = hCol.classList.contains('swapped');
                
                if (isSwapped) {
                    // Columns are swapped, so first input is C, second is H
                    hsqcData.push([cShift, hShift, intensity]);
                } else {
                    // Normal order: H, C, intensity
                    hsqcData.push([hShift, cShift, intensity]);
                }
            }
        }
    });
    
    if (hsqcData.length > 0) {
        data.hsqc = hsqcData.flat();
    }
    
    // H NMR data
    const hNmrTable = document.getElementById('h_nmr-table');
    const hNmrRows = hNmrTable.querySelectorAll('tbody tr');
    const hNmrData = [];
    
    hNmrRows.forEach(row => {
        const input = row.querySelector('input');
        if (input && input.value.trim()) {
            const value = parseFloat(input.value);
            if (!isNaN(value)) {
                hNmrData.push(value);
            }
        }
    });
    
    if (hNmrData.length > 0) {
        data.h_nmr = hNmrData;
    }
    
    // C NMR data
    const cNmrTable = document.getElementById('c_nmr-table');
    const cNmrRows = cNmrTable.querySelectorAll('tbody tr');
    const cNmrData = [];
    
    cNmrRows.forEach(row => {
        const input = row.querySelector('input');
        if (input && input.value.trim()) {
            const value = parseFloat(input.value);
            if (!isNaN(value)) {
                cNmrData.push(value);
            }
        }
    });
    
    if (cNmrData.length > 0) {
        data.c_nmr = cNmrData;
    }
    
    // Mass spec data
    const massSpecTable = document.getElementById('mass_spec-table');
    const massSpecRows = massSpecTable.querySelectorAll('tbody tr');
    const massSpecData = [];
    
    massSpecRows.forEach(row => {
        const inputs = row.querySelectorAll('input');
        if (inputs.length >= 2) {
            const mz = parseFloat(inputs[0].value);
            const intensity = parseFloat(inputs[1].value);
            
            if (!isNaN(mz) && !isNaN(intensity)) {
                massSpecData.push([mz, intensity]);
            }
        }
    });
    
    if (massSpecData.length > 0) {
        data.mass_spec = massSpecData.flat();
    }
    
    // Molecular weight
    const mwInput = document.getElementById('mw-input');
    if (mwInput && mwInput.value.trim()) {
        const mw = parseFloat(mwInput.value);
        if (!isNaN(mw)) {
            data.mw = mw;
        }
    }
    
    return data;
}

// Validate input data
function validateInputData(data) {
    const errors = [];
    
    if (Object.keys(data).length === 0) {
        errors.push('Please enter at least one type of spectral data.');
    }
    
    // Check for required fields in each data type
    if (data.hsqc && data.hsqc.length % 3 !== 0) {
        errors.push('HSQC data must have complete triplets (C-shift, H-shift, intensity).');
    }
    
    if (data.mass_spec && data.mass_spec.length % 2 !== 0) {
        errors.push('Mass spec data must have complete pairs (m/z, intensity).');
    }
    
    if (data.mw && (data.mw <= 0 || data.mw > 10000)) {
        errors.push('Molecular weight should be between 0 and 10000 g/mol.');
    }
    
    return errors;
}

// Run prediction
async function runPrediction() {
    const data = collectInputData();
    const errors = validateInputData(data);
    
    if (errors.length > 0) {
        errors.forEach(error => showMessage(error, 'error'));
        return;
    }
    
    const k = parseInt(document.getElementById('top-k').value);
    
    // Show loading state
    const resultsSection = document.getElementById('results-section');
    const loading = document.getElementById('loading');
    const resultsGrid = document.getElementById('results-grid');
    
    resultsSection.style.display = 'block';
    loading.style.display = 'flex';
    resultsGrid.innerHTML = '';
    
    try {
        // First check if backend is available
        const healthResponse = await fetch('/health', { 
            method: 'GET',
            timeout: 5000 
        });
        
        if (!healthResponse.ok) {
            throw new Error('Backend service is not available. Please check if the server is running.');
        }
        
        const healthData = await healthResponse.json();
        
        // Check if model is ready
        if (!healthData.model_loaded) {
            if (healthData.error) {
                throw new Error(`Model initialization failed: ${healthData.error}`);
            } else {
                throw new Error('Model is still initializing. Please wait a moment and try again.');
            }
        }
        
        const requestBody = {
            raw: data,
            k: k
        };
        
        const response = await fetch('/predict', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify(requestBody)
        });
        
        if (!response.ok) {
            let errorMessage = 'Prediction request failed.';
            try {
                const errorData = await response.json();
                if (errorData.error) {
                    errorMessage = errorData.error;
                }
            } catch (e) {
                // If we can't parse the error response, use default message
                if (response.status === 404) {
                    errorMessage = 'Prediction endpoint not found. Please check server configuration.';
                } else if (response.status === 500) {
                    errorMessage = 'Internal server error. The model may not be loaded correctly.';
                } else if (response.status === 503) {
                    errorMessage = 'Service temporarily unavailable. Please try again in a moment.';
                } else {
                    errorMessage = `Server error (${response.status}). Please contact support.`;
                }
            }
            throw new Error(errorMessage);
        }
        
        const result = await response.json();
        
        if (result.error) {
            throw new Error(result.error);
        }
        
        // Hide loading, show results
        loading.style.display = 'none';
        displayResults(result);
        
    } catch (error) {
        loading.style.display = 'none';
        
        // Handle different types of errors with user-friendly messages
        let userMessage = error.message;
        
        if (error.name === 'TypeError' && error.message.includes('fetch')) {
            userMessage = 'Cannot connect to the prediction server. Please ensure the backend is running and accessible.';
        } else if (error.message.includes('timeout')) {
            userMessage = 'Request timed out. The server may be overloaded. Please try again.';
        } else if (error.message.includes('NetworkError')) {
            userMessage = 'Network error. Please check your internet connection and try again.';
        }
        
        // Show detailed error message
        showDetailedError(userMessage, error);
        console.error('Prediction error:', error);
    }
}

// Display prediction results
async function displayResults(result) {
    const resultsGrid = document.getElementById('results-grid');
    resultsGrid.innerHTML = '';
    
    const indices = result.topk_indices || result.topk || [];
    const scores = result.topk_scores || [];
    const predFp = result.pred_fp || result.prediction || null;
    
    if (!indices || indices.length === 0) {
        resultsGrid.innerHTML = '<p class="error-message">No results returned from the model.</p>';
      return;
    }

    // Process results
    for (let i = 0; i < Math.min(indices.length, 20); i++) {
      const idx = indices[i];
      const similarity = scores[i] || 0.0; // Use backend-computed similarity score
        
        try {
            // Fetch metadata
            const metaResponse = await fetch(`/meta/${idx}`);
            if (!metaResponse.ok) continue;
            
            const meta = await metaResponse.json();
      const smiles = meta.smiles;
            
            if (!smiles) continue;
            
            // Create result card
            const card = createResultCard(idx, smiles, similarity);
            resultsGrid.appendChild(card);
            
        } catch (error) {
            console.warn(`Failed to process result ${idx}:`, error);
        }
    }
    
    if (resultsGrid.children.length === 0) {
        resultsGrid.innerHTML = '<p class="error-message">No valid results could be retrieved.</p>';
    }
}

// Create result card element
function createResultCard(idx, smiles, tanimoto) {
    const card = document.createElement('div');
    card.className = 'result-card';
    
    card.innerHTML = `
        <div class="card-header">
            <div class="card-title">Result #${idx}</div>
            <div class="similarity-score">${(tanimoto * 100).toFixed(1)}%</div>
        </div>
        <div class="card-body">
            <div id="molecule-${idx}" class="molecule-container"></div>
            <div><strong>SMILES:</strong></div>
            <div class="smiles-code">${smiles}</div>
            <div><strong>Similarity:</strong> ${tanimoto.toFixed(3)}</div>
        </div>
    `;
    
    // Render molecular structure
    renderMolecule(smiles, `molecule-${idx}`);
    
    return card;
}

// Update backend status indicator
async function updateBackendStatus() {
    const statusElement = document.getElementById('backend-status');
    const indicator = statusElement.querySelector('.status-indicator');
    const statusText = statusElement.querySelector('.status-text');
    
    try {
        const response = await fetch('/health');
        const data = await response.json();
        
        if (response.ok && data.model_loaded) {
            indicator.className = 'fas fa-circle status-indicator ready';
            statusText.textContent = 'Backend ready - Model loaded';
        } else if (data.status === 'loading') {
            indicator.className = 'fas fa-circle status-indicator loading';
            statusText.textContent = 'Model loading in background...';
        } else if (data.status === 'initializing') {
            indicator.className = 'fas fa-circle status-indicator initializing';
            statusText.textContent = 'Model initializing...';
        } else if (data.error) {
            indicator.className = 'fas fa-circle status-indicator error';
            statusText.textContent = 'Model initialization failed';
        } else {
            indicator.className = 'fas fa-circle status-indicator error';
            statusText.textContent = 'Backend unavailable';
        }
    } catch (error) {
        indicator.className = 'fas fa-circle status-indicator error';
        statusText.textContent = 'Cannot connect to backend';
    }
}


// SMILES search function
async function runSmilesSearch() {
    const smilesInput = document.getElementById('smiles-input');
    const smiles = smilesInput.value.trim();
    
    if (!smiles) {
        showMessage('Please enter a SMILES string', 'error');
        return;
    }
    
    const k = parseInt(document.getElementById('top-k').value);
    
    // Show loading state
    const resultsSection = document.getElementById('results-section');
    const loading = document.getElementById('loading');
    const resultsGrid = document.getElementById('results-grid');
    
    resultsSection.style.display = 'block';
    loading.style.display = 'flex';
    resultsGrid.innerHTML = '';
    
    try {
        const response = await fetch('/smiles-search', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({
                smiles: smiles,
                k: k
            })
        });
        
        if (!response.ok) {
            let errorMessage = 'SMILES search failed.';
            try {
                const errorData = await response.json();
                if (errorData.error) {
                    errorMessage = errorData.error;
                }
            } catch (e) {
                errorMessage = `Server error (${response.status}). Please try again.`;
            }
            throw new Error(errorMessage);
        }
        
        const result = await response.json();
        
        if (result.error) {
            throw new Error(result.error);
        }
        
        // Hide loading, show results
        loading.style.display = 'none';
        displayResults(result);
        
    } catch (error) {
        loading.style.display = 'none';
        showMessage(error.message, 'error');
    }
}

// Initialize the application
document.addEventListener('DOMContentLoaded', function() {
    // Check backend status on load
    updateBackendStatus();
    
    // Update input summary on any input change, paste, or table modification
    document.addEventListener('input', updateInputSummary);
    document.addEventListener('paste', updateInputSummary);
    document.addEventListener('change', updateInputSummary);
    
    // Also update when rows are added/removed
    const observer = new MutationObserver(function(mutations) {
        mutations.forEach(function(mutation) {
            if (mutation.type === 'childList' && 
                (mutation.target.tagName === 'TBODY' || mutation.target.classList.contains('data-table'))) {
                updateInputSummary();
            }
        });
    });
    
    // Observe all table bodies for changes
    document.querySelectorAll('tbody').forEach(tbody => {
        observer.observe(tbody, { childList: true, subtree: true });
    });
    
    // Initial input summary update
    updateInputSummary();
    
    // Check status every 5 seconds if not ready
    const statusCheckInterval = setInterval(async () => {
        try {
            const response = await fetch('/health');
            const data = await response.json();
            
            if (response.ok && data.model_loaded) {
                clearInterval(statusCheckInterval);
                updateBackendStatus();
            } else {
                updateBackendStatus();
            }
        } catch (error) {
            updateBackendStatus();
        }
    }, 5000);
});