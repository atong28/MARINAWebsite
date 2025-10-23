// Global state
let currentTab = 'hsqc';
let helpVisible = false;

// Collapsible sections state
const collapsibleState = {
    'visualizations': true,
    'example-molecules': true
};

// Sample data for example molecules
const exampleMolecules = {
    ruticarponineA: {
        name: 'Ruticarponine A',
        smiles: 'COC1=C(CCC(C)(C)O)C(=O)NC2=C(O)C=CC=C12',
        hsqc: [
            7.14, 112.8, 1,
            7.01, 121.7, 1,
            6.92, 114.0, 1,
            2.55, 19.3, -1,
            1.56, 42.3, -1,
            1.15, 29.2, 1,
            1.15, 29.2, 1,
            3.86, 61.5, 1
        ],
        h_nmr: [7.14, 7.01, 6.92, 2.55, 1.56, 1.15, 1.15, 10.31, 10.22, 4.24, 3.86],
        c_nmr: [162.8, 123.8, 160.8, 112.8, 121.7, 114.0, 143.8, 126.9, 117.2, 19.3, 42.3, 68.9, 29.2, 29.2, 61.5],
        mass_spec: [],
        mw: 277.131
    },
    homatropineHydrobromide: {
        name: 'Homatropine hydrobromide',
        smiles: 'CN1C2CCC1CC(OC(=O)C(O)c1ccccc1)C2',
        hsqc: [
            7.19, 128.33, 10867,
            7.29, 128.27, 4768.7,
            7.4895, 126.39, 10559,
            5.1294, 76.94, 5882.4,
            5.035, 68.19, 7546.1,
            2.9696, 59.31, 10559,
            2.1943, 40.25, 22638,
            2.03, 36.24, -15739,
            1.8285, 36.24, -13941,
            2.2352, 25.28, -10559,
            1.8792, 25.28, -12971
        ],
        h_nmr: [2.13, 2.25, 3.18, 2.03, 5.14, 7.41, 7.32],
        c_nmr: [128.58, 128.47, 66.91, 35.73, 39.2, 72.96, 138.32, 172.05, 25.37, 61.55, 126.58],
        mass_spec: [
            68.0495, 2.1604,
            70.0651, 1.9298,
            82.0651, 2.1584,
            83.073, 1.3015,
            91.0542, 3.9915,
            93.0699, 8.8415,
            107.0491, 7.3605,
            110.0964, 1.6722,
            124.1121, 62.9276,
            125.1199, 2.7327,
            140.107, 1.3497,
            142.1226, 32.9808,
            338.075, 1.7876,
            356.0856, 62.8202
        ],
        mw: 355.0783
    }
};

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
    
    // Update visualizations after adding row
    if (typeof updateVisualizations === 'function') {
        updateVisualizations();
    }
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
        
        // Update visualizations after removing row
        if (typeof updateVisualizations === 'function') {
            updateVisualizations();
        }
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
    
    // Update visualizations after clearing
    if (typeof updateVisualizations === 'function') {
        updateVisualizations();
    }
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

// Load example molecule data
function loadExampleMolecule(moleculeKey) {
    const molecule = exampleMolecules[moleculeKey];
    if (!molecule) return;
    
    // Clear all existing data
    clearTable('hsqc');
    clearTable('h_nmr');
    clearTable('c_nmr');
    clearTable('mass_spec');
    
    // Load HSQC data
    if (molecule.hsqc && molecule.hsqc.length > 0) {
        loadDataIntoTable('hsqc', molecule.hsqc, 3);
    }
    
    // Load H NMR data
    if (molecule.h_nmr && molecule.h_nmr.length > 0) {
        loadDataIntoTable('h_nmr', molecule.h_nmr, 1);
    }
    
    // Load C NMR data
    if (molecule.c_nmr && molecule.c_nmr.length > 0) {
        loadDataIntoTable('c_nmr', molecule.c_nmr, 1);
    }
    
    // Load Mass Spec data
    if (molecule.mass_spec && molecule.mass_spec.length > 0) {
        loadDataIntoTable('mass_spec', molecule.mass_spec, 2);
    }
    
    // Load Molecular Weight
    if (molecule.mw) {
        const mwInput = document.getElementById('mw-input');
        if (mwInput) {
            mwInput.value = molecule.mw;
        }
    }
    
    // Update summary and visualizations
    updateInputSummary();
    updateVisualizations();
    
    // Show success message
    showMessage(`Loaded ${molecule.name} example data`, 'success');
}

// Helper function to load data into a specific table
function loadDataIntoTable(tableId, data, columnsPerRow) {
    const table = document.getElementById(`${tableId}-table`);
    const tbody = table.querySelector('tbody');
    
    // Clear existing rows except the first one
    while (tbody.children.length > 1) {
        tbody.removeChild(tbody.lastChild);
    }
    
    // Add rows as needed
    const numRows = Math.ceil(data.length / columnsPerRow);
    for (let i = 1; i < numRows; i++) {
        addRow(tableId);
    }
    
    // Fill data
    for (let i = 0; i < data.length; i += columnsPerRow) {
        const rowIndex = Math.floor(i / columnsPerRow);
        const row = tbody.children[rowIndex];
        const inputs = row.querySelectorAll('input');
        
        for (let j = 0; j < columnsPerRow && (i + j) < data.length; j++) {
            if (inputs[j]) {
                inputs[j].value = data[i + j];
            }
        }
    }
    
    // Update visualizations after loading data
    if (typeof updateVisualizations === 'function') {
        updateVisualizations();
    }
    
    // Update visualization tab visibility
    updateVisualizationTabVisibility();
}

// Toggle collapsible sections
function toggleSection(sectionId) {
    const content = document.getElementById(`${sectionId}-content`);
    const icon = document.getElementById(`${sectionId}-icon`);
    
    if (content && icon) {
        collapsibleState[sectionId] = !collapsibleState[sectionId];
        
        if (collapsibleState[sectionId]) {
            content.classList.remove('collapsed');
            icon.classList.remove('rotated');
        } else {
            content.classList.add('collapsed');
            icon.classList.add('rotated');
        }
    }
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
    
    // Handle new optimized response format
    if (result.results && Array.isArray(result.results)) {
        const results = result.results;
        
        if (results.length === 0) {
            resultsGrid.innerHTML = '<p class="error-message">No results returned from the model.</p>';
            return;
        }

        // Store results globally for analysis page access
        currentResults = results;
        
        // Clear and display all results immediately
        resultsGrid.innerHTML = '';
        const maxResults = Math.min(10, results.length);
        
        for (let i = 0; i < maxResults; i++) {
            const resultData = results[i];
            const card = createResultCard(i + 1, resultData);
            resultsGrid.appendChild(card);
        }
        
        return;
    }
    
    // Fallback for old response format (backward compatibility)
    const indices = result.topk_indices || result.topk || [];
    const scores = result.topk_scores || [];
    
    if (!indices || indices.length === 0) {
        resultsGrid.innerHTML = '<p class="error-message">No results returned from the model.</p>';
        return;
    }

    // Create 10 blank result cards first
    resultsGrid.innerHTML = '';
    const maxResults = Math.min(10, indices.length);
    
    for (let i = 0; i < maxResults; i++) {
        const card = createBlankResultCard(i + 1);
        resultsGrid.appendChild(card);
    }
    
    // Process results and fill them in gradually
    for (let i = 0; i < maxResults; i++) {
        const idx = indices[i];
        const similarity = scores[i] || 0.0;
        
        try {
            // Fetch metadata
            const metaResponse = await fetch(`/meta/${idx}`);
            if (!metaResponse.ok) {
                updateResultCard(i, null, null, null, 'Failed to load metadata');
                continue;
            }
            
            const meta = await metaResponse.json();
            const smiles = meta.smiles;
            
            if (!smiles) {
                updateResultCard(i, null, null, null, 'No SMILES data available');
                continue;
            }
            
            // Update the card with real data
            updateResultCard(i, idx, smiles, similarity, null);
            
            // Add a small delay to show gradual loading
            await new Promise(resolve => setTimeout(resolve, 200));
            
        } catch (error) {
            console.warn(`Failed to process result ${idx}:`, error);
            updateResultCard(i, null, null, null, 'Error loading data');
        }
    }
}

// Create result card with complete data
function createResultCard(position, resultData) {
    const card = document.createElement('div');
    card.className = 'result-card';
    
    // Robust similarity handling
    let similarity = 0.0;
    if (resultData && typeof resultData.similarity === 'number' && !isNaN(resultData.similarity)) {
        similarity = resultData.similarity;
    }
    
    const smiles = (resultData && resultData.smiles) ? resultData.smiles : '';
    const svg = (resultData && resultData.svg) ? resultData.svg : '';
    const name = (resultData && resultData.name) ? resultData.name : '';
    const primaryLink = (resultData && resultData.primary_link) ? resultData.primary_link : '';
    const databaseLinks = (resultData && resultData.database_links) ? resultData.database_links : {};
    const npPathway = (resultData && resultData.np_pathway) ? resultData.np_pathway : null;
    const npSuperclass = (resultData && resultData.np_superclass) ? resultData.np_superclass : null;
    const npClass = (resultData && resultData.np_class) ? resultData.np_class : null;
    
    // Ensure similarity is a valid number and convert to percentage
    const similarityPercent = (similarity * 100).toFixed(1);
    
    // Handle both SVG and base64 image formats
    let moleculeDisplay = '<div class="no-molecule">No structure available</div>';
    if (svg) {
        if (svg.startsWith('data:image/')) {
            // Base64 image format (enhanced visualization)
            moleculeDisplay = `<img src="${svg}" alt="Molecular structure" class="molecule-image" style="max-width: 100%; max-height: 200px; border-radius: 4px;">`;
        } else if (svg.startsWith('<svg') || svg.startsWith('<')) {
            // SVG format (basic visualization)
            moleculeDisplay = svg;
        }
    }
    
    // Create database links section
    let databaseLinksHtml = '';
    if (Object.keys(databaseLinks).length > 0) {
        const linkItems = [];
        if (databaseLinks.coconut) {
            linkItems.push(`<a href="${databaseLinks.coconut}" target="_blank" class="db-link coconut">COCONUT</a>`);
        }
        if (databaseLinks.lotus) {
            linkItems.push(`<a href="${databaseLinks.lotus}" target="_blank" class="db-link lotus">LOTUS</a>`);
        }
        if (databaseLinks.npmrd) {
            linkItems.push(`<a href="${databaseLinks.npmrd}" target="_blank" class="db-link npmrd">NP-MRD</a>`);
        }
        databaseLinksHtml = `<div class="database-links">${linkItems.join(' | ')}</div>`;
    }
    
    // Create name section with primary link
    let nameHtml = '';
    if (name) {
        if (primaryLink) {
            nameHtml = `<div class="molecule-name"><a href="${primaryLink}" target="_blank" class="primary-link">${name}</a></div>`;
        } else {
            nameHtml = `<div class="molecule-name">${name}</div>`;
        }
    }
    
    // Create placeholder fields section
    let placeholderFieldsHtml = '';
    if (npPathway !== null || npSuperclass !== null || npClass !== null) {
        const fields = [];
        if (npPathway !== null) fields.push(`<div class="placeholder-field"><strong>NP Pathway:</strong> <span class="placeholder-value">${npPathway || 'N/A'}</span></div>`);
        if (npSuperclass !== null) fields.push(`<div class="placeholder-field"><strong>NP Superclass:</strong> <span class="placeholder-value">${npSuperclass || 'N/A'}</span></div>`);
        if (npClass !== null) fields.push(`<div class="placeholder-field"><strong>NP Class:</strong> <span class="placeholder-value">${npClass || 'N/A'}</span></div>`);
        placeholderFieldsHtml = `<div class="placeholder-fields">${fields.join('')}</div>`;
    }
    
    card.innerHTML = `
        <div class="card-header">
            <div class="card-title">Result #${position}</div>
            <div class="similarity-score">${similarityPercent}%</div>
        </div>
        <div class="card-body">
            <div id="molecule-${position}" class="molecule-container">
                ${moleculeDisplay}
            </div>
            ${nameHtml}
            ${databaseLinksHtml}
            <div><strong>SMILES:</strong></div>
            <div class="smiles-code">${smiles}</div>
            <div><strong>Similarity:</strong> <span>${similarityPercent}%</span></div>
            ${placeholderFieldsHtml}
            <div class="card-actions">
                <button class="btn btn-primary btn-sm" onclick="openAnalysis(${position - 1})">
                    <i class="fas fa-microscope"></i> Analyze
                </button>
            </div>
        </div>
    `;
    
    return card;
}

// Create blank result card element (for backward compatibility)
function createBlankResultCard(position) {
    const card = document.createElement('div');
    card.className = 'result-card loading';
    
    card.innerHTML = `
        <div class="card-header">
            <div class="card-title">Result #${position}</div>
            <div class="similarity-score loading">Loading...</div>
        </div>
        <div class="card-body">
            <div id="molecule-${position}" class="molecule-container loading">
                <div class="loading-spinner"></div>
            </div>
            <div><strong>SMILES:</strong></div>
            <div class="smiles-code loading">Loading...</div>
            <div><strong>Similarity:</strong> <span class="loading">Loading...</span></div>
        </div>
    `;
    
    return card;
}

// Update result card with actual data
function updateResultCard(position, idx, smiles, tanimoto, errorMessage) {
    const card = document.querySelector(`#results-grid .result-card:nth-child(${position + 1})`);
    if (!card) return;
    
    if (errorMessage) {
        // Show error state
        card.innerHTML = `
            <div class="card-header">
                <div class="card-title">Result #${position + 1}</div>
                <div class="similarity-score error">Error</div>
            </div>
            <div class="card-body">
                <div class="molecule-container error">
                    <i class="fas fa-exclamation-triangle"></i>
                </div>
                <div><strong>Error:</strong></div>
                <div class="error-message">${errorMessage}</div>
            </div>
        `;
        card.className = 'result-card error';
        return;
    }
    
    // Update with real data
    card.className = 'result-card';
    card.innerHTML = `
        <div class="card-header">
            <div class="card-title">Result #${position + 1}</div>
            <div class="similarity-score">${(tanimoto * 100).toFixed(1)}%</div>
        </div>
        <div class="card-body">
            <div id="molecule-${position + 1}" class="molecule-container">
                <div class="no-molecule">Loading enhanced visualization...</div>
            </div>
            <div><strong>SMILES:</strong></div>
            <div class="smiles-code">${smiles}</div>
            <div><strong>Similarity:</strong> ${tanimoto.toFixed(3)}</div>
        </div>
    `;
    
}


// Update backend status indicator
async function updateBackendStatus() {
    const statusElement = document.getElementById('backend-status');
    if (!statusElement) return; // Exit early if element doesn't exist
    
    const indicator = statusElement.querySelector('.status-indicator');
    const statusText = statusElement.querySelector('.status-text');
    
    if (!indicator || !statusText) return; // Exit early if sub-elements don't exist
    
    // Show loading state immediately
    indicator.className = 'fas fa-circle status-indicator loading';
    statusText.textContent = 'Checking backend status...';
    
    try {
        // Add timeout to prevent hanging on slow backend responses
        const controller = new AbortController();
        const timeoutId = setTimeout(() => controller.abort(), 5000); // 5 second timeout
        
        const response = await fetch('/health', {
            method: 'GET',
            headers: {
                'Accept': 'application/json',
                'Cache-Control': 'no-cache'
            },
            signal: controller.signal
        });
        
        clearTimeout(timeoutId);
        const data = await response.json();
        
        if (response.ok && data.model_loaded) {
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

// Global variables for analysis state
let currentAnalysisResult = null;
let originalAnalysisData = null;
let currentResults = [];

// Page navigation functions
function openAnalysis(resultIndex) {
    console.log('openAnalysis called with index:', resultIndex);
    console.log('currentResults:', currentResults);
    
    if (!currentResults || !currentResults[resultIndex]) {
        console.error('No result data available for analysis');
        showMessage('No result data available for analysis', 'error');
        return;
    }
    
    currentAnalysisResult = currentResults[resultIndex];
    originalAnalysisData = collectInputData();
    
    console.log('currentAnalysisResult:', currentAnalysisResult);
    console.log('originalAnalysisData:', originalAnalysisData);
    
    // Show analysis page
    const mainPage = document.getElementById('main-page');
    const analysisPage = document.getElementById('analysis-page');
    
    console.log('mainPage element:', mainPage);
    console.log('analysisPage element:', analysisPage);
    
    if (mainPage) {
        mainPage.style.display = 'none';
    }
    if (analysisPage) {
        analysisPage.style.display = 'block';
    }
    
    // Populate selected molecule info
    populateSelectedMolecule(currentAnalysisResult);
    
    // Copy original data to analysis page
    copyDataToAnalysisPage(originalAnalysisData);
    
    // Initially hide visualizations section
    hideAnalysisVisualizations();
    
    // Update visualizations
    if (typeof updateVisualizationsAnalysis === 'function') {
        updateVisualizationsAnalysis();
    }
    
    console.log('Analysis page should now be visible');
}

function goBackToMain() {
    // Show main page
    document.getElementById('main-page').style.display = 'block';
    document.getElementById('analysis-page').style.display = 'none';
    
    // Reset analysis state
    currentAnalysisResult = null;
    originalAnalysisData = null;
}

function populateSelectedMolecule(result) {
    const container = document.getElementById('selected-molecule-container');
    const info = document.getElementById('selected-molecule-info');
    const originalContainer = document.getElementById('original-molecule-visualization');
    
    // Display molecule structure
    let moleculeDisplay = '<div class="no-molecule">No structure available</div>';
    if (result.svg) {
        if (result.svg.startsWith('data:image/')) {
            moleculeDisplay = `<img src="${result.svg}" alt="Molecular structure" class="molecule-image" style="max-width: 100%; max-height: 200px; border-radius: 4px;">`;
        } else if (result.svg.startsWith('<svg') || result.svg.startsWith('<')) {
            moleculeDisplay = result.svg;
        }
    }
    container.innerHTML = moleculeDisplay;
    
    // Also populate the original molecule visualization
    originalContainer.innerHTML = moleculeDisplay;
    
    // Initialize zoom/pan for both images
    if (typeof panzoom !== 'undefined') {
        const topImage = container.querySelector('.molecule-image');
        const originalImage = originalContainer.querySelector('.molecule-image');
        
        if (topImage) {
            panzoom(topImage, {
                maxZoom: 5,
                minZoom: 0.5,
                initialZoom: 1
            });
        }
        
        if (originalImage) {
            panzoom(originalImage, {
                maxZoom: 5,
                minZoom: 0.5,
                initialZoom: 1
            });
        }
    }
    
    // Display molecule info
    let nameHtml = '';
    if (result.name) {
        if (result.primary_link) {
            nameHtml = `<a href="${result.primary_link}" target="_blank" class="primary-link">${result.name}</a>`;
        } else {
            nameHtml = result.name;
        }
    }
    
    let databaseLinksHtml = '';
    if (result.database_links && Object.keys(result.database_links).length > 0) {
        const linkItems = [];
        if (result.database_links.coconut) {
            linkItems.push(`<a href="${result.database_links.coconut}" target="_blank" class="db-link coconut">COCONUT</a>`);
        }
        if (result.database_links.lotus) {
            linkItems.push(`<a href="${result.database_links.lotus}" target="_blank" class="db-link lotus">LOTUS</a>`);
        }
        if (result.database_links.npmrd) {
            linkItems.push(`<a href="${result.database_links.npmrd}" target="_blank" class="db-link npmrd">NP-MRD</a>`);
        }
        databaseLinksHtml = linkItems.join(' | ');
    }
    
    info.innerHTML = `
        <div class="molecule-name-display">${nameHtml}</div>
        <div class="molecule-smiles-display">${result.smiles}</div>
        <div class="molecule-database-links">${databaseLinksHtml}</div>
    `;
}

function copyDataToAnalysisPage(data) {
    // Copy HSQC data
    if (data.hsqc) {
        loadDataIntoTableAnalysis('hsqc', data.hsqc, 3);
    }
    
    // Copy H NMR data
    if (data.h_nmr) {
        loadDataIntoTableAnalysis('h_nmr', data.h_nmr, 1);
    }
    
    // Copy C NMR data
    if (data.c_nmr) {
        loadDataIntoTableAnalysis('c_nmr', data.c_nmr, 1);
    }
    
    // Copy Mass Spec data
    if (data.mass_spec) {
        loadDataIntoTableAnalysis('mass_spec', data.mass_spec, 2);
    }
    
    // Copy Molecular Weight
    if (data.mw) {
        const mwInput = document.getElementById('analysis-mw-input');
        if (mwInput) {
            mwInput.value = data.mw;
        }
    }
}

// Analysis page table functions
function addRowAnalysis(tableId) {
    const table = document.getElementById(`analysis-${tableId}-table`);
    const tbody = table.querySelector('tbody');
    
    let newRow;
    
    switch (tableId) {
        case 'hsqc':
            newRow = createRowAnalysis([
                { type: 'number', step: '0.01', placeholder: 'e.g., 7.2' },
                { type: 'number', step: '0.01', placeholder: 'e.g., 120.5' },
                { type: 'number', step: '0.01', placeholder: 'e.g., 1.0' },
                { type: 'button', onclick: 'removeRowAnalysis(this)', icon: 'fas fa-times', title: 'Remove row' }
            ]);
            break;
        case 'h_nmr':
            newRow = createRowAnalysis([
                { type: 'number', step: '0.01', placeholder: 'e.g., 7.2' },
                { type: 'button', onclick: 'removeRowAnalysis(this)', icon: 'fas fa-times', title: 'Remove row' }
            ]);
            break;
        case 'c_nmr':
            newRow = createRowAnalysis([
                { type: 'number', step: '0.01', placeholder: 'e.g., 120.5' },
                { type: 'button', onclick: 'removeRowAnalysis(this)', icon: 'fas fa-times', title: 'Remove row' }
            ]);
            break;
        case 'mass_spec':
            newRow = createRowAnalysis([
                { type: 'number', step: '0.01', placeholder: 'e.g., 180.5' },
                { type: 'number', step: '0.01', placeholder: 'e.g., 1000' },
                { type: 'button', onclick: 'removeRowAnalysis(this)', icon: 'fas fa-times', title: 'Remove row' }
            ]);
            break;
    }
    
    tbody.appendChild(newRow);
    
    // Update visualizations after adding row
    if (typeof updateVisualizationsAnalysis === 'function') {
        updateVisualizationsAnalysis();
    }
}

function createRowAnalysis(cellConfigs) {
    const row = document.createElement('tr');
    
    cellConfigs.forEach(config => {
        const cell = document.createElement('td');
        
        if (config.type === 'number') {
            const input = document.createElement('input');
            input.type = 'number';
            if (config.step) input.step = config.step;
            if (config.placeholder) input.placeholder = config.placeholder;
            input.className = 'spreadsheet-input small';
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

function removeRowAnalysis(button) {
    const row = button.closest('tr');
    const tbody = row.closest('tbody');
    
    // Don't remove if it's the last row
    if (tbody.children.length > 1) {
        row.remove();
        
        // Update visualizations after removing row
        if (typeof updateVisualizationsAnalysis === 'function') {
            updateVisualizationsAnalysis();
        }
    }
}

function clearTableAnalysis(tableId) {
    const table = document.getElementById(`analysis-${tableId}-table`);
    const tbody = table.querySelector('tbody');
    
    // Keep only the first row
    while (tbody.children.length > 1) {
        tbody.removeChild(tbody.lastChild);
    }
    
    // Clear all inputs in the remaining row
    tbody.querySelectorAll('input').forEach(input => {
        input.value = '';
    });
    
    // Update visualizations after clearing
    if (typeof updateVisualizationsAnalysis === 'function') {
        updateVisualizationsAnalysis();
    }
}

function swapHSQCColumnsAnalysis() {
    const table = document.getElementById('analysis-hsqc-table');
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
    
    // Update visualizations
    if (typeof updateVisualizationsAnalysis === 'function') {
        updateVisualizationsAnalysis();
    }
}

function toggleSectionAnalysis(sectionId) {
    const content = document.getElementById(`analysis-${sectionId}-content`);
    const icon = document.getElementById(`analysis-${sectionId}-icon`);
    
    if (content && icon) {
        const isCollapsed = content.classList.contains('collapsed');
        
        if (isCollapsed) {
            content.classList.remove('collapsed');
            icon.classList.remove('rotated');
        } else {
            content.classList.add('collapsed');
            icon.classList.add('rotated');
        }
    }
}

// Run analysis function
async function runAnalysis() {
    if (!currentAnalysisResult) {
        showMessage('No molecule selected for analysis', 'error');
        return;
    }
    
    const data = collectInputDataAnalysis();
    const errors = validateInputData(data);
    
    if (errors.length > 0) {
        errors.forEach(error => showMessage(error, 'error'));
        return;
    }
    
    // Hide visualizations during analysis
    hideAnalysisVisualizations();
    
    // Show loading state
    const btn = document.getElementById('analysis-predict-btn');
    const originalText = btn.innerHTML;
    btn.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Running Analysis...';
    btn.disabled = true;
    
    try {
        const requestBody = {
            current_data: data,
            original_data: originalAnalysisData,
            target_smiles: currentAnalysisResult.smiles,
            target_index: currentAnalysisResult.index
        };
        
        const response = await fetch('/analyze', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify(requestBody)
        });
        
        if (!response.ok) {
            let errorMessage = 'Analysis request failed.';
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
        
        // Process analysis results
        await processAnalysisResults(result);
        
        // Show visualizations after analysis is complete
        showAnalysisVisualizations();
        
    } catch (error) {
        showMessage(error.message, 'error');
        console.error('Analysis error:', error);
        // Show visualizations even if there was an error
        showAnalysisVisualizations();
    } finally {
        btn.innerHTML = originalText;
        btn.disabled = false;
    }
}

function hideAnalysisVisualizations() {
    const visualizationsSection = document.querySelector('.analysis-visualizations-section');
    if (visualizationsSection) {
        visualizationsSection.style.display = 'none';
    }
}

function showAnalysisVisualizations() {
    const visualizationsSection = document.querySelector('.analysis-visualizations-section');
    if (visualizationsSection) {
        visualizationsSection.style.display = 'block';
    }
}

async function processAnalysisResults(result) {
    // Update similarity visualization
    await updateSimilarityVisualization(result.similarity_visualization);
    
    // Update change visualization
    await updateChangeVisualization(result.change_visualization);
    
    // Update fingerprint differences
    console.log('Processing fingerprint differences:', result.fingerprint_differences);
    if (result.fingerprint_differences) {
        updateFingerprintDifferences(result.fingerprint_differences);
    } else {
        console.log('No fingerprint differences found in result');
    }
}

async function updateSimilarityVisualization(imageData) {
    const container = document.getElementById('similarity-visualization');
    
    if (imageData) {
        const moleculeDisplay = `<img src="${imageData}" alt="Similarity highlighting" class="molecule-image">`;
        container.innerHTML = moleculeDisplay;
        
        // Initialize zoom/pan for the image
        const img = container.querySelector('.molecule-image');
        if (img && typeof panzoom !== 'undefined') {
            panzoom(img, {
                maxZoom: 5,
                minZoom: 0.5,
                initialZoom: 1
            });
        }
    } else {
        container.innerHTML = '<div class="no-molecule">Failed to generate similarity visualization</div>';
    }
}

async function updateChangeVisualization(imageData) {
    const container = document.getElementById('change-visualization');
    
    if (imageData) {
        const moleculeDisplay = `<img src="${imageData}" alt="Change highlighting" class="molecule-image">`;
        container.innerHTML = moleculeDisplay;
        
        // Initialize zoom/pan for the image
        const img = container.querySelector('.molecule-image');
        if (img && typeof panzoom !== 'undefined') {
            panzoom(img, {
                maxZoom: 5,
                minZoom: 0.5,
                initialZoom: 1
            });
        }
    } else {
        container.innerHTML = '<div class="no-molecule">Failed to generate change visualization</div>';
    }
}

function updateFingerprintDifferences(differences) {
    console.log('updateFingerprintDifferences called with:', differences);
    
    // Show the fingerprint differences section
    const fingerprintSection = document.querySelector('.fingerprint-differences-section');
    console.log('Fingerprint section found:', fingerprintSection);
    if (fingerprintSection) {
        fingerprintSection.style.display = 'block';
        console.log('Fingerprint section made visible');
    }
    
    // Update added bits list
    const addedContainer = document.getElementById('added-bits-list');
    if (differences.added && differences.added.length > 0) {
        addedContainer.innerHTML = differences.added.map(bit => `
            <div class="bit-item">
                <div class="bit-header">Bit ${bit.bit_id}</div>
                <div class="bit-details">
                    Atom: ${bit.atom_symbol} | Radius: ${bit.radius}
                </div>
                <div class="fragment-structure">
                    ${bit.substructure_image ? `<img src="${bit.substructure_image}" alt="Substructure" class="substructure-image">` : ''}
                    <span class="smiles-text">${bit.fragment_smiles}</span>
                </div>
            </div>
        `).join('');
    } else {
        addedContainer.innerHTML = '<div class="no-differences">No added bits detected</div>';
    }
    
    // Update removed bits list
    const removedContainer = document.getElementById('removed-bits-list');
    if (differences.removed && differences.removed.length > 0) {
        removedContainer.innerHTML = differences.removed.map(bit => `
            <div class="bit-item">
                <div class="bit-header">Bit ${bit.bit_id}</div>
                <div class="bit-details">
                    Atom: ${bit.atom_symbol} | Radius: ${bit.radius}
                </div>
                <div class="fragment-structure">
                    ${bit.substructure_image ? `<img src="${bit.substructure_image}" alt="Substructure" class="substructure-image">` : ''}
                    <span class="smiles-text">${bit.fragment_smiles}</span>
                </div>
            </div>
        `).join('');
    } else {
        removedContainer.innerHTML = '<div class="no-differences">No removed bits detected</div>';
    }
}

function collectInputDataAnalysis() {
    const data = {};
    
    // HSQC data
    const hsqcTable = document.getElementById('analysis-hsqc-table');
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
                const table = document.getElementById('analysis-hsqc-table');
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
    const hNmrTable = document.getElementById('analysis-h_nmr-table');
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
    const cNmrTable = document.getElementById('analysis-c_nmr-table');
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
    const massSpecTable = document.getElementById('analysis-mass_spec-table');
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
    const mwInput = document.getElementById('analysis-mw-input');
    if (mwInput && mwInput.value.trim()) {
        const mw = parseFloat(mwInput.value);
        if (!isNaN(mw)) {
            data.mw = mw;
        }
    }
    
    return data;
}

function loadDataIntoTableAnalysis(tableId, data, columnsPerRow) {
    const table = document.getElementById(`analysis-${tableId}-table`);
    const tbody = table.querySelector('tbody');
    
    // Clear existing rows except the first one
    while (tbody.children.length > 1) {
        tbody.removeChild(tbody.lastChild);
    }
    
    // Add rows as needed
    const numRows = Math.ceil(data.length / columnsPerRow);
    for (let i = 1; i < numRows; i++) {
        addRowAnalysis(tableId);
    }
    
    // Fill data
    for (let i = 0; i < data.length; i += columnsPerRow) {
        const rowIndex = Math.floor(i / columnsPerRow);
        const row = tbody.children[rowIndex];
        const inputs = row.querySelectorAll('input');
        
        for (let j = 0; j < columnsPerRow && (i + j) < data.length; j++) {
            if (inputs[j]) {
                inputs[j].value = data[i + j];
            }
        }
    }
    
    // Update visualizations after loading data
    if (typeof updateVisualizationsAnalysis === 'function') {
        updateVisualizationsAnalysis();
    }
}

// Update spectral visualizations for analysis page (removed - no longer needed)
function updateVisualizationsAnalysis() {
    // Analysis page visualizations have been removed
    return;
}

// Show/hide visualization tabs based on data availability
function showVisualizationTab(sectionId) {
    const section = document.getElementById(sectionId);
    if (section) {
        section.style.display = 'block';
    }
}

function hideVisualizationTab(sectionId) {
    const section = document.getElementById(sectionId);
    if (section) {
        section.style.display = 'none';
    }
}

// Update visualization tab visibility for main page
function updateVisualizationTabVisibility() {
    const data = collectInputData();
    
    // Check if any data is present
    const hasAnyData = (data.hsqc && data.hsqc.length > 0) ||
                      (data.h_nmr && data.h_nmr.length > 0) ||
                      (data.c_nmr && data.c_nmr.length > 0) ||
                      (data.mass_spec && data.mass_spec.length > 0);
    
    // Show/hide the entire visualizations section
    const visualizationsSection = document.querySelector('.visualizations-section');
    if (visualizationsSection) {
        visualizationsSection.style.display = hasAnyData ? 'block' : 'none';
    }
    
    // Update each visualization tab visibility
    if (data.hsqc && data.hsqc.length > 0) {
        showVisualizationTab('hsqc-section');
    } else {
        hideVisualizationTab('hsqc-section');
    }
    
    if (data.h_nmr && data.h_nmr.length > 0) {
        showVisualizationTab('h_nmr-section');
    } else {
        hideVisualizationTab('h_nmr-section');
    }
    
    if (data.c_nmr && data.c_nmr.length > 0) {
        showVisualizationTab('c_nmr-section');
    } else {
        hideVisualizationTab('c_nmr-section');
    }
    
    if (data.mass_spec && data.mass_spec.length > 0) {
        showVisualizationTab('mass_spec-section');
    } else {
        hideVisualizationTab('mass_spec-section');
    }
}


// Parse clipboard data and fill analysis table
function parseAndFillTableAnalysis(tableId, text) {
    const lines = text.trim().split('\n');
    const table = document.getElementById(`analysis-${tableId}-table`);
    const tbody = table.querySelector('tbody');
    
    // Clear existing data
    clearTableAnalysis(tableId);
    
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
        addRowAnalysis(tableId);
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
    
    // Analysis page visualizations have been removed
}

// Tab switching functionality
function switchMainTab(tabName) {
    // Hide all tab contents
    const tabContents = document.querySelectorAll('.tab-content');
    tabContents.forEach(tab => tab.classList.remove('active'));
    
    // Remove active class from all tab buttons
    const tabButtons = document.querySelectorAll('.tab-button');
    tabButtons.forEach(button => button.classList.remove('active'));
    
    // Show selected tab content
    const selectedTab = document.getElementById(`${tabName}-tab`);
    if (selectedTab) {
        selectedTab.classList.add('active');
    }
    
    // Add active class to clicked button
    const clickedButton = event.target.closest('.tab-button');
    if (clickedButton) {
        clickedButton.classList.add('active');
    }
}

// Initialize the application
document.addEventListener('DOMContentLoaded', function() {
    // Check backend status after a short delay to allow page to render first
    setTimeout(() => {
        updateBackendStatus();
    }, 100);
    
    // Update input summary on any input change, paste, or table modification
    document.addEventListener('input', function() {
        updateInputSummary();
        updateVisualizationTabVisibility();
    });
    document.addEventListener('paste', function() {
        updateInputSummary();
        updateVisualizationTabVisibility();
    });
    document.addEventListener('change', function() {
        updateInputSummary();
        updateVisualizationTabVisibility();
    });
    
    // Add copy-paste functionality for analysis page tables
    document.addEventListener('paste', function(event) {
        if (event.target.closest('.analysis-input-section')) {
            const targetInput = event.target;
            const table = targetInput.closest('table');
            if (table) {
                event.preventDefault();
                const text = (event.clipboardData || window.clipboardData).getData('text');
                const tableId = table.id.replace('analysis-', '').replace('-table', '');
                parseAndFillTableAnalysis(tableId, text);
            }
        }
    });
    
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
    
    // Check status every 10 seconds if not ready (reduced frequency to prevent multiple model setups)
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
    }, 10000);
});