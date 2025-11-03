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
    
    summaryContent.textContent = '';
    if (inputs.length === 0) {
        const p = document.createElement('p');
        p.className = 'no-inputs';
        p.textContent = 'No spectral data entered yet';
        summaryContent.appendChild(p);
    } else {
        inputs.forEach(input => {
            const item = document.createElement('div');
            item.className = 'input-item';
            const iconEl = document.createElement('i');
            iconEl.className = input.icon;
            const span = document.createElement('span');
            const countText = `${input.count} ${input.count === 1 ? 'entry' : 'entries'}`;
            span.textContent = `${input.type}: ${countText}${input.value ? ` (${input.value})` : ''}`;
            item.appendChild(iconEl);
            item.appendChild(span);
            summaryContent.appendChild(item);
        });
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
            if (config.title) button.title = config.title;
            const icon = document.createElement('i');
            icon.className = config.icon;
            button.appendChild(icon);
            // Bind known handlers without eval
            if (config.onclick === 'removeRow(this)') {
                button.addEventListener('click', function() { removeRow(button); });
            } else if (typeof config.onclick === 'function') {
                button.addEventListener('click', () => config.onclick(button));
            }
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
    if (window.Toggles && typeof window.Toggles.toggleSection === 'function') {
        return window.Toggles.toggleSection(sectionId);
    }
    const content = document.getElementById(`${sectionId}-content`);
    const icon = document.getElementById(`${sectionId}-icon`);
    if (content && icon) {
        if (content.classList.contains('collapsed')) {
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
    if (window.ErrorUI && typeof window.ErrorUI.showMessage === 'function') {
        return window.ErrorUI.showMessage(message, type);
    }
    const messageDiv = document.createElement('div');
    messageDiv.className = `${type}-message`;
    messageDiv.textContent = message;
    const container = document.querySelector('.container');
    container.insertBefore(messageDiv, container.firstChild);
    setTimeout(() => {
        messageDiv.remove();
    }, 5000);
}

// Show detailed error message with troubleshooting
function showDetailedError(userMessage, error) {
    if (window.ErrorUI && typeof window.ErrorUI.showDetailedError === 'function') {
        return window.ErrorUI.showDetailedError(userMessage, error);
    }
    const resultsGrid = document.getElementById('results-grid');
    const errorDiv = document.createElement('div');
    errorDiv.className = 'error-details';
    errorDiv.textContent = userMessage || 'Unexpected error';
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
    if (window.Prediction && typeof window.Prediction.run === 'function') {
        return await window.Prediction.run();
    }
}

// Display prediction results
async function displayResults(result) {
    if (window.Results && typeof window.Results.displayResults === 'function') {
        return await window.Results.displayResults(result);
    }
    // Fallback minimal rendering
    const resultsGrid = document.getElementById('results-grid');
    resultsGrid.innerHTML = '<p class="error-message">Renderer unavailable.</p>';
}

// Create result card with complete data
function createResultCard(position, resultData) {
    if (window.ResultRenderer && typeof window.ResultRenderer.createCard === 'function') {
        return window.ResultRenderer.createCard(position, resultData);
    }
    // Fallback: simple card
    const card = document.createElement('div');
    card.className = 'result-card';
    card.textContent = `Result #${position}`;
    return card;
}

function createBlankResultCard(position) {
    if (window.Results && typeof window.Results.createBlankResultCard === 'function') {
        return window.Results.createBlankResultCard(position);
    }
    const card = document.createElement('div');
    card.className = 'result-card';
    card.textContent = `Result #${position}`;
    return card;
}

function updateResultCard(position, idx, smiles, tanimoto, errorMessage) {
    if (window.Results && typeof window.Results.updateResultCard === 'function') {
        return window.Results.updateResultCard(position, idx, smiles, tanimoto, errorMessage);
    }
}


// Update backend status indicator
async function updateBackendStatus() {
    if (window.Health && typeof window.Health.updateBackendStatus === 'function') {
        return await window.Health.updateBackendStatus();
    }
}


// SMILES search function
async function runSmilesSearch() {
    if (window.SmilesSearch && typeof window.SmilesSearch.runSmilesSearch === 'function') {
        return await window.SmilesSearch.runSmilesSearch();
    }
}

// Global variables for analysis state
let currentAnalysisResult = null;
let originalAnalysisData = null;
let currentResults = [];
let currentPredictedFp = null;
let currentRetrievedFp = null;

// Page navigation functions
function openAnalysis(resultIndex) {
    let resultsSource = (window.State && State.get) ? State.get('currentResults', currentResults) : currentResults;
    
    if (!resultsSource || !resultsSource[resultIndex]) {
        showMessage('No result data available for analysis', 'error');
        return;
    }
    
    // Reset per-analysis caches/state
    if (window.State && State.set) {
        try { State.set('selectedBits', new Set()); } catch (e) {}
    } else {
        try { window.selectedBits = new Set(); } catch (e) {}
    }
    currentAnalysisResult = resultsSource[resultIndex];
    if (window.State && State.set) State.set('currentAnalysisResult', currentAnalysisResult);
    
    // Check if this is from SMILES search (check if analysisSource flag exists)
    const isFromSmilesSearch = (window.State && State.get && State.get('analysisSource')) === 'smiles-search';
    
    originalAnalysisData = collectInputData();
    
    // Show analysis page
    const mainPage = document.getElementById('main-page');
    const analysisPage = document.getElementById('analysis-page');
    
    if (mainPage) {
        mainPage.style.display = 'none';
    }
    if (analysisPage) {
        analysisPage.style.display = 'block';
    }
    
    // Hide spectral search and secondary search sections if from SMILES search
    if (isFromSmilesSearch) {
        const analysisInputSection = document.querySelector('.analysis-input-section');
        const secondaryRetrievalSection = document.getElementById('secondary-retrieval-section');
        if (analysisInputSection) analysisInputSection.style.display = 'none';
        if (secondaryRetrievalSection) secondaryRetrievalSection.style.display = 'none';
    } else {
        // Show sections for prediction-based analysis
        const analysisInputSection = document.querySelector('.analysis-input-section');
        if (analysisInputSection) analysisInputSection.style.display = 'block';
    }
    
    // Populate selected molecule info
    populateSelectedMolecule(currentAnalysisResult);
    
    // Copy original data to analysis page (only if not from SMILES search)
    if (!isFromSmilesSearch) {
        copyDataToAnalysisPage(originalAnalysisData);
    }
    
    // Initially hide visualizations section
    hideAnalysisVisualizations();
    
    // Update visualizations
    if (typeof updateVisualizationsAnalysis === 'function') {
        updateVisualizationsAnalysis();
    }
    
    // Instantly compute and show predicted/retrieved bit indices
    try {
        // Clear stale retrievedFpIndices from State to prevent using wrong molecule's data
        if (window.State && State.remove) {
            State.remove('retrievedFpIndices');
        }
        const pred = (window.currentPredictedFp) || (window.State && State.get && State.get('currentPredictedFp'));
        let predictedIdx = [];
        if (pred && Array.isArray(pred)) {
            for (let i = 0; i < pred.length; i++) if (pred[i] > 0.5) predictedIdx.push(i);
        }
        const retrievedIdx = (currentAnalysisResult && currentAnalysisResult.retrieved_molecule_fp_indices) || [];
        try { console.debug('[openAnalysis] calling updateFingerprintIndices', { predictedCount: predictedIdx.length, retrievedCount: retrievedIdx.length }); } catch (e) {}
        updateFingerprintIndices(predictedIdx, retrievedIdx);
        
        // Restore analysis state if same molecule was analyzed before
        try {
            if (window.State && State.restoreAnalysis) {
                const restored = State.restoreAnalysis();
                if (restored && restored.currentAnalysisResult) {
                    const currentSmiles = (currentAnalysisResult.smiles || currentAnalysisResult.target_smiles);
                    const restoredSmiles = (restored.currentAnalysisResult.smiles || restored.currentAnalysisResult.target_smiles);
                    if (currentSmiles && restoredSmiles && currentSmiles === restoredSmiles) {
                        // Same molecule - restore selections and unavailable flags
                        if (restored.selectedBits && Array.isArray(restored.selectedBits)) {
                            const sel = new Set(restored.selectedBits);
                            if (window.State && State.set) State.set('selectedBits', sel);
                            // Mark badges as selected
                            document.querySelectorAll('.index-badge').forEach(el => {
                                const bit = parseInt(el.textContent, 10);
                                if (!isNaN(bit) && sel.has(bit)) {
                                    el.classList.add('badge-selected');
                                }
                            });
                        }
                        // Mark unavailable bits (yellow) based on bit_environments
                        const bitEnvs = restored.currentAnalysisResult && restored.currentAnalysisResult.bit_environments;
                        if (bitEnvs && typeof bitEnvs === 'object') {
                            const allBits = new Set([...(restored.predictedFpIndices || []), ...(restored.retrievedFpIndices || [])]);
                            allBits.forEach(bit => {
                                if (!bitEnvs[bit]) {
                                    document.querySelectorAll('.index-badge').forEach(el => {
                                        if (parseInt(el.textContent, 10) === bit) {
                                            el.classList.add('badge-unavailable');
                                            el.title = 'No substructure mapping available for this bit';
                                        }
                                    });
                                }
                            });
                        }
                        // Rebuild overlay if there are selected bits
                        if (typeof rebuildSelectedOverlay === 'function') {
                            setTimeout(() => rebuildSelectedOverlay(), 100);
                        }
                    }
                }
            }
        } catch (e) {
            console.warn('Failed to restore analysis state:', e);
        }
    } catch (e) {}

    // Kick off secondary retrieval asynchronously with spinner
    try {
        const secondarySection = document.getElementById('secondary-retrieval-section');
        const secondaryLoading = document.getElementById('secondary-retrieval-loading');
        const secondaryGrid = document.getElementById('secondary-results-grid');
        if (secondarySection) secondarySection.style.display = 'block';
        if (secondaryLoading) secondaryLoading.style.display = 'flex';
        if (secondaryGrid) secondaryGrid.innerHTML = '';
        const predFp = (window.currentPredictedFp) || (window.State && State.get && State.get('currentPredictedFp')) || (currentAnalysisResult && currentAnalysisResult.predicted_fp) || null;
        const retrievedFp = (currentAnalysisResult && currentAnalysisResult.retrieved_molecule_fp) || null;
        if (predFp && retrievedFp) {
            ApiClient.postSecondaryRetrieval({ predicted_fp: predFp, retrieved_fp: retrievedFp, k: 10 })
                .then(resp => {
                    if (secondaryLoading) secondaryLoading.style.display = 'none';
                    if (!resp || !resp.results) return;
                    if (secondaryGrid) {
                        secondaryGrid.innerHTML = '';
                        const arr = resp.results;
                        if (!arr.length) {
                            secondaryGrid.innerHTML = '<p class="no-results-message">No results found for remaining substructure.</p>';
                        } else {
                            for (let i = 0; i < arr.length; i++) {
                                const card = (window.ResultRenderer && window.ResultRenderer.createCard)
                                  ? window.ResultRenderer.createCard(i + 1, arr[i], { showAnalyze: false })
                                  : createResultCard(i + 1, arr[i]);
                                secondaryGrid.appendChild(card);
                            }
                        }
                    }
                })
                .catch(() => { if (secondaryLoading) secondaryLoading.style.display = 'none'; });
        } else {
            if (secondaryLoading) secondaryLoading.style.display = 'none';
        }
    } catch (e) {}

    // DO NOT auto-run analysis - let users manually trigger it when ready
    // Auto-running would overwrite currentAnalysisResult and break highlighting
}

function goBackToMain() {
    // Show main page
    document.getElementById('main-page').style.display = 'block';
    document.getElementById('analysis-page').style.display = 'none';
    
    // Reset analysis state
    currentAnalysisResult = null;
    originalAnalysisData = null;
    // Clear overlay and indices UI to avoid stale DOM on next analysis
    try {
        const overlay = document.getElementById('selected-molecule-overlay-svg');
        if (overlay && overlay.parentElement) overlay.parentElement.removeChild(overlay);
    } catch (e) {}
    try {
        const pred = document.getElementById('predicted-fp-indices');
        const retr = document.getElementById('retrieved-fp-indices');
        if (pred) pred.textContent = '';
        if (retr) retr.textContent = '';
        const indicesSection = document.getElementById('fingerprint-indices-section');
        if (indicesSection) {
            // Move indices section back to original location (after analysis-visualizations-section)
            // if it was moved to selected-molecule-info
            const analysisLayout = document.querySelector('.analysis-layout');
            const visualizationsSection = document.querySelector('.analysis-visualizations-section');
            if (analysisLayout && visualizationsSection && indicesSection.parentElement !== analysisLayout) {
                // Remove from current location (e.g., selected-molecule-info)
                if (indicesSection.parentElement) {
                    indicesSection.parentElement.removeChild(indicesSection);
                }
                // Insert after analysis-visualizations-section in analysis-layout
                if (visualizationsSection.nextSibling) {
                    analysisLayout.insertBefore(indicesSection, visualizationsSection.nextSibling);
                } else {
                    analysisLayout.appendChild(indicesSection);
                }
            }
            indicesSection.style.display = 'none';
        }
    } catch (e) {}
}

function populateSelectedMolecule(result) {
    const container = document.getElementById('selected-molecule-container');
    const info = document.getElementById('selected-molecule-info');
    const originalContainer = document.getElementById('original-molecule-visualization');
    
    // Display molecule structure (plain, non-highlighted preferred)
    let moleculeDisplay = '<div class="no-molecule">No structure available</div>';
    const pick = result.plain_svg || result.svg;
    if (pick) {
        if (pick.startsWith && pick.startsWith('data:image/')) {
            moleculeDisplay = `<img src="${pick}" alt="Molecular structure" class="molecule-image" style="max-width: 100%; max-height: 200px; border-radius: 4px;">`;
        } else if (pick.startsWith && (pick.startsWith('<svg') || pick.startsWith('<'))) {
            moleculeDisplay = pick;
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
            nameHtml = { link: result.primary_link, text: result.name };
        } else {
            nameHtml = { text: result.name };
        }
    }
    
    // Preserve fingerprint-indices-section if it exists as a child (from previous analysis)
    const indicesSection = document.getElementById('fingerprint-indices-section');
    let preservedIndicesSection = null;
    if (indicesSection && indicesSection.parentElement === info) {
        // Remove it from parent before clearing, so it's not destroyed
        preservedIndicesSection = indicesSection;
        info.removeChild(indicesSection);
    }
    
    // Clear and rebuild info (now safe to clear since indices section is removed)
    info.textContent = '';
    const nameDiv = document.createElement('div');
    nameDiv.className = 'molecule-name-display';
    if (nameHtml) {
        if (nameHtml.link) {
            const a = document.createElement('a');
            a.href = nameHtml.link;
            a.target = '_blank';
            a.className = 'primary-link';
            a.textContent = nameHtml.text;
            nameDiv.appendChild(a);
        } else {
            nameDiv.textContent = nameHtml.text;
        }
    }
    
    info.appendChild(nameDiv);
    const smilesDiv = document.createElement('div');
    smilesDiv.className = 'molecule-smiles-display';
    smilesDiv.textContent = result.smiles || '';
    info.appendChild(smilesDiv);
    const linksDiv = document.createElement('div');
    linksDiv.className = 'molecule-database-links';
    if (result.database_links) {
        if (result.database_links.coconut) {
            const a = document.createElement('a');
            a.href = result.database_links.coconut;
            a.target = '_blank';
            a.className = 'db-link coconut';
            a.textContent = 'COCONUT';
            linksDiv.appendChild(a);
        }
        if (result.database_links.lotus) {
            const a = document.createElement('a');
            a.href = result.database_links.lotus;
            a.target = '_blank';
            a.className = 'db-link lotus';
            a.textContent = 'LOTUS';
            linksDiv.appendChild(a);
        }
        if (result.database_links.npmrd) {
            const a = document.createElement('a');
            a.href = result.database_links.npmrd;
            a.target = '_blank';
            a.className = 'db-link npmrd';
            a.textContent = 'NP-MRD';
            linksDiv.appendChild(a);
        }
    }
    info.appendChild(linksDiv);
    
    // Re-append preserved indices section if it was saved
    if (preservedIndicesSection) {
        info.appendChild(preservedIndicesSection);
    }
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
            if (config.title) button.title = config.title;
            const icon = document.createElement('i');
            icon.className = config.icon;
            button.appendChild(icon);
            if (config.onclick === 'removeRowAnalysis(this)') {
                button.addEventListener('click', function() { removeRowAnalysis(button); });
            } else if (typeof config.onclick === 'function') {
                button.addEventListener('click', () => config.onclick(button));
            }
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
    if (window.Toggles && typeof window.Toggles.toggleSectionAnalysis === 'function') {
        return window.Toggles.toggleSectionAnalysis(sectionId);
    }
    const content = document.getElementById(`analysis-${sectionId}-content`);
    const icon = document.getElementById(`analysis-${sectionId}-icon`);
    if (content && icon) {
        if (content.classList.contains('collapsed')) {
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
        
        const result = await ApiClient.postAnalyze(requestBody);
        
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
    // Delegate core processing to feature module
    if (window.Analysis && typeof window.Analysis.processAnalysisResults === 'function') {
        await window.Analysis.processAnalysisResults(result);
    }
    
    // Display backend-provided secondary retrieval or skip with note
    const secondarySection = document.getElementById('secondary-retrieval-section');
    const secondaryGrid = document.getElementById('secondary-results-grid');
    const secondaryNote = document.getElementById('secondary-note');
    
    if (result.secondary_skipped) {
        if (secondarySection) secondarySection.style.display = 'none';
        if (secondaryNote) {
            secondaryNote.style.display = 'block';
            secondaryNote.textContent = result.secondary_message || 'Predicted and selected fingerprints are identical.';
        }
    } else if (result.secondary_results && Array.isArray(result.secondary_results)) {
        if (secondaryNote) secondaryNote.style.display = 'none';
        if (secondarySection) secondarySection.style.display = 'block';
        if (secondaryGrid) {
            secondaryGrid.innerHTML = '';
            if (result.secondary_results.length === 0) {
                secondaryGrid.innerHTML = '<p class="no-results-message">No results found for remaining substructure.</p>';
    } else {
                const secondaryResults = result.secondary_results;
                for (let i = 0; i < secondaryResults.length; i++) {
                    const resultData = secondaryResults[i];
                    const card = (window.ResultRenderer && window.ResultRenderer.createCard)
                      ? window.ResultRenderer.createCard(i + 1, resultData, { showAnalyze: false })
                      : createResultCard(i + 1, resultData);
                    secondaryGrid.appendChild(card);
                }
            }
        }
    } else {
        // Fallback: hide section if nothing provided
        if (secondarySection) secondarySection.style.display = 'none';
    }
}

function updateFingerprintIndices(predictedIndices, retrievedIndices) {
    try { 
        console.debug('[updateFingerprintIndices] start', { 
            predictedCount: (predictedIndices||[]).length, 
            retrievedCount: (retrievedIndices||[]).length,
            predictedIndices: (predictedIndices||[]).slice(0, 10),
            retrievedIndices: (retrievedIndices||[]).slice(0, 10)
        }); 
    } catch (e) {}
    const indicesSection = document.getElementById('fingerprint-indices-section');
    if (!indicesSection) {
        try { console.error('[updateFingerprintIndices] fingerprint-indices-section not found in DOM'); } catch (e) {}
        return; // Cannot proceed without the section
    }
    
    indicesSection.style.display = 'block';
    // Move indices section under selected molecule info to keep always visible
    const infoHost = document.getElementById('selected-molecule-info');
    if (infoHost && indicesSection.parentElement !== infoHost) {
        infoHost.appendChild(indicesSection);
    }
    
    // Verify containers exist before building lists
    const predictedContainer = document.getElementById('predicted-fp-indices');
    const retrievedContainer = document.getElementById('retrieved-fp-indices');
    if (!predictedContainer || !retrievedContainer) {
        try { 
            console.error('[updateFingerprintIndices] containers not found', { 
                predicted: !!predictedContainer, 
                retrieved: !!retrievedContainer,
                indicesSectionExists: !!indicesSection,
                indicesSectionParent: indicesSection.parentElement ? indicesSection.parentElement.id : 'none'
            }); 
        } catch (e) {}
        return; // Cannot proceed without containers
    }
    
    // Helper: rebuild combined overlay for all selected bits
    const rebuildSelectedOverlay = () => {
        try { console.debug('[overlay] rebuildSelectedOverlay invoked'); } catch (e) {}
        const container = document.getElementById('selected-molecule-container');
        if (!container) { try { console.warn('[overlay] container not found'); } catch (e) {} return; }
        // Ensure positioning so overlay can sit absolutely
        try { if (getComputedStyle(container).position === 'static') { container.style.position = 'relative'; } } catch (e) {}
        let svg = document.getElementById('selected-molecule-overlay-svg');
        if (!svg) {
            svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
            svg.setAttribute('id', 'selected-molecule-overlay-svg');
            svg.style.position = 'absolute';
            svg.style.left = '0';
            svg.style.top = '0';
            svg.style.width = '100%';
            svg.style.height = '100%';
            svg.style.pointerEvents = 'none';
            container.appendChild(svg);
            try { console.debug('[overlay] created overlay svg'); } catch (e) {}
        }
        while (svg.firstChild) svg.removeChild(svg.firstChild);
        const selectedSet = (window.State && State.get && State.get('selectedBits')) || new Set();
        if (!selectedSet || selectedSet.size === 0) { try { console.debug('[overlay] no selected bits; removing overlay'); } catch (e) {} svg.remove(); return; }
        const stateRes = (window.State && State.get && State.get('currentAnalysisResult')) || (window.currentAnalysisResult);
        const smiles = stateRes && (stateRes.smiles || stateRes.target_smiles);
        if (!smiles) { try { console.warn('[overlay] no smiles present'); } catch (e) {} return; }
        const bitEnvironments = stateRes && stateRes.bit_environments;
        if (!bitEnvironments || typeof bitEnvironments !== 'object') {
            try { console.warn('[overlay] no bit_environments available'); } catch (e) {}
            return;
        }
        // Compute drawing area based on actual rendered molecule element bounds (to avoid padding/aspect mismatch)
        const molEl = container.querySelector('svg, img');
        const cRect = container.getBoundingClientRect();
        let offsetX = 0, offsetY = 0, drawW = container.clientWidth, drawH = container.clientHeight;
        if (molEl && molEl.tagName.toLowerCase() === 'svg') {
            const svgRect = molEl.getBoundingClientRect();
            const viewBox = (molEl.getAttribute('viewBox') || '').split(/\s+/).map(Number);
            let vbX=0, vbY=0, vbW=svgRect.width, vbH=svgRect.height;
            if (viewBox.length === 4 && viewBox.every(v => !isNaN(v))) {
                [vbX, vbY, vbW, vbH] = viewBox;
            }
            // Union bbox of all paths (bonds) to estimate content box (padding-aware)
            let minX=Infinity, minY=Infinity, maxX=-Infinity, maxY=-Infinity;
            const paths = molEl.querySelectorAll('path');
            paths.forEach(p => {
                try {
                    const b = p.getBBox();
                    if (b && isFinite(b.x) && isFinite(b.y) && isFinite(b.width) && isFinite(b.height)) {
                        minX = Math.min(minX, b.x);
                        minY = Math.min(minY, b.y);
                        maxX = Math.max(maxX, b.x + b.width);
                        maxY = Math.max(maxY, b.y + b.height);
                    }
                } catch (e) {}
            });
            if (isFinite(minX) && isFinite(minY) && isFinite(maxX) && isFinite(maxY)) {
                const pxPerUnitX = svgRect.width / vbW;
                const pxPerUnitY = svgRect.height / vbH;
                offsetX = (svgRect.left - cRect.left) + (minX - vbX) * pxPerUnitX;
                offsetY = (svgRect.top - cRect.top) + (minY - vbY) * pxPerUnitY;
                drawW = Math.max(1, (maxX - minX) * pxPerUnitX);
                drawH = Math.max(1, (maxY - minY) * pxPerUnitY);
            } else {
                offsetX = Math.max(0, svgRect.left - cRect.left);
                offsetY = Math.max(0, svgRect.top - cRect.top);
                drawW = Math.max(1, Math.floor(svgRect.width));
                drawH = Math.max(1, Math.floor(svgRect.height));
            }
        } else if (molEl) {
            const mRect = molEl.getBoundingClientRect();
            offsetX = Math.max(0, mRect.left - cRect.left);
            offsetY = Math.max(0, mRect.top - cRect.top);
            drawW = Math.max(1, Math.floor(mRect.width));
            drawH = Math.max(1, Math.floor(mRect.height));
        }
        try { console.debug('[overlay] drawing with dims', { w: drawW, h: drawH, offsetX, offsetY, selectedBits: Array.from(selectedSet) }); } catch (e) {}
        const drawEnv = (env) => {
            if (!(env && env.coords && Array.isArray(env.coords.atoms))) return;
            const coords = env.coords.atoms;
            const pos = (i) => {
                const c = coords.find(a => a.id === i);
                if (!c) return null;
                return { x: offsetX + (c.x * drawW), y: offsetY + ((1 - c.y) * drawH) };
            };
            (env.bonds || []).forEach(pair => {
                const a = pair[0], b = pair[1];
                const pa = pos(a), pb = pos(b);
                if (!pa || !pb) return;
                const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                line.setAttribute('x1', String(pa.x));
                line.setAttribute('y1', String(pa.y));
                line.setAttribute('x2', String(pb.x));
                line.setAttribute('y2', String(pb.y));
                line.setAttribute('stroke', 'rgba(16,185,129,0.7)');
                line.setAttribute('stroke-width', '4');
                svg.appendChild(line);
            });
            (env.atoms || []).forEach(a => {
                const p = pos(a);
                if (!p) return;
                const circ = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
                circ.setAttribute('cx', String(p.x));
                circ.setAttribute('cy', String(p.y));
                circ.setAttribute('r', '10');
                circ.setAttribute('fill', 'rgba(16,185,129,0.25)');
                circ.setAttribute('stroke', 'rgba(16,185,129,0.7)');
                circ.setAttribute('stroke-width', '2');
                svg.appendChild(circ);
            });
        };
        // For each selected bit, read env from bit_environments and draw
        const bits = Array.from(selectedSet);
        bits.forEach(bit => {
            const env = bitEnvironments[bit];
            if (env) {
                try { console.debug('[overlay] using env for bit', bit); } catch (e) {}
                drawEnv(env);
            } else {
                try { console.debug('[overlay] no env available for bit', bit); } catch (e) {}
            }
        });
    };

    // Rebuild overlay on resize (debounced)
    try {
        if (!window._overlayResizeBound) {
            window._overlayResizeBound = true;
            let _t = null;
            window.addEventListener('resize', () => {
                if (_t) clearTimeout(_t);
                _t = setTimeout(() => {
                    if (typeof rebuildSelectedOverlay === 'function') rebuildSelectedOverlay();
                }, 150);
            });
        }
    } catch (e) {}

    const buildIndicesList = (container, indices, typeKey) => {
        if (!container) {
            try { console.warn('[buildIndicesList] container not found for', typeKey); } catch (e) {}
            return;
        }
        container.textContent = '';
        try { console.debug('[buildIndicesList] building', typeKey, { indicesCount: (indices||[]).length, retrievedIndicesCount: (retrievedIndices||[]).length }); } catch (e) {}
        const smilesForPreview = (window.currentAnalysisResult && window.currentAnalysisResult.smiles) || (window.State && State.get && State.get('currentAnalysisResult') && State.get('currentAnalysisResult').smiles) || null;
        const presentSet = new Set(indices || []);
        const selectedSet = (window.State && State.get && State.get('selectedBits')) || new Set();
        if (window.State && State.set) State.set('selectedBits', selectedSet);
        const isPredicted = typeKey === 'predicted';
        // Use retrievedIndices from closure (updateFingerprintIndices parameter), not from State
        // This ensures we use the current molecule's indices, not stale data from previous analysis
        const retrievedSet = new Set(retrievedIndices);
        // Update State when building retrieved list, and when building predicted list if State is empty
        if (typeKey === 'retrieved') {
            // Update State with current molecule's retrieved indices for restore logic
            if (window.State && State.set) {
                State.set('retrievedFpIndices', retrievedIndices);
            }
        } else if (isPredicted && window.State && !State.get('retrievedFpIndices')) {
            // Cache retrieved indices for disable logic if not already set
            State.set('retrievedFpIndices', retrievedIndices);
        }
        if (indices && Array.isArray(indices) && indices.length > 0) {
            const sorted = [...indices].sort((a, b) => a - b);
            const count = sorted.length;
            const previewCount = Math.min(50, count);
            const preview = sorted.slice(0, previewCount);
            const hasMore = count > previewCount;
            
            const summary = document.createElement('div');
            summary.className = 'indices-summary';
            summary.textContent = `Total: ${count} bits`;
            container.appendChild(summary);
            
            const valuesDiv = document.createElement('div');
            valuesDiv.className = 'indices-values';
            const makeBadge = (idx) => {
                const badge = document.createElement('span');
                badge.className = 'index-badge';
                badge.textContent = String(idx);
                const isInRetrieved = retrievedSet.has(idx);
                if (isPredicted && !isInRetrieved) {
                    badge.classList.add('badge-disabled');
                    badge.title = 'Not present in selected molecule';
                } else {
                    if (selectedSet.has(idx)) badge.classList.add('badge-selected');
                    badge.addEventListener('click', () => {
                        try { console.debug('[bits] click', { type: typeKey, bit: idx }); } catch (e) {}
                        // Check if bit has env available (from embedded data)
                        const stateRes = (window.State && State.get && State.get('currentAnalysisResult')) || (window.currentAnalysisResult);
                        const bitEnvs = stateRes && stateRes.bit_environments;
                        const hasEnv = bitEnvs && bitEnvs[idx];
                        if (!hasEnv) {
                            // Bit has no env - mark as unavailable but don't allow selection
                            badge.classList.add('badge-unavailable');
                            badge.title = 'No substructure mapping available for this bit';
                            return;
                        }
                        const currentlySelected = selectedSet.has(idx);
                        if (currentlySelected) {
                            selectedSet.delete(idx);
                            badge.classList.remove('badge-selected');
                        } else {
                            selectedSet.add(idx);
                            badge.classList.add('badge-selected');
                        }
                        if (window.State && State.set) State.set('selectedBits', selectedSet);
                        // Mirror selection across both lists if bit exists in both
                        document.querySelectorAll(`.index-badge`).forEach(el => {
                            if (parseInt(el.textContent, 10) === idx) {
                                if (selectedSet.has(idx)) el.classList.add('badge-selected'); else el.classList.remove('badge-selected');
                            }
                        });
                        // Rebuild overlay for all selected bits
                        if (typeof rebuildSelectedOverlay === 'function') {
                            try { console.debug('[overlay] rebuilding after click'); } catch (e) {}
                            rebuildSelectedOverlay();
                        }
                    });
                }
                return badge;
            };
            preview.forEach((idx) => { valuesDiv.appendChild(makeBadge(idx)); });
            if (hasMore) {
                const more = document.createElement('span');
                more.className = 'index-more';
                more.textContent = `... and ${count - previewCount} more`;
                valuesDiv.appendChild(more);
            }
            container.appendChild(valuesDiv);
            
            if (hasMore) {
                const fullDiv = document.createElement('div');
                fullDiv.className = 'indices-full';
                fullDiv.style.display = 'none';
                sorted.forEach((idx) => { fullDiv.appendChild(makeBadge(idx)); });
                container.appendChild(fullDiv);
                
                const button = document.createElement('button');
                button.className = 'btn btn-sm btn-secondary';
                button.textContent = 'Show All';
                button.addEventListener('click', () => {
                    if (fullDiv.style.display === 'none') {
                        fullDiv.style.display = 'block';
                        valuesDiv.style.display = 'none';
                        button.textContent = 'Show Less';
                    } else {
                        fullDiv.style.display = 'none';
                        valuesDiv.style.display = 'block';
                        button.textContent = 'Show All';
                    }
                });
                container.appendChild(button);
            }
        } else {
            const empty = document.createElement('div');
            empty.className = 'no-indices';
            empty.textContent = 'No indices available';
            container.appendChild(empty);
        }
    };
    
    // Build indices lists (containers already verified above)
    buildIndicesList(predictedContainer, predictedIndices, 'predicted');
    buildIndicesList(retrievedContainer, retrievedIndices, 'retrieved');

    // Mark unavailable bits in yellow based on bit_environments
    try {
        const stateRes = (window.State && State.get && State.get('currentAnalysisResult')) || (window.currentAnalysisResult);
        const bitEnvs = stateRes && stateRes.bit_environments;
        if (bitEnvs && typeof bitEnvs === 'object') {
            const allBits = new Set([...(predictedIndices||[]), ...(retrievedIndices||[])]);
            allBits.forEach(bit => {
                if (!bitEnvs[bit]) {
                    // Mark badges for this bit as unavailable (yellow)
                    document.querySelectorAll('.index-badge').forEach(el => {
                        if (parseInt(el.textContent, 10) === bit) {
                            el.classList.add('badge-unavailable');
                            el.title = 'No substructure mapping available for this bit';
                        }
                    });
                }
            });
        }
        // Persist analysis state
        try {
            if (window.State && State.persistAnalysis) {
                const selected = Array.from((window.State && State.get && State.get('selectedBits')) || new Set());
                State.persistAnalysis({
                    currentAnalysisResult: stateRes,
                    predictedFpIndices: predictedIndices || [],
                    retrievedFpIndices: retrievedIndices || [],
                    selectedBits: selected
                });
            }
        } catch (e) {}
    } catch (e) {}
}

async function updateSimilarityVisualization(imageData) {
    const container = document.getElementById('similarity-visualization');
    
    container.textContent = '';
    if (imageData) {
        const img = document.createElement('img');
        img.src = imageData;
        img.alt = 'Similarity highlighting';
        img.className = 'molecule-image';
        container.appendChild(img);
        
        if (img && typeof panzoom !== 'undefined') {
            panzoom(img, {
                maxZoom: 5,
                minZoom: 0.5,
                initialZoom: 1
            });
        }
    } else {
        const noMol = document.createElement('div');
        noMol.className = 'no-molecule';
        noMol.textContent = 'Failed to generate similarity visualization';
        container.appendChild(noMol);
    }
}

async function updateChangeVisualization(imageData) {
    const container = document.getElementById('change-visualization');
    
    container.textContent = '';
    if (imageData) {
        const img = document.createElement('img');
        img.src = imageData;
        img.alt = 'Change highlighting';
        img.className = 'molecule-image';
        container.appendChild(img);
        
        if (img && typeof panzoom !== 'undefined') {
            panzoom(img, {
                maxZoom: 5,
                minZoom: 0.5,
                initialZoom: 1
            });
        }
    } else {
        const noMol = document.createElement('div');
        noMol.className = 'no-molecule';
        noMol.textContent = 'Failed to generate change visualization';
        container.appendChild(noMol);
    }
}

function updateFingerprintDifferences(differences) {
    
    const fingerprintSection = document.querySelector('.fingerprint-differences-section');
    if (fingerprintSection) {
        fingerprintSection.style.display = 'block';
    }
    
    const addedContainer = document.getElementById('added-bits-list');
    const removedContainer = document.getElementById('removed-bits-list');
    
    // Added bits
    addedContainer.textContent = '';
    if (differences.added && differences.added.length > 0) {
        differences.added.forEach(bit => {
            const item = document.createElement('div');
            item.className = 'bit-item';
            const header = document.createElement('div');
            header.className = 'bit-header';
            header.textContent = `Bit ${bit.bit_id}`;
            const details = document.createElement('div');
            details.className = 'bit-details';
            details.textContent = `Atom: ${bit.atom_symbol} | Radius: ${bit.radius}`;
            const frag = document.createElement('div');
            frag.className = 'fragment-structure';
            if (bit.substructure_image) {
                const img = document.createElement('img');
                img.src = bit.substructure_image;
                img.alt = 'Substructure';
                img.className = 'substructure-image';
                frag.appendChild(img);
            }
            const smiles = document.createElement('span');
            smiles.className = 'smiles-text';
            smiles.textContent = bit.fragment_smiles || '';
            frag.appendChild(smiles);
            item.appendChild(header);
            item.appendChild(details);
            item.appendChild(frag);
            addedContainer.appendChild(item);
        });
    } else {
        const empty = document.createElement('div');
        empty.className = 'no-differences';
        empty.textContent = 'No added bits detected';
        addedContainer.appendChild(empty);
    }
    
    // Removed bits
    removedContainer.textContent = '';
    if (differences.removed && differences.removed.length > 0) {
        differences.removed.forEach(bit => {
            const item = document.createElement('div');
            item.className = 'bit-item';
            const header = document.createElement('div');
            header.className = 'bit-header';
            header.textContent = `Bit ${bit.bit_id}`;
            const details = document.createElement('div');
            details.className = 'bit-details';
            details.textContent = `Atom: ${bit.atom_symbol} | Radius: ${bit.radius}`;
            const frag = document.createElement('div');
            frag.className = 'fragment-structure';
            if (bit.substructure_image) {
                const img = document.createElement('img');
                img.src = bit.substructure_image;
                img.alt = 'Substructure';
                img.className = 'substructure-image';
                frag.appendChild(img);
            }
            const smiles = document.createElement('span');
            smiles.className = 'smiles-text';
            smiles.textContent = bit.fragment_smiles || '';
            frag.appendChild(smiles);
            item.appendChild(header);
            item.appendChild(details);
            item.appendChild(frag);
            removedContainer.appendChild(item);
        });
    } else {
        const empty = document.createElement('div');
        empty.className = 'no-differences';
        empty.textContent = 'No removed bits detected';
        removedContainer.appendChild(empty);
    }
}

// Secondary retrieval function
async function runSecondaryRetrieval(predictedFp, retrievedFp, k = 10) {
    if (!predictedFp || !retrievedFp) {
        console.warn('Cannot run secondary retrieval: fingerprints not available');
        return;
    }
    
    // Show loading state
    const secondarySection = document.getElementById('secondary-retrieval-section');
    if (secondarySection) {
        secondarySection.style.display = 'block';
    }
    
    const secondaryLoading = document.getElementById('secondary-retrieval-loading');
    const secondaryGrid = document.getElementById('secondary-results-grid');
    
    if (secondaryLoading) secondaryLoading.style.display = 'flex';
    if (secondaryGrid) secondaryGrid.innerHTML = '';
    
    try {
        const response = await fetch('/secondary-retrieval', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({
                predicted_fp: predictedFp,
                retrieved_fp: retrievedFp,
                k: k
            })
        });
        
        if (!response.ok) {
            let errorMessage = 'Secondary retrieval failed.';
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
        if (secondaryLoading) secondaryLoading.style.display = 'none';
        
        // Display results using the same display format as main page
        if (result.results && Array.isArray(result.results)) {
            if (secondaryGrid) {
                secondaryGrid.innerHTML = '';
                if (result.results.length === 0) {
                    secondaryGrid.innerHTML = '<p class="no-results-message">No results found for remaining substructure.</p>';
                } else {
                    // Store results for potential analysis (same format as main page)
                    const secondaryResults = result.results;
                    for (let i = 0; i < secondaryResults.length; i++) {
                        const resultData = secondaryResults[i];
                        const card = (window.ResultRenderer && window.ResultRenderer.createCard)
                          ? window.ResultRenderer.createCard(i + 1, resultData, { showAnalyze: false })
                          : createResultCard(i + 1, resultData);
                        secondaryGrid.appendChild(card);
                    }
                }
            }
        } else {
            if (secondaryGrid) {
                secondaryGrid.innerHTML = '<p class="error-message">No results returned from secondary retrieval.</p>';
            }
        }
        
    } catch (error) {
        console.error('Secondary retrieval error:', error);
        if (secondaryLoading) secondaryLoading.style.display = 'none';
        if (secondaryGrid) {
            const p = document.createElement('p');
            p.className = 'error-message';
            p.textContent = `Secondary retrieval failed: ${error.message}`;
            secondaryGrid.textContent = '';
            secondaryGrid.appendChild(p);
        }
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
    // Restore main page state if available
    try {
        if (window.State && State.restoreMain) {
            const snap = State.restoreMain();
            if (snap) {
                // Restore inputs
                if (snap.inputs && typeof loadDataIntoTable === 'function') {
                    if (snap.inputs.hsqc && Array.isArray(snap.inputs.hsqc)) {
                        loadDataIntoTable('hsqc', snap.inputs.hsqc, 3);
                    }
                    if (snap.inputs.h_nmr && Array.isArray(snap.inputs.h_nmr)) {
                        loadDataIntoTable('h_nmr', snap.inputs.h_nmr, 1);
                    }
                    if (snap.inputs.c_nmr && Array.isArray(snap.inputs.c_nmr)) {
                        loadDataIntoTable('c_nmr', snap.inputs.c_nmr, 1);
                    }
                    if (snap.inputs.mass_spec && Array.isArray(snap.inputs.mass_spec)) {
                        loadDataIntoTable('mass_spec', snap.inputs.mass_spec, 2);
                    }
                }
                const mwInput = document.getElementById('mw-input');
                if (mwInput && snap.inputs && snap.inputs.mw) {
                    mwInput.value = snap.inputs.mw;
                }
                // Restore SMILES input
                const smilesInput = document.getElementById('smiles-input');
                if (smilesInput && snap.smilesInput) {
                    smilesInput.value = snap.smilesInput;
                }
                // Restore results if available
                if (snap.currentResults && Array.isArray(snap.currentResults) && snap.currentResults.length > 0) {
                    if (window.Results && Results.displayResults) {
                        window.Results.displayResults({ results: snap.currentResults });
                    }
                    if (window.State && State.set) {
                        State.set('currentResults', snap.currentResults);
                    }
                }
            }
        }
    } catch (e) {
        console.warn('Failed to restore main page state:', e);
    }
    
    // Check backend status after a short delay to allow page to render first
    setTimeout(() => {
        if (window.Health && Health.updateBackendStatus) Health.updateBackendStatus();
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
    if (window.Health && typeof window.Health.startPolling === 'function') {
        window.Health.startPolling(10000);
    } else {
    const statusCheckInterval = setInterval(async () => {
        try {
                if (window.Health && Health.updateBackendStatus) await Health.updateBackendStatus();
            } catch (e) {
                if (window.Health && Health.updateBackendStatus) Health.updateBackendStatus();
        }
    }, 10000);
    }
});