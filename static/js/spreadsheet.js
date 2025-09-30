// Spreadsheet-like functionality for data tables
class SpreadsheetTable {
    constructor(tableId) {
        this.table = document.getElementById(tableId);
        this.currentCell = null;
        this.selectedRange = null;
        this.isSelecting = false;
        this.startCell = null;
        
        this.init();
    }
    
    init() {
        // Add spreadsheet classes
        this.table.classList.add('spreadsheet-table');
        
        // Make table focusable
        this.table.setAttribute('tabindex', '0');
        
        // Make all cells focusable
        this.makeCellsFocusable();
        
        // Add event listeners
        this.addEventListeners();
        
        // Add keyboard shortcuts
        this.addKeyboardShortcuts();
    }
    
    makeCellsFocusable() {
        this.addPasteListenersToInputs();
    }
    
    addPasteListenersToInputs() {
        const inputs = this.table.querySelectorAll('input');
        inputs.forEach(input => {
            input.classList.add('spreadsheet-input');
            
                // Add paste event listener directly to each input
                input.addEventListener('paste', (e) => {
                    e.preventDefault();
                    e.stopPropagation();
                    const pastedText = e.clipboardData.getData('text/plain');
                    if (pastedText) {
                        // Temporarily disable input validation during paste
                        input.setAttribute('data-pasting', 'true');
                        this.pasteData(pastedText);
                        // Re-enable validation after a short delay
                        setTimeout(() => {
                            input.removeAttribute('data-pasting');
                        }, 100);
                    }
                });
        });
        
        // Also handle cells without inputs
        const cells = this.table.querySelectorAll('td');
        cells.forEach(cell => {
            if (!cell.querySelector('input')) {
                cell.setAttribute('tabindex', '0');
                cell.classList.add('spreadsheet-cell');
            }
        });
    }
    
    addEventListeners() {
        // Click to select cell
        this.table.addEventListener('click', (e) => {
            const cell = e.target.closest('td');
            if (cell) {
                this.selectCell(cell);
            }
        });
        
        // Double-click to edit
        this.table.addEventListener('dblclick', (e) => {
            const cell = e.target.closest('td');
            if (cell && cell.querySelector('input')) {
                const input = cell.querySelector('input');
                input.focus();
                input.select();
            }
        });
        
        // Mouse down for range selection
        this.table.addEventListener('mousedown', (e) => {
            const cell = e.target.closest('td');
            if (cell) {
                this.isSelecting = true;
                this.startCell = cell;
                this.selectCell(cell);
                e.preventDefault();
            }
        });
        
        // Mouse over for range selection
        this.table.addEventListener('mouseover', (e) => {
            if (this.isSelecting && this.startCell) {
                const cell = e.target.closest('td');
                if (cell) {
                    this.selectRange(this.startCell, cell);
                }
            }
        });
        
        // Mouse up to end selection
        this.table.addEventListener('mouseup', () => {
            this.isSelecting = false;
        });
        
        // Input focus handling
        this.table.addEventListener('focusin', (e) => {
            if (e.target.tagName === 'INPUT') {
                const cell = e.target.closest('td');
                this.selectCell(cell);
            }
        });
        
        // Prevent default tab behavior on inputs
        this.table.addEventListener('keydown', (e) => {
            if (e.target.tagName === 'INPUT' && e.key === 'Tab') {
                e.preventDefault();
                this.navigate(e.shiftKey ? 'left' : 'right');
            }
        });
        
        // Add paste event listener to all inputs in the table
        this.table.addEventListener('paste', (e) => {
            e.preventDefault();
            e.stopPropagation();
            const pastedText = e.clipboardData.getData('text/plain');
            if (pastedText) {
                // Temporarily disable validation on all inputs during paste
                const inputs = this.table.querySelectorAll('input');
                inputs.forEach(input => {
                    input.setAttribute('data-pasting', 'true');
                });
                
                this.pasteData(pastedText);
                
                // Re-enable validation after paste is complete
                setTimeout(() => {
                    inputs.forEach(input => {
                        input.removeAttribute('data-pasting');
                    });
                }, 100);
            }
        });
    }
    
    addKeyboardShortcuts() {
        document.addEventListener('keydown', (e) => {
            // Only handle shortcuts when table is focused or input is focused within table
            const activeElement = document.activeElement;
            const isTableFocused = this.table.contains(activeElement) || 
                                 (activeElement.tagName === 'INPUT' && this.table.contains(activeElement));
            
            if (!isTableFocused) return;
            
            // Prevent default for navigation keys
            if (['ArrowUp', 'ArrowDown', 'ArrowLeft', 'ArrowRight', 'Tab'].includes(e.key)) {
                e.preventDefault();
            }
            
            switch (e.key) {
                case 'ArrowUp':
                    this.navigate('up');
                    break;
                case 'ArrowDown':
                    this.navigate('down');
                    break;
                case 'ArrowLeft':
                    this.navigate('left');
                    break;
                case 'ArrowRight':
                    this.navigate('right');
                    break;
                case 'Tab':
                    this.navigate(e.shiftKey ? 'left' : 'right');
                    break;
                case 'Enter':
                    this.navigate('down');
                    break;
                case 'Escape':
                    this.clearSelection();
                    break;
                case 'Delete':
                case 'Backspace':
                    this.clearCell();
                    break;
                case 'a':
                    if (e.ctrlKey || e.metaKey) {
                        e.preventDefault();
                        this.selectAll();
                    }
                    break;
                case 'c':
                    if (e.ctrlKey || e.metaKey) {
                        e.preventDefault();
                        this.copySelection();
                    }
                    break;
                case 'v':
                    if (e.ctrlKey || e.metaKey) {
                        // Let the paste event handlers handle this
                        // No need to prevent default or call pasteSelection
                    }
                    break;
            }
        });
    }
    
    selectCell(cell) {
        // Clear previous selection
        this.clearSelection();
        
        // Select new cell
        cell.classList.add('selected');
        this.currentCell = cell;
        
        // Focus the input if it exists, otherwise focus the cell
        const input = cell.querySelector('input');
        if (input) {
            input.focus();
        } else {
            cell.focus();
        }
    }
    
    selectRange(startCell, endCell) {
        this.clearSelection();
        
        const startPos = this.getCellPosition(startCell);
        const endPos = this.getCellPosition(endCell);
        
        const minRow = Math.min(startPos.row, endPos.row);
        const maxRow = Math.max(startPos.row, endPos.row);
        const minCol = Math.min(startPos.col, endPos.col);
        const maxCol = Math.max(startPos.col, endPos.col);
        
        const rows = this.table.querySelectorAll('tbody tr');
        for (let row = minRow; row <= maxRow; row++) {
            if (rows[row]) {
                const cells = rows[row].querySelectorAll('td');
                for (let col = minCol; col <= maxCol; col++) {
                    if (cells[col]) {
                        cells[col].classList.add('selected');
                    }
                }
            }
        }
        
        this.selectedRange = { startCell, endCell, minRow, maxRow, minCol, maxCol };
        this.currentCell = endCell;
    }
    
    getCellPosition(cell) {
        const row = cell.closest('tr');
        const rows = this.table.querySelectorAll('tbody tr');
        const cells = row.querySelectorAll('td');
        
        return {
            row: Array.from(rows).indexOf(row),
            col: Array.from(cells).indexOf(cell)
        };
    }
    
    navigate(direction) {
        if (!this.currentCell) return;
        
        const currentPos = this.getCellPosition(this.currentCell);
        const rows = this.table.querySelectorAll('tbody tr');
        const currentRowCells = rows[currentPos.row].querySelectorAll('td');
        
        let newRow = currentPos.row;
        let newCol = currentPos.col;
        
        switch (direction) {
            case 'up':
                newRow = Math.max(0, currentPos.row - 1);
                break;
            case 'down':
                newRow = Math.min(rows.length - 1, currentPos.row + 1);
                break;
            case 'left':
                newCol = Math.max(0, currentPos.col - 1);
                break;
            case 'right':
                newCol = Math.min(currentRowCells.length - 1, currentPos.col + 1);
                break;
        }
        
        // Handle row changes for left/right navigation
        if (direction === 'left' || direction === 'right') {
            const targetRowCells = rows[newRow].querySelectorAll('td');
            newCol = Math.min(newCol, targetRowCells.length - 1);
        }
        
        const targetRow = rows[newRow];
        const targetCells = targetRow.querySelectorAll('td');
        const targetCell = targetCells[newCol];
        
        if (targetCell) {
            this.selectCell(targetCell);
        }
    }
    
    clearSelection() {
        const selectedCells = this.table.querySelectorAll('.selected');
        selectedCells.forEach(cell => cell.classList.remove('selected'));
        this.selectedRange = null;
    }
    
    clearCell() {
        if (this.selectedRange) {
            // Clear all selected cells
            const selectedCells = this.table.querySelectorAll('.selected');
            selectedCells.forEach(cell => {
                const input = cell.querySelector('input');
                if (input) {
                    input.value = '';
                }
            });
        } else if (this.currentCell) {
            // Clear current cell
            const input = this.currentCell.querySelector('input');
            if (input) {
                input.value = '';
            }
        }
    }
    
    selectAll() {
        const firstRow = this.table.querySelector('tbody tr');
        const lastRow = this.table.querySelector('tbody tr:last-child');
        if (firstRow && lastRow) {
            this.selectRange(firstRow.querySelector('td'), lastRow.querySelector('td:last-child'));
        }
    }
    
    copySelection() {
        let text = '';
        
        if (this.selectedRange) {
            // Copy selected range
            const rows = this.table.querySelectorAll('tbody tr');
            for (let row = this.selectedRange.minRow; row <= this.selectedRange.maxRow; row++) {
                const rowData = [];
                const cells = rows[row].querySelectorAll('td');
                for (let col = this.selectedRange.minCol; col <= this.selectedRange.maxCol; col++) {
                    const cell = cells[col];
                    const input = cell.querySelector('input');
                    rowData.push(input ? input.value : cell.textContent.trim());
                }
                text += rowData.join('\t');
                if (row < this.selectedRange.maxRow) {
                    text += '\n';
                }
            }
        } else if (this.currentCell) {
            // Copy single cell
            const input = this.currentCell.querySelector('input');
            text = input ? input.value : this.currentCell.textContent.trim();
        }
        
        if (text) {
            // Use the reliable document.execCommand approach first
            this.showFallbackCopyDialog(text);
        }
    }
    
    
    // Test function to manually trigger paste with sample data
    testPaste() {
        const testData = "7.2\t120.5\t1.0\n7.5\t125.0\t0.8\n8.1\t130.2\t1.2";
        this.pasteData(testData);
    }
    
    showFallbackCopyDialog(text) {
        // Create a temporary textarea to select and copy text
        const textarea = document.createElement('textarea');
        textarea.value = text;
        textarea.style.position = 'fixed';
        textarea.style.left = '-999999px';
        textarea.style.top = '-999999px';
        document.body.appendChild(textarea);
        textarea.focus();
        textarea.select();
        
        try {
            const successful = document.execCommand('copy');
            // Silent copy - no popup messages
        } catch (err) {
            // Silent fallback - no popup messages
        }
        
        document.body.removeChild(textarea);
    }
    
    pasteData(text) {
        const lines = text.trim().split('\n');
        const data = lines.map(line => line.split(/[\t,;]/).map(v => v.trim().replace(/\r/g, ''))); // Handle various separators and remove \r
        
        const tableId = this.table.id.replace('-table', '');
        const tbody = this.table.querySelector('tbody');
        
        // Determine how many columns this table expects
        let expectedCols;
        switch (tableId) {
            case 'hsqc':
                expectedCols = 3; // H, C, intensity
                break;
            case 'mass_spec':
                expectedCols = 2; // m/z, intensity
                break;
            case 'h_nmr':
            case 'c_nmr':
                expectedCols = 1; // single value
                break;
            default:
                expectedCols = 1;
        }
        
        // Process data and create rows as needed
        let currentRowIndex = 0;
        let currentColIndex = 0;
        
        if (this.selectedRange) {
            currentRowIndex = this.selectedRange.minRow;
            currentColIndex = this.selectedRange.minCol;
        } else if (this.currentCell) {
            const pos = this.getCellPosition(this.currentCell);
            currentRowIndex = pos.row;
            currentColIndex = pos.col;
        }
        
        // Flatten data for easier processing
        const flatData = [];
        data.forEach(row => {
            if (tableId === 'hsqc') {
                // For HSQC, group every 3 values into a row
                for (let i = 0; i < row.length; i += 3) {
                    if (i + 2 < row.length) {
                        flatData.push([row[i], row[i + 1], row[i + 2]]);
                    }
                }
            } else if (tableId === 'mass_spec') {
                // For mass spec, group every 2 values into a row
                for (let i = 0; i < row.length; i += 2) {
                    if (i + 1 < row.length) {
                        flatData.push([row[i], row[i + 1]]);
                    }
                }
            } else {
                // For single column data, each value is a row
                row.forEach(value => {
                    if (value) flatData.push([value]);
                });
            }
        });
        
        // Ensure we have enough rows
        while (tbody.children.length <= currentRowIndex + flatData.length - 1) {
            // Add new row using the global addRow function
            if (typeof addRow === 'function') {
                addRow(tableId);
            } else {
                // Fallback: create a simple row
                const newRow = document.createElement('tr');
                for (let i = 0; i < expectedCols; i++) {
                    const cell = document.createElement('td');
                    const input = document.createElement('input');
                    input.type = 'number';
                    input.step = '0.01';
                    input.className = 'spreadsheet-input small';
                    cell.appendChild(input);
                    newRow.appendChild(cell);
                }
                // Add remove button
                const removeCell = document.createElement('td');
                const removeBtn = document.createElement('button');
                removeBtn.className = 'btn-icon';
                removeBtn.onclick = function() { removeRow(this); };
                removeBtn.title = 'Remove';
                removeBtn.innerHTML = '<i class="fas fa-times"></i>';
                removeCell.appendChild(removeBtn);
                newRow.appendChild(removeCell);
                tbody.appendChild(newRow);
            }
        }
        
        // Fill the data
        flatData.forEach((rowData, dataIndex) => {
            const targetRowIndex = currentRowIndex + dataIndex;
            if (targetRowIndex < tbody.children.length) {
                const row = tbody.children[targetRowIndex];
                const inputs = row.querySelectorAll('input');
                
                rowData.forEach((value, colIndex) => {
                    if (inputs[colIndex] && value) {
                        // Clean the value and validate it's a number
                        const cleanValue = value.replace(/[^\d.-]/g, '');
                        if (cleanValue && !isNaN(parseFloat(cleanValue))) {
                            // Set value directly without triggering validation
                            const input = inputs[colIndex];
                            input.setAttribute('data-pasting', 'true');
                            input.value = cleanValue;
                            // Trigger input event to update any dependent logic
                            input.dispatchEvent(new Event('input', { bubbles: true }));
                        }
                    }
                });
            }
        });
        
        // Re-add paste listeners to any new inputs that were created
        this.addPasteListenersToInputs();
        
        // Clean up data-pasting attributes after paste is complete
        setTimeout(() => {
            const allInputs = this.table.querySelectorAll('input[data-pasting]');
            allInputs.forEach(input => {
                input.removeAttribute('data-pasting');
            });
        }, 200);
        
        // Update input summary if the function exists
        if (typeof updateInputSummary === 'function') {
            updateInputSummary();
        }
        
        // Update visualizations after paste
        if (typeof updateVisualizations === 'function') {
            updateVisualizations();
        }
    }
    
    showMessage(message, type) {
        // Create a temporary message element
        const messageEl = document.createElement('div');
        messageEl.className = `spreadsheet-message ${type}`;
        messageEl.textContent = message;
        messageEl.style.cssText = `
            position: fixed;
            top: 20px;
            right: 20px;
            padding: 10px 20px;
            border-radius: 4px;
            color: white;
            font-weight: bold;
            z-index: 1000;
            opacity: 0;
            transition: opacity 0.3s ease;
        `;
        
        if (type === 'success') {
            messageEl.style.backgroundColor = '#10b981';
        } else if (type === 'error') {
            messageEl.style.backgroundColor = '#ef4444';
        } else if (type === 'info') {
            messageEl.style.backgroundColor = '#3b82f6';
        }
        
        document.body.appendChild(messageEl);
        
        // Fade in
        setTimeout(() => {
            messageEl.style.opacity = '1';
        }, 10);
        
        // Remove after 3 seconds
        setTimeout(() => {
            messageEl.style.opacity = '0';
            setTimeout(() => {
                document.body.removeChild(messageEl);
            }, 300);
        }, 3000);
    }
}

// Initialize spreadsheet functionality for all tables
document.addEventListener('DOMContentLoaded', function() {
    // Initialize spreadsheet for all data tables
    const tableIds = ['hsqc-table', 'h_nmr-table', 'c_nmr-table', 'mass_spec-table'];
    const spreadsheetTables = [];
    
    tableIds.forEach(tableId => {
        if (document.getElementById(tableId)) {
            const table = new SpreadsheetTable(tableId);
            spreadsheetTables.push(table);
            
            // Make test function globally accessible
            window[`testPaste_${tableId.replace('-table', '')}`] = () => table.testPaste();
        }
    });
    
    // Add document-level paste handler for better coverage
    document.addEventListener('paste', function(e) {
        // Check if the paste is happening in one of our tables
        const activeElement = document.activeElement;
        let targetTable = null;
        
        // First check if the active element is directly in a table
        for (const table of spreadsheetTables) {
            if (table.table.contains(activeElement)) {
                targetTable = table;
                break;
            }
        }
        
        // If no direct match, check if we're in an input within a table
        if (!targetTable && activeElement.tagName === 'INPUT') {
            for (const table of spreadsheetTables) {
                if (table.table.contains(activeElement)) {
                    targetTable = table;
                    break;
                }
            }
        }
        
        // If still no match, check if the event target is in a table
        if (!targetTable) {
            for (const table of spreadsheetTables) {
                if (table.table.contains(e.target)) {
                    targetTable = table;
                    break;
                }
            }
        }
        
        if (targetTable) {
            e.preventDefault();
            e.stopPropagation();
            const pastedText = e.clipboardData.getData('text/plain');
            if (pastedText) {
                // Temporarily disable validation on all inputs during paste
                const inputs = targetTable.table.querySelectorAll('input');
                inputs.forEach(input => {
                    input.setAttribute('data-pasting', 'true');
                });
                
                targetTable.pasteData(pastedText);
                
                // Re-enable validation after paste is complete
                setTimeout(() => {
                    inputs.forEach(input => {
                        input.removeAttribute('data-pasting');
                    });
                }, 100);
            }
        }
    });
});
