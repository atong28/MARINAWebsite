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
        
        // Make all cells focusable
        this.makeCellsFocusable();
        
        // Add event listeners
        this.addEventListeners();
        
        // Add keyboard shortcuts
        this.addKeyboardShortcuts();
    }
    
    makeCellsFocusable() {
        const cells = this.table.querySelectorAll('td');
        cells.forEach(cell => {
            if (!cell.querySelector('input')) {
                cell.setAttribute('tabindex', '0');
                cell.classList.add('spreadsheet-cell');
            } else {
                const input = cell.querySelector('input');
                input.classList.add('spreadsheet-input');
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
                        e.preventDefault();
                        this.pasteSelection();
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
            navigator.clipboard.writeText(text).then(() => {
                this.showMessage('Copied to clipboard', 'success');
            }).catch(() => {
                this.showMessage('Failed to copy to clipboard', 'error');
            });
        }
    }
    
    pasteSelection() {
        navigator.clipboard.readText().then(text => {
            this.pasteData(text);
        }).catch(() => {
            this.showMessage('Failed to read from clipboard', 'error');
        });
    }
    
    pasteData(text) {
        const lines = text.trim().split('\n');
        const data = lines.map(line => line.split('\t'));
        
        if (this.selectedRange) {
            // Paste into selected range
            const rows = this.table.querySelectorAll('tbody tr');
            for (let i = 0; i < data.length && (this.selectedRange.minRow + i) < rows.length; i++) {
                const row = rows[this.selectedRange.minRow + i];
                const cells = row.querySelectorAll('td');
                for (let j = 0; j < data[i].length && (this.selectedRange.minCol + j) < cells.length; j++) {
                    const cell = cells[this.selectedRange.minCol + j];
                    const input = cell.querySelector('input');
                    if (input && !input.closest('th')) { // Don't paste into header cells
                        input.value = data[i][j];
                    }
                }
            }
        } else if (this.currentCell) {
            // Paste starting from current cell
            const startPos = this.getCellPosition(this.currentCell);
            const rows = this.table.querySelectorAll('tbody tr');
            
            for (let i = 0; i < data.length && (startPos.row + i) < rows.length; i++) {
                const row = rows[startPos.row + i];
                const cells = row.querySelectorAll('td');
                for (let j = 0; j < data[i].length && (startPos.col + j) < cells.length; j++) {
                    const cell = cells[startPos.col + j];
                    const input = cell.querySelector('input');
                    if (input && !input.closest('th')) {
                        input.value = data[i][j];
                    }
                }
            }
        }
        
        this.showMessage('Data pasted', 'success');
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
    tableIds.forEach(tableId => {
        if (document.getElementById(tableId)) {
            new SpreadsheetTable(tableId);
        }
    });
});
