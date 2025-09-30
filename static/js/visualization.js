// Client-side spectral visualization functionality
class SpectrumVisualizer {
    constructor(canvasId, type) {
        this.canvas = document.getElementById(canvasId);
        if (!this.canvas) {
            console.warn(`Canvas element with id '${canvasId}' not found`);
            return;
        }
        this.ctx = this.canvas.getContext('2d');
        this.type = type;
        this.data = [];
        this.isHsqcSwapped = false;
    }

    updateData(data) {
        this.data = data;
        this.draw();
    }

    draw() {
        const { ctx, canvas } = this;
        const width = canvas.width;
        const height = canvas.height;
        
        // Clear canvas
        ctx.clearRect(0, 0, width, height);
        
        if (this.data.length === 0) {
            this.drawEmpty();
            return;
        }

        switch (this.type) {
            case 'hsqc':
                this.drawHSQC();
                break;
            case 'h_nmr':
                this.drawH1D();
                break;
            case 'c_nmr':
                this.drawC1D();
                break;
            case 'mass_spec':
                this.drawMassSpec();
                break;
        }
    }

    drawEmpty() {
        const { ctx, canvas } = this;
        const width = canvas.width;
        const height = canvas.height;
        
        ctx.fillStyle = '#f8fafc';
        ctx.fillRect(0, 0, width, height);
        
        ctx.fillStyle = '#94a3b8';
        ctx.font = '16px -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText('No visualizations', width/2, height/2);
    }

    drawHSQC() {
        const { ctx, canvas } = this;
        const width = canvas.width;
        const height = canvas.height;
        const margin = 60;
        const plotWidth = width - 2 * margin;
        const plotHeight = height - 2 * margin;
        
        // Background
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, width, height);
        
        // Grid
        ctx.strokeStyle = '#f1f5f9';
        ctx.lineWidth = 1;
        for (let i = 0; i <= 10; i++) {
            const x = margin + (i * plotWidth / 10);
            const y = margin + (i * plotHeight / 10);
            
            ctx.beginPath();
            ctx.moveTo(x, margin);
            ctx.lineTo(x, margin + plotHeight);
            ctx.stroke();
            
            ctx.beginPath();
            ctx.moveTo(margin, y);
            ctx.lineTo(margin + plotWidth, y);
            ctx.stroke();
        }
        
        // Data points
        if (this.data.length >= 3 && this.data.length % 3 === 0) {
            const points = [];
            for (let i = 0; i < this.data.length; i += 3) {
                const hShift = this.data[i];
                const cShift = this.data[i + 1];
                const intensity = this.data[i + 2];
                points.push({ hShift, cShift, intensity });
            }
            
            // Find ranges
            const hShifts = points.map(p => p.hShift);
            const cShifts = points.map(p => p.cShift);
            const intensities = points.map(p => p.intensity);
            
            const hMin = Math.min(...hShifts);
            const hMax = Math.max(...hShifts);
            const cMin = Math.min(...cShifts);
            const cMax = Math.max(...cShifts);
            const intensityMax = Math.max(...intensities);
            
            // Draw points
            points.forEach(point => {
                // Map high values to left, low values to right
                const x = margin + ((hMax - point.hShift) / (hMax - hMin)) * plotWidth;
                // Map low C values to top, high C values to bottom (lowest at top, increasing downward)
                const y = margin + ((point.cShift - cMin) / (cMax - cMin)) * plotHeight;
                const radius = Math.max(3, (Math.abs(point.intensity) / intensityMax) * 12);
                
                // Use red for positive peaks, blue for negative peaks
                ctx.fillStyle = point.intensity >= 0 ? '#ef4444' : '#3b82f6';
                ctx.beginPath();
                ctx.arc(x, y, radius, 0, 2 * Math.PI);
                ctx.fill();
            });
            
            // Draw axis numbers
            ctx.fillStyle = '#374151';
            ctx.font = '12px -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif';
            ctx.textAlign = 'center';
            
            // X-axis numbers (<sup>1</sup>H shifts) - high to low (right to left)
            for (let i = 0; i <= 5; i++) {
                const x = margin + (i * plotWidth / 5);
                const value = hMax - (i * (hMax - hMin) / 5);
                ctx.fillText(value.toFixed(1), x, height - 15);
            }
            
            // Y-axis numbers (¹³C shifts) - low to high (top to bottom)
            ctx.textAlign = 'right';
            for (let i = 0; i <= 5; i++) {
                const y = margin + (i * plotHeight / 5);
                const value = cMin + (i * (cMax - cMin) / 5);
                ctx.fillText(value.toFixed(0), margin - 10, y + 4);
            }
        }
        
        // Axis labels
        ctx.fillStyle = '#374151';
        ctx.font = '14px -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('¹H Shift (ppm)', width/2, height - 5);
        
        ctx.save();
        ctx.translate(15, height/2);
        ctx.rotate(-Math.PI/2);
        ctx.fillText('¹³C Shift (ppm)', 0, 0);
        ctx.restore();
    }

    drawH1D() {
        const { ctx, canvas } = this;
        const width = canvas.width;
        const height = canvas.height;
        const margin = 60;
        const plotWidth = width - 2 * margin;
        const plotHeight = height - 2 * margin;
        
        // Background
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, width, height);
        
        // Grid
        ctx.strokeStyle = '#f1f5f9';
        ctx.lineWidth = 1;
        for (let i = 0; i <= 10; i++) {
            const x = margin + (i * plotWidth / 10);
            ctx.beginPath();
            ctx.moveTo(x, margin);
            ctx.lineTo(x, margin + plotHeight);
            ctx.stroke();
        }
        
        // Data points
        if (this.data.length > 0) {
            const maxShift = Math.max(...this.data);
            const minShift = Math.min(...this.data);
            const range = maxShift - minShift || 1;
            
            this.data.forEach((shift, index) => {
                // Map high values to left, low values to right
                const x = margin + ((maxShift - shift) / range) * plotWidth;
                
                ctx.strokeStyle = '#ef4444'; // Red for positive peaks
                ctx.lineWidth = 3;
                ctx.beginPath();
                ctx.moveTo(x, margin + plotHeight);
                ctx.lineTo(x, margin + plotHeight - 40);
                ctx.stroke();
            });
            
            // Draw axis numbers
            ctx.fillStyle = '#374151';
            ctx.font = '12px -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif';
            ctx.textAlign = 'center';
            
            for (let i = 0; i <= 5; i++) {
                const x = margin + (i * plotWidth / 5);
                const value = maxShift - (i * range / 5);
                ctx.fillText(value.toFixed(1), x, height - 15);
            }
        }
        
        // Axis label
        ctx.fillStyle = '#374151';
        ctx.font = '14px -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('¹H Shift (ppm)', width/2, height - 5);
    }

    drawC1D() {
        const { ctx, canvas } = this;
        const width = canvas.width;
        const height = canvas.height;
        const margin = 60;
        const plotWidth = width - 2 * margin;
        const plotHeight = height - 2 * margin;
        
        // Background
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, width, height);
        
        // Grid
        ctx.strokeStyle = '#f1f5f9';
        ctx.lineWidth = 1;
        for (let i = 0; i <= 10; i++) {
            const x = margin + (i * plotWidth / 10);
            ctx.beginPath();
            ctx.moveTo(x, margin);
            ctx.lineTo(x, margin + plotHeight);
            ctx.stroke();
        }
        
        // Data points
        if (this.data.length > 0) {
            const maxShift = Math.max(...this.data);
            const minShift = Math.min(...this.data);
            const range = maxShift - minShift || 1;
            
            this.data.forEach((shift, index) => {
                // Map high values to left, low values to right
                const x = margin + ((maxShift - shift) / range) * plotWidth;
                
                ctx.strokeStyle = '#3b82f6'; // Blue for C NMR
                ctx.lineWidth = 3;
                ctx.beginPath();
                ctx.moveTo(x, margin + plotHeight);
                ctx.lineTo(x, margin + plotHeight - 40);
                ctx.stroke();
            });
            
            // Draw axis numbers
            ctx.fillStyle = '#374151';
            ctx.font = '12px -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif';
            ctx.textAlign = 'center';
            
            for (let i = 0; i <= 5; i++) {
                const x = margin + (i * plotWidth / 5);
                const value = maxShift - (i * range / 5);
                ctx.fillText(value.toFixed(0), x, height - 15);
            }
        }
        
        // Axis label
        ctx.fillStyle = '#374151';
        ctx.font = '14px -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('¹³C Shift (ppm)', width/2, height - 5);
    }

    drawMassSpec() {
        const { ctx, canvas } = this;
        const width = canvas.width;
        const height = canvas.height;
        const margin = 60;
        const plotWidth = width - 2 * margin;
        const plotHeight = height - 2 * margin;
        
        // Background
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, width, height);
        
        // Grid
        ctx.strokeStyle = '#f1f5f9';
        ctx.lineWidth = 1;
        for (let i = 0; i <= 10; i++) {
            const x = margin + (i * plotWidth / 10);
            const y = margin + (i * plotHeight / 10);
            
            ctx.beginPath();
            ctx.moveTo(x, margin);
            ctx.lineTo(x, margin + plotHeight);
            ctx.stroke();
            
            ctx.beginPath();
            ctx.moveTo(margin, y);
            ctx.lineTo(margin + plotWidth, y);
            ctx.stroke();
        }
        
        // Data points
        if (this.data.length >= 2 && this.data.length % 2 === 0) {
            const points = [];
            for (let i = 0; i < this.data.length; i += 2) {
                const mz = this.data[i];
                const intensity = this.data[i + 1];
                points.push({ mz, intensity });
            }
            
            // Find ranges
            const mzValues = points.map(p => p.mz);
            const intensities = points.map(p => p.intensity);
            
            const mzMin = Math.min(...mzValues);
            const mzMax = Math.max(...mzValues);
            const intensityMax = Math.max(...intensities);
            
            // Draw peaks
            points.forEach(point => {
                // Map high m/z to left, low m/z to right
                const x = margin + ((mzMax - point.mz) / (mzMax - mzMin)) * plotWidth;
                const barHeight = (point.intensity / intensityMax) * plotHeight;
                
                ctx.fillStyle = '#8b5cf6';
                ctx.fillRect(x - 3, margin + plotHeight - barHeight, 6, barHeight);
            });
            
            // Draw axis numbers
            ctx.fillStyle = '#374151';
            ctx.font = '12px -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif';
            ctx.textAlign = 'center';
            
            // X-axis numbers (m/z) - high to low (right to left)
            for (let i = 0; i <= 5; i++) {
                const x = margin + (i * plotWidth / 5);
                const value = mzMax - (i * (mzMax - mzMin) / 5);
                ctx.fillText(value.toFixed(0), x, height - 15);
            }
            
            // Y-axis numbers (intensity)
            ctx.textAlign = 'right';
            for (let i = 0; i <= 5; i++) {
                const y = margin + (i * plotHeight / 5);
                const value = (i * intensityMax / 5);
                ctx.fillText(value.toFixed(0), margin - 10, y + 4);
            }
        }
        
        // Axis labels
        ctx.fillStyle = '#374151';
        ctx.font = '14px -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('m/z', width/2, height - 5);
        
        ctx.save();
        ctx.translate(15, height/2);
        ctx.rotate(-Math.PI/2);
        ctx.fillText('Intensity', 0, 0);
        ctx.restore();
    }

    setHsqcSwapped(swapped) {
        this.isHsqcSwapped = swapped;
        if (this.type === 'hsqc') {
            this.draw();
        }
    }
}

// Global visualizers
let visualizers = {};

// Initialize visualizers
function initializeVisualizers() {
    visualizers.hsqc = new SpectrumVisualizer('hsqc-canvas', 'hsqc');
    visualizers.h_nmr = new SpectrumVisualizer('h_nmr-canvas', 'h_nmr');
    visualizers.c_nmr = new SpectrumVisualizer('c_nmr-canvas', 'c_nmr');
    visualizers.mass_spec = new SpectrumVisualizer('mass_spec-canvas', 'mass_spec');
}

// Update visualizations when data changes
function updateVisualizations() {
    // HSQC
    const hsqcData = collectTableData('hsqc');
    if (hsqcData.length > 0) {
        visualizers.hsqc.updateData(hsqcData);
    } else {
        visualizers.hsqc.updateData([]);
    }
    
    // H NMR
    const hNmrData = collectTableData('h_nmr');
    if (hNmrData.length > 0) {
        visualizers.h_nmr.updateData(hNmrData);
    } else {
        visualizers.h_nmr.updateData([]);
    }
    
    // C NMR
    const cNmrData = collectTableData('c_nmr');
    if (cNmrData.length > 0) {
        visualizers.c_nmr.updateData(cNmrData);
    } else {
        visualizers.c_nmr.updateData([]);
    }
    
    // Mass Spec
    const massSpecData = collectTableData('mass_spec');
    if (massSpecData.length > 0) {
        visualizers.mass_spec.updateData(massSpecData);
    } else {
        visualizers.mass_spec.updateData([]);
    }
}

// Collect data from table
function collectTableData(tableId) {
    const table = document.getElementById(`${tableId}-table`);
    const rows = table.querySelectorAll('tbody tr');
    const data = [];
    
    rows.forEach(row => {
        const inputs = row.querySelectorAll('input');
        inputs.forEach(input => {
            const value = parseFloat(input.value);
            if (!isNaN(value) && input.value.trim()) {
                data.push(value);
            }
        });
    });
    
    return data;
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
    
    // Update visualization
    visualizers.hsqc.setHsqcSwapped(!visualizers.hsqc.isHsqcSwapped);
    
    // Update visualization data
    updateVisualizations();
}

// Initialize when DOM is loaded
document.addEventListener('DOMContentLoaded', function() {
    initializeVisualizers();
    
    // Update visualizations when inputs change
    document.addEventListener('input', function(e) {
        if (e.target.classList.contains('spreadsheet-input')) {
            updateVisualizations();
        }
    });
});
