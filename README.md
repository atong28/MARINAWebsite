# MARINA: Molecular Structure Annotator

A deep learning system for predicting molecular structures from spectral data using PyTorch Lightning and modern web technologies.

## Overview

MARINA (Multimodal Annotation and Retrieval for Identification of NAtural products) is a machine learning model that predicts molecular structures from various types of spectral data including:

- **HSQC NMR**: 2D heteronuclear single quantum coherence spectroscopy
- **H NMR**: Proton nuclear magnetic resonance
- **C NMR**: Carbon-13 nuclear magnetic resonance  
- **Mass Spectrometry**: Mass-to-charge ratio and intensity data
- **Molecular Weight**: Molecular mass information

The system uses a transformer-based architecture with cross-attention mechanisms to process multi-modal spectral inputs and predict molecular fingerprints for structure retrieval.

## Features

### Modern Web Interface
- **Tabbed Input Interface**: Separate tabs for each spectral data type
- **Editable Tables**: Dynamic tables for entering spectral peaks with validation
- **Clipboard Support**: Paste data directly from spreadsheets (CSV/TSV)
- **Real-time Validation**: Input validation with helpful error messages
- **Responsive Design**: Works on desktop and mobile devices

### Data Input Formats
- **HSQC NMR**: C-13 shift, H-1 shift, and intensity values (3 columns)
- **H NMR**: H-1 chemical shift values (1 column)
- **C NMR**: C-13 chemical shift values (1 column)
- **Mass Spec**: m/z and intensity pairs (2 columns)
- **Molecular Weight**: Single scalar value in g/mol

### Prediction Results
- **Top-K Retrieval**: Configurable number of similar structures (5-50)
- **Similarity Scoring**: Tanimoto similarity between predicted and retrieved fingerprints
- **PubChem Integration**: Direct links to view structures on PubChem
- **SMILES Display**: Chemical structure notation for each result

## Quick Start

### Prerequisites

- Docker and Docker Compose
- At least 4GB RAM available
- Model checkpoint files in `data/` directory:
  - `best.ckpt`: PyTorch Lightning model checkpoint (~2GB)
  - `count_hashes_under_radius_10.pkl`: Dataset count for entropy fp generation
  - `params.json`: Model configuration parameters
  - `rankingset.pt`: Sparse tensor with molecular fingerprints
  - `rankingset_meta.pkl`: Metadata mapping indices to SMILES strings
  - `index.pkl`: Dataset index mapping (if using index-based prediction)

### Using Docker (Recommended)

1. **Build the Docker image:**
   ```bash
   docker build -t marina-predictor:latest .
   ```

2. **Run with Docker Compose:**
   ```bash
   PORT=[PORT] docker-compose up -d
   ```

3. **Quick restart with logs (recommended):**
   ```bash
   ./restart.sh
   ```

4. **Access the web interface:**
   Open your browser to `http://localhost:[PORT]`

### Logs

Runtime logs are written to the `logs/` directory in the project root (mounted into the container):

- `logs/app.log`: application logs (warnings and errors)
- `logs/gunicorn.access.log`: HTTP access logs
- `logs/gunicorn.error.log`: Gunicorn errors

Logs rotate automatically (app log ~2MB x 5 backups). Ensure the `logs/` directory exists or is mounted (already configured in `docker-compose.yml`).

### Using Pixi Package Manager

1. **Install Pixi:**
   ```bash
   curl -fsSL https://pixi.sh/install.sh | bash
   ```

2. **Install dependencies:**
   ```bash
   pixi install
   ```

3. **Run the server:**
   ```bash
   pixi run gunicorn -w 1 -b 0.0.0.0:5000 app:app
   ```

### Using pip

1. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Run the server:**
   ```bash
   gunicorn -w 1 -b 0.0.0.0:5000 app:app
   ```

## Usage Guide

### Web Interface

1. **Select Input Type**: Click on the appropriate tab for your spectral data
2. **Enter Data**: 
   - Use the "Add Row" button to add more data points
   - Paste data from spreadsheets using "Paste from Clipboard" or Ctrl+V
   - Use keyboard shortcuts: Arrow keys to navigate, Ctrl+C to copy, Ctrl+V to paste
   - Select multiple cells by clicking and dragging
   - Enter molecular weight in the dedicated field
3. **Configure Results**: Choose the number of results to retrieve (5-50)
4. **Run Prediction**: Click "Predict Structure" to start the analysis
5. **View Results**: Browse the ranked list of similar molecular structures

#### Spreadsheet Features
- **Keyboard Navigation**: Use arrow keys to move between cells
- **Range Selection**: Click and drag to select multiple cells
- **Copy/Paste**: Ctrl+C to copy, Ctrl+V to paste data
- **Clear Cells**: Delete or Backspace to clear selected cells
- **Select All**: Ctrl+A to select all cells

### API Endpoints

#### Health Check
```bash
curl http://localhost:5000/health
```

#### Predict from Raw Data
```bash
curl -X POST -H "Content-Type: application/json" \
  -d '{
    "raw": {
      "hsqc": [120.5, 7.2, 1.0, 110.3, 6.8, 0.8],
      "h_nmr": [7.2, 6.8, 3.5],
      "mw": 180.16
    },
    "k": 10
  }' \
  http://localhost:5000/predict
```

#### Get Metadata
```bash
curl http://localhost:5000/meta/123
```

## Data Format Specifications

### HSQC NMR Data
- **Format**: Array of triplets [H_shift, C_shift, intensity]
- **Example**: `[7.2, 120.5, 1.0, 6.8, 110.3, 0.8]`
- **Units**: Chemical shifts in ppm, intensities as relative values

### H NMR Data  
- **Format**: Array of chemical shift values
- **Example**: `[7.2, 6.8, 3.5, 2.1]`
- **Units**: Chemical shifts in ppm

### C NMR Data
- **Format**: Array of chemical shift values  
- **Example**: `[120.5, 110.3, 55.2, 21.7]`
- **Units**: Chemical shifts in ppm

### Mass Spectrometry Data
- **Format**: Array of pairs [m/z, intensity]
- **Example**: `[180.5, 1000, 150.3, 800, 120.1, 600]`
- **Units**: m/z in atomic mass units, intensities as relative values

### Molecular Weight
- **Format**: Single scalar value
- **Example**: `180.16`
- **Units**: g/mol

### Environment Variables

- `PORT`: Server port (default: 5000)
- `DATA_DIR`: Data directory path

## Citation

If you use this software in your research, please cite:

```

```
