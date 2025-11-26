import { useMainPageStore } from '../../store/store'
import UnifiedSpreadsheetTable from '../spreadsheet/UnifiedSpreadsheetTable'
import './SpectralInputTabs.css'
import RetrievalMWFilter from './RetrievalMWFilter'

interface SpectralInputTabsProps {
  onValidationChange?: (summary: { hsqcInvalid: number; hInvalid: number; cInvalid: number; msInvalid: number; anyInvalid: boolean }) => void
}

function SpectralInputTabs({ onValidationChange }: SpectralInputTabsProps) {
  const { 
    hsqc, 
    h_nmr, 
    c_nmr, 
    mass_spec, 
    mw, 
    retrievalMwMin,
    retrievalMwMax,
    setHSQC, 
    setHNMR, 
    setCNMR, 
    setMassSpec, 
    setMW,
    setRetrievalMwRange,
  } = useMainPageStore()
  
  return (
    <div className="spectral-input-tabs">
      <div className="unified-input-header">
        <h3>Spectral Data Input</h3>
        <p>Enter all spectral data in the spreadsheet below. Columns are organized by data type:</p>
        <ul>
          <li><strong>HSQC:</strong> H Shift, C Shift, Intensity (3 columns)</li>
          <li><strong>H NMR:</strong> H Shift (1 column)</li>
          <li><strong>C NMR:</strong> C Shift (1 column)</li>
          <li><strong>Mass Spec:</strong> m/z, Intensity (2 columns)</li>
        </ul>
      </div>
      
      <div className="tab-content">
        <UnifiedSpreadsheetTable
          hsqc={hsqc}
          h_nmr={h_nmr}
          c_nmr={c_nmr}
          mass_spec={mass_spec}
          onHSQCChange={setHSQC}
          onHNMRChange={setHNMR}
          onCNMRChange={setCNMR}
          onMassSpecChange={setMassSpec}
          onValidationChange={onValidationChange}
        />
      </div>
      
      <div className="mw-input">
        <label>
          Molecular Weight (g/mol):
          <input
            type="number"
            step="0.01"
            value={mw || ''}
            onChange={(e) => setMW(e.target.value ? parseFloat(e.target.value) : null)}
            placeholder="Enter molecular weight"
          />
        </label>
      </div>

      <RetrievalMWFilter
        minMw={retrievalMwMin}
        maxMw={retrievalMwMax}
        onChange={setRetrievalMwRange}
      />
    </div>
  )
}

export default SpectralInputTabs

