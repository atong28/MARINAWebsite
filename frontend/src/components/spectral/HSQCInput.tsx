import SpreadsheetTable from '../spreadsheet/SpreadsheetTable'
import { useMainPageStore } from '../../store/store'

function HSQCInput() {
  const { hsqc, setHSQC } = useMainPageStore()
  
  // HSQC has 3 columns: H shift, C shift, intensity
  const columns = ['H Shift (ppm)', 'C Shift (ppm)', 'Intensity']
  
  return (
    <div className="hsqc-input">
      <h3>HSQC NMR Data</h3>
      <p>Enter H shift, C shift, and intensity values (3 columns)</p>
      <SpreadsheetTable
        columns={columns}
        data={hsqc}
        onChange={setHSQC}
        rowsPerValue={3}
        maxColumns={3}
      />
    </div>
  )
}

export default HSQCInput

