import SpreadsheetTable from '../spreadsheet/SpreadsheetTable'
import { useMainPageStore } from '../../store/store'

function HNMRInput() {
  const { h_nmr, setHNMR } = useMainPageStore()
  
  const columns = ['H Shift (ppm)']
  
  return (
    <div className="hnmr-input">
      <h3>H NMR Data</h3>
      <p>Enter H chemical shift values in ppm</p>
      <SpreadsheetTable
        columns={columns}
        data={h_nmr}
        onChange={setHNMR}
        rowsPerValue={1}
        maxColumns={1}
      />
    </div>
  )
}

export default HNMRInput

