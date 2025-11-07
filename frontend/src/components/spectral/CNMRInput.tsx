import SpreadsheetTable from '../spreadsheet/SpreadsheetTable'
import { useMainPageStore } from '../../store/store'

function CNMRInput() {
  const { c_nmr, setCNMR } = useMainPageStore()
  
  const columns = ['C Shift (ppm)']
  
  return (
    <div className="cnmr-input">
      <h3>C NMR Data</h3>
      <p>Enter C chemical shift values in ppm</p>
      <SpreadsheetTable
        columns={columns}
        data={c_nmr}
        onChange={setCNMR}
        rowsPerValue={1}
        maxColumns={1}
      />
    </div>
  )
}

export default CNMRInput

