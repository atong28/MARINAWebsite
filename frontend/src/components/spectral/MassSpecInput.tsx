import SpreadsheetTable from '../spreadsheet/SpreadsheetTable'
import { useMainPageStore } from '../../store/store'

function MassSpecInput() {
  const { mass_spec, setMassSpec } = useMainPageStore()
  
  const columns = ['m/z', 'Intensity']
  
  return (
    <div className="mass-spec-input">
      <h3>Mass Spectrometry Data</h3>
      <p>Enter m/z and intensity values (2 columns)</p>
      <SpreadsheetTable
        columns={columns}
        data={mass_spec}
        onChange={setMassSpec}
        rowsPerValue={2}
        maxColumns={2}
      />
    </div>
  )
}

export default MassSpecInput

