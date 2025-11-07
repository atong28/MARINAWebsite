/**
 * Unified spreadsheet table component for all spectral input types.
 * Combines HSQC, H NMR, C NMR, and Mass Spec into one spreadsheet with separate columns.
 */
import { useMemo, useCallback, useRef } from 'react'
import { HotTable, type HotTableRef } from '@handsontable/react-wrapper'
import { registerAllModules } from 'handsontable/registry'
import 'handsontable/styles/handsontable.min.css';
import 'handsontable/styles/ht-theme-main.min.css';
import './SpreadsheetTable.css'

registerAllModules()

interface UnifiedSpreadsheetTableProps {
  hsqc: number[]
  h_nmr: number[]
  c_nmr: number[]
  mass_spec: number[]
  onHSQCChange: (data: number[]) => void
  onHNMRChange: (data: number[]) => void
  onCNMRChange: (data: number[]) => void
  onMassSpecChange: (data: number[]) => void
}

function UnifiedSpreadsheetTable({
  hsqc,
  h_nmr,
  c_nmr,
  mass_spec,
  onHSQCChange,
  onHNMRChange,
  onCNMRChange,
  onMassSpecChange,
}: UnifiedSpreadsheetTableProps) {
  const hotTableRef = useRef<HotTableRef>(null)

  // Column definitions: HSQC (3), H NMR (1), C NMR (1), Mass Spec (2) = 7 total
  const columns = useMemo(() => [
    'HSQC H Shift (ppm)',
    'HSQC C Shift (ppm)',
    'HSQC Intensity',
    'H NMR Shift (ppm)',
    'C NMR Shift (ppm)',
    'Mass Spec m/z',
    'Mass Spec Intensity',
  ], [])

  // Fixed maximum rows: 400 rows with scrolling, always show 10 rows
  const maxRows = 400

  // Convert all data types to 2D array format - always 400 rows
  const tableData = useMemo(() => {
    const result: (number | string)[][] = []

    for (let i = 0; i < maxRows; i++) {
      const row: (number | string)[] = []

      // HSQC columns (0, 1, 2)
      const hsqcIdx = i * 3
      row.push(hsqcIdx < hsqc.length ? hsqc[hsqcIdx] ?? '' : '')
      row.push(hsqcIdx + 1 < hsqc.length ? hsqc[hsqcIdx + 1] ?? '' : '')
      row.push(hsqcIdx + 2 < hsqc.length ? hsqc[hsqcIdx + 2] ?? '' : '')

      // H NMR column (3)
      row.push(i < h_nmr.length ? h_nmr[i] ?? '' : '')

      // C NMR column (4)
      row.push(i < c_nmr.length ? c_nmr[i] ?? '' : '')

      // Mass Spec columns (5, 6)
      const massSpecIdx = i * 2
      row.push(massSpecIdx < mass_spec.length ? mass_spec[massSpecIdx] ?? '' : '')
      row.push(massSpecIdx + 1 < mass_spec.length ? mass_spec[massSpecIdx + 1] ?? '' : '')

      result.push(row)
    }

    return result
  }, [hsqc, h_nmr, c_nmr, mass_spec])

  // Create column definitions
  const columnDefs = useMemo(() => {
    return columns.map((header, index) => ({
      data: index,
      type: 'numeric' as const,
      editor: 'numeric' as const,
      readOnly: false,
      title: header,
      width: 150,
      allowEmpty: true,
      numericFormat: {
        pattern: '0.00',
        culture: 'en-US'
      }
    }))
  }, [columns])

  // Calculate width from column definitions (row header + all columns)
  const calculatedWidth = useMemo(() => {
    const rowHeaderWidth = 50
    const totalColumnWidth = columnDefs.reduce((sum, col) => sum + (col.width || 150), 0)
    return rowHeaderWidth + totalColumnWidth
  }, [columnDefs])

  // Fixed height for 10 visible rows with scrolling
  const calculatedHeight = 30 + (10 * 24) // header (30px) + 10 rows (24px each) = 270px

  // Convert 2D array back to separate data arrays on changes
  const handleAfterChange = useCallback((changes: any[] | null, source: string) => {
    if (!changes || source === 'loadData' || !hotTableRef.current?.hotInstance) {
      return
    }

    const currentData = hotTableRef.current.hotInstance.getData() as (number | string)[][]

    // Extract HSQC data (columns 0, 1, 2)
    const hsqcData: number[] = []
    currentData.forEach((row) => {
      for (let j = 0; j < 3; j++) {
        const val = row[j]
        if (val !== '' && val !== null && val !== undefined) {
          const num = typeof val === 'number' ? val : parseFloat(String(val))
          if (!isNaN(num)) {
            hsqcData.push(num)
          }
        }
      }
    })
    onHSQCChange(hsqcData)

    // Extract H NMR data (column 3)
    const hNmrData: number[] = []
    currentData.forEach((row) => {
      const val = row[3]
      if (val !== '' && val !== null && val !== undefined) {
        const num = typeof val === 'number' ? val : parseFloat(String(val))
        if (!isNaN(num)) {
          hNmrData.push(num)
        }
      }
    })
    onHNMRChange(hNmrData)

    // Extract C NMR data (column 4)
    const cNmrData: number[] = []
    currentData.forEach((row) => {
      const val = row[4]
      if (val !== '' && val !== null && val !== undefined) {
        const num = typeof val === 'number' ? val : parseFloat(String(val))
        if (!isNaN(num)) {
          cNmrData.push(num)
        }
      }
    })
    onCNMRChange(cNmrData)

    // Extract Mass Spec data (columns 5, 6)
    const massSpecData: number[] = []
    currentData.forEach((row) => {
      for (let j = 5; j < 7; j++) {
        const val = row[j]
        if (val !== '' && val !== null && val !== undefined) {
          const num = typeof val === 'number' ? val : parseFloat(String(val))
          if (!isNaN(num)) {
            massSpecData.push(num)
          }
        }
      }
    })
    onMassSpecChange(massSpecData)
  }, [onHSQCChange, onHNMRChange, onCNMRChange, onMassSpecChange])

  return (
    <HotTable
      ref={hotTableRef}
      data={tableData}
      columns={columnDefs}
      colHeaders={columns}
      rowHeaders={true}
      height={calculatedHeight}
      width={calculatedWidth}
      afterChange={handleAfterChange}
      licenseKey="non-commercial-and-evaluation"
      themeName="ht-theme-main"
    />
  )
}

export default UnifiedSpreadsheetTable

