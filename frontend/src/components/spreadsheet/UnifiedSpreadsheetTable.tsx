/**
 * Unified spreadsheet table component for all spectral input types.
 * Combines HSQC, H NMR, C NMR, and Mass Spec into one spreadsheet with separate columns.
 */
import { useMemo, useCallback, useRef, useState, useEffect } from 'react'

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
  onValidationChange?: (summary: { hsqcInvalid: number; hInvalid: number; cInvalid: number; msInvalid: number; anyInvalid: boolean }) => void
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
  onValidationChange,
}: UnifiedSpreadsheetTableProps) {
  const hotTableRef = useRef<HotTableRef>(null)
  const [validationSummary, setValidationSummary] = useState<{ hsqcInvalid: number; hInvalid: number; cInvalid: number; msInvalid: number }>({ hsqcInvalid: 0, hInvalid: 0, cInvalid: 0, msInvalid: 0 })
  
  // Track when we're doing an internal update (from our own extractAndEmit)
  // This prevents the feedback loop where our update triggers a re-render which reloads Handsontable
  const isInternalUpdateRef = useRef(false)
  
  // Track the last data we sent to Handsontable to detect external changes
  const lastSyncedDataRef = useRef<(number | string)[][] | null>(null)

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

  // Track previous props to detect changes
  const prevPropsRef = useRef<{ hsqc: number[]; h_nmr: number[]; c_nmr: number[]; mass_spec: number[] }>({
    hsqc: [],
    h_nmr: [],
    c_nmr: [],
    mass_spec: [],
  })

  // Convert all data types to 2D array format - always 400 rows
  const tableData = useMemo(() => {
    // Update prev props
    prevPropsRef.current = { hsqc, h_nmr, c_nmr, mass_spec }

    const result: (number | string)[][] = []

    const getNum = (arr: number[], idx: number): number | '' => {
      if (!arr || idx < 0 || idx >= arr.length) return ''
      const v = arr[idx]
      return typeof v === 'number' && isFinite(v) ? v : ''
    }

    for (let i = 0; i < maxRows; i++) {
      const row: (number | string)[] = []

      // HSQC columns (0, 1, 2)
      const hsqcIdx = i * 3
      row.push(getNum(hsqc, hsqcIdx))
      row.push(getNum(hsqc, hsqcIdx + 1))
      row.push(getNum(hsqc, hsqcIdx + 2))

      // H NMR column (3)
      row.push(getNum(h_nmr, i))

      // C NMR column (4)
      row.push(getNum(c_nmr, i))

      // Mass Spec columns (5, 6)
      const massSpecIdx = i * 2
      row.push(getNum(mass_spec, massSpecIdx))
      row.push(getNum(mass_spec, massSpecIdx + 1))

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


  // Helpers: row validation per data type
  const parseNumber = (val: any): number | null => {
    if (val === '' || val === null || val === undefined) return null
    const num = typeof val === 'number' ? val : parseFloat(String(val))
    return isNaN(num) ? null : num
  }

  const isHSQCRowComplete = (row: (number | string)[]) => {
    const h = parseNumber(row[0]); const c = parseNumber(row[1]); const i = parseNumber(row[2])
    const filled = h !== null || c !== null || i !== null
    const allFilled = h !== null && c !== null && i !== null
    return !filled || allFilled
  }
  const isHNMRRowComplete = (row: (number | string)[]) => {
    const h = parseNumber(row[3])
    return h === null || h !== null
  }
  const isCNMRRowComplete = (row: (number | string)[]) => {
    const c = parseNumber(row[4])
    return c === null || c !== null
  }
  const isMSRowComplete = (row: (number | string)[]) => {
    const mz = parseNumber(row[5]); const inten = parseNumber(row[6])
    const filled = mz !== null || inten !== null
    const bothFilled = mz !== null && inten !== null
    return !filled || bothFilled
  }

  const extractAndEmit = useCallback((gridData: (number | string)[][]) => {
    if (!hotTableRef.current?.hotInstance) {
      return
    }

    const currentData = gridData
    
    // Mark that we're doing an internal update to prevent feedback loop
    isInternalUpdateRef.current = true
    
    // Update the last synced data to match what's in Handsontable
    lastSyncedDataRef.current = currentData.map(row => [...row])

    // Position-preserving emission: fixed-size arrays using NaN for blanks
    const hsqcData: number[] = new Array(maxRows * 3).fill(NaN)
    const hNmrData: number[] = new Array(maxRows).fill(NaN)
    const cNmrData: number[] = new Array(maxRows).fill(NaN)
    const massSpecData: number[] = new Array(maxRows * 2).fill(NaN)

    let hsqcInvalid = 0, hInvalid = 0, cInvalid = 0, msInvalid = 0

    currentData.forEach((row, rIdx) => {
      // HSQC: always persist per-cell values; validate separately
      const h = parseNumber(row[0]); const c = parseNumber(row[1]); const i = parseNumber(row[2])
      const baseHSQC = rIdx * 3
      hsqcData[baseHSQC] = h === null ? NaN : h
      hsqcData[baseHSQC + 1] = c === null ? NaN : c
      hsqcData[baseHSQC + 2] = i === null ? NaN : i
      if (!isHSQCRowComplete(row)) hsqcInvalid++

      // H NMR
      const hv = parseNumber(row[3])
      hNmrData[rIdx] = hv === null ? NaN : hv
      if (!isHNMRRowComplete(row)) hInvalid++

      // C NMR
      const cv = parseNumber(row[4])
      cNmrData[rIdx] = cv === null ? NaN : cv
      if (!isCNMRRowComplete(row)) cInvalid++

      // MS
      const mz = parseNumber(row[5]); const inten = parseNumber(row[6])
      const baseMS = rIdx * 2
      massSpecData[baseMS] = mz === null ? NaN : mz
      massSpecData[baseMS + 1] = inten === null ? NaN : inten
      if (!isMSRowComplete(row)) msInvalid++
    })

    setValidationSummary({ hsqcInvalid, hInvalid, cInvalid, msInvalid })
    if (typeof onValidationChange === 'function') {
      onValidationChange({ hsqcInvalid, hInvalid, cInvalid, msInvalid, anyInvalid: (hsqcInvalid + hInvalid + cInvalid + msInvalid) > 0 })
    }
    
    // Emit changes to store - this will trigger re-renders
    onHSQCChange(hsqcData)
    onHNMRChange(hNmrData)
    onCNMRChange(cNmrData)
    onMassSpecChange(massSpecData)
    
    // Clear the internal update flag in the next tick to allow re-render cycle to complete
    setTimeout(() => {
      isInternalUpdateRef.current = false
    }, 0)
  }, [maxRows, onHSQCChange, onHNMRChange, onCNMRChange, onMassSpecChange, onValidationChange])

  // Convert 2D array back to separate data arrays on changes
  const handleAfterChange = useCallback((changes: any[] | null, source: string) => {
    // Ignore changes from loadData source to prevent loops
    if (!changes || source === 'loadData' || !hotTableRef.current?.hotInstance) {
      return
    }

    const currentData = hotTableRef.current.hotInstance.getData() as (number | string)[][]
    extractAndEmit(currentData)
  }, [extractAndEmit])

  // Sync Handsontable data when props change externally (not from our own update)
  useEffect(() => {
    // Skip if we're in the middle of an internal update (our own extractAndEmit call)
    if (isInternalUpdateRef.current) {
      return
    }
    
    // Skip if Handsontable isn't ready
    if (!hotTableRef.current?.hotInstance) {
      return
    }
    
    const hot = hotTableRef.current.hotInstance
    
    // Compare current tableData with what we last synced
    // Only update if they're different (external change)
    const dataChanged = lastSyncedDataRef.current === null || 
      JSON.stringify(tableData) !== JSON.stringify(lastSyncedDataRef.current)
    
    if (dataChanged) {
      // Update Handsontable with the new data from props
      hot.loadData(tableData)
      // Update our tracking ref
      lastSyncedDataRef.current = tableData.map(row => [...row])
    }
  }, [tableData])

  const handleCondense = useCallback(() => {
    if (!hotTableRef.current?.hotInstance) return
    const hot = hotTableRef.current.hotInstance
    const currentData = (hot.getData() as (number | string | undefined)[][]).map((r) => r ?? [])
    const newData: (number | string)[][] = currentData.map((r) => r.map((v) => (v === undefined ? '' : v)) as (number | string)[])

    // Remove fully empty rows (per data type) by shifting up only over empty rows; partial rows stay
    const condenseSegment = (cols: number[]) => {
      let write = 0
      for (let read = 0; read < newData.length; read++) {
        const row = newData[read] || []
        const vals = cols.map((c) => parseNumber(row[c] ?? ''))
        const allEmpty = vals.every((v) => v === null)
        if (!allEmpty) {
          if (write !== read) {
            cols.forEach((c) => {
              const fromRow = newData[read] || []
              const toRow = newData[write] || []
              const val = (fromRow[c] === undefined ? '' : fromRow[c]) as string | number
              toRow[c] = val
              newData[write] = toRow
            })
            cols.forEach((c) => {
              const fromRow = newData[read] || []
              fromRow[c] = ''
              newData[read] = fromRow
            })
          }
          write++
        }
      }
    }

    // Condense per data type
    condenseSegment([0,1,2]) // HSQC
    condenseSegment([3])     // H NMR
    condenseSegment([4])     // C NMR
    condenseSegment([5,6])   // Mass Spec
    
    // For condense, we explicitly want to update both store and table
    // Mark as internal update to prevent the useEffect from interfering
    isInternalUpdateRef.current = true
    extractAndEmit(newData)
    hot.loadData(newData)
    lastSyncedDataRef.current = newData.map(row => [...row])
    setTimeout(() => {
      isInternalUpdateRef.current = false
    }, 0)
  }, [extractAndEmit])

  return (
    <div>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 8 }}>
        <div style={{ fontSize: 12, color: '#666' }}>
          HSQC invalid: {validationSummary.hsqcInvalid} · H invalid: {validationSummary.hInvalid} · C invalid: {validationSummary.cInvalid} · MS invalid: {validationSummary.msInvalid}
        </div>
        <button onClick={handleCondense} style={{ padding: '6px 10px' }}>Condense rows</button>
      </div>
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
        copyPaste={{ pasteMode: 'overwrite' } as any}
        fillHandle={false}
    />
    </div>
  )
}

export default UnifiedSpreadsheetTable

