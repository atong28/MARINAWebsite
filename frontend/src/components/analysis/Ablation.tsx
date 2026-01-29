import { useMemo, useCallback, useState, useRef, useEffect } from 'react'
import UnifiedSpreadsheetTable from '../spreadsheet/UnifiedSpreadsheetTable'
import { useAnalysisPageStore, useMainPageStore } from '../../store/store'
import { useAblation } from '../../services/api'
import SimilarityMap from './SimilarityMap'
import BitComparison from './BitComparison'
import MoleculeOverlay from './MoleculeOverlay'
import './Ablation.css'

// Import type from store
type SpectralDataSnapshot = {
  hsqc: number[]
  h_nmr: number[]
  c_nmr: number[]
  mass_spec: number[]
  mw: number | null
}

const DEFAULT_THRESHOLD = 0.5

// Helper function to filter out NaN, null, undefined, and non-finite values
function filterNaNValues(arr: number[]): number[] {
  return arr.filter((val) => {
    return val !== null && val !== undefined && !Number.isNaN(val) && isFinite(val)
  })
}

function Ablation() {
  const {
    selectedMolecule,
    setAblationSpectralData,
    ablationPredictedFp,
    setAblationPredictedFp,
    ablationSimilarityMap,
    setAblationSimilarityMap,
    originalPredictedFp,
    originalSimilarityMap,
    originalSpectralData,
  } = useAnalysisPageStore()
  const selectedModelId = useMainPageStore((s) => s.selectedModelId)

  // Ensure we always use the latest store state for props
  // Using a selector ensures we always get the latest state and component re-renders on changes
  const propsToUse = useAnalysisPageStore((state) => state.ablationSpectralData)

  const [changeOverlaySvg, setChangeOverlaySvg] = useState<string | null>(null)
  const [validationSummary, setValidationSummary] = useState<{ anyInvalid: boolean }>({ anyInvalid: false })
  const hasResults = Boolean(ablationPredictedFp)
  
  // Use a ref to track the latest propsToUse to avoid stale closures
  const ablationSpectralDataRef = useRef(propsToUse)
  useEffect(() => {
    ablationSpectralDataRef.current = propsToUse
  }, [propsToUse])

  const ablationMutation = useAblation({
    onSuccess: (data) => {
      setAblationPredictedFp(data.pred_fp)
      setAblationSimilarityMap(data.similarity_map ?? null)
      setChangeOverlaySvg(data.change_overlay_svg ?? null)
    },
  })

  const handleSpectralUpdate = useCallback(
    (partial: Partial<SpectralDataSnapshot>) => {
      // Read the absolute latest store state right before updating to avoid race conditions
      const absoluteLatestStoreState = useAnalysisPageStore.getState().ablationSpectralData
      
      // Merge with absolute latest store state to ensure we don't lose any concurrent updates
      const newState = {
        ...absoluteLatestStoreState,
        ...partial,
      }
      
      setAblationSpectralData(newState)
      
      // Update the ref immediately so next callback has latest state
      ablationSpectralDataRef.current = newState
    },
    [setAblationSpectralData]
  )

  const handleResetToOriginal = () => {
    if (!originalSpectralData) return
    setAblationSpectralData(originalSpectralData)
    setAblationPredictedFp(null)
    setAblationSimilarityMap(null)
    setChangeOverlaySvg(null)
  }

  const handleRunAblation = () => {
    if (!selectedMolecule) return

    // Sanitize arrays by filtering out NaN, null, undefined, and non-finite values
    const sanitizedHSQC = filterNaNValues(propsToUse.hsqc)
    const sanitizedHNMR = filterNaNValues(propsToUse.h_nmr)
    const sanitizedCNMR = filterNaNValues(propsToUse.c_nmr)
    const sanitizedMassSpec = filterNaNValues(propsToUse.mass_spec)

    // Build request payload - only include arrays with at least one valid value
    const raw: any = {}
    if (sanitizedHSQC.length > 0) raw.hsqc = sanitizedHSQC
    if (sanitizedHNMR.length > 0) raw.h_nmr = sanitizedHNMR
    if (sanitizedCNMR.length > 0) raw.c_nmr = sanitizedCNMR
    if (sanitizedMassSpec.length > 0) raw.mass_spec = sanitizedMassSpec
    if (propsToUse.mw !== null && propsToUse.mw !== undefined && !Number.isNaN(propsToUse.mw) && isFinite(propsToUse.mw)) {
      raw.mw = propsToUse.mw
    }

    const payload = {
      raw,
      smiles: selectedMolecule.smiles,
      bit_threshold: DEFAULT_THRESHOLD,
      reference_fp: originalPredictedFp ?? undefined,
      model_id: selectedModelId ?? undefined,
    }

    ablationMutation.mutate(payload)
  }

  const bitDiff = useMemo(() => {
    const gained: number[] = []
    const lost: number[] = []

    if (!originalPredictedFp || !ablationPredictedFp) {
      return { gained, lost, unchanged: [] as number[] }
    }

    const originalBits = new Set<number>()
    const newBits = new Set<number>()

    originalPredictedFp.forEach((value, idx) => {
      if (value !== undefined && value > DEFAULT_THRESHOLD) {
        originalBits.add(idx)
      }
    })

    ablationPredictedFp.forEach((value, idx) => {
      if (value !== undefined && value > DEFAULT_THRESHOLD) {
        newBits.add(idx)
      }
    })

    const unchanged: number[] = []

    originalBits.forEach((bit) => {
      if (newBits.has(bit)) {
        unchanged.push(bit)
      } else {
        lost.push(bit)
      }
    })

    newBits.forEach((bit) => {
      if (!originalBits.has(bit)) {
        gained.push(bit)
      }
    })

    gained.sort((a, b) => a - b)
    lost.sort((a, b) => a - b)
    unchanged.sort((a, b) => a - b)

    return { gained, lost, unchanged }
  }, [originalPredictedFp, ablationPredictedFp])

  const changeOverlayBits = useMemo(() => {
    if (!hasResults) {
      return new Set<number>()
    }
    if (!bitDiff.gained.length && !bitDiff.lost.length) {
      return new Set<number>()
    }
    return new Set<number>([...bitDiff.gained, ...bitDiff.lost])
  }, [hasResults, bitDiff.gained, bitDiff.lost])

  const mwDisplay = propsToUse.mw ?? ''

  return (
    <section className="ablation-section">
      <div className="ablation-header">
        <div>
          <h3>Ablation Analysis</h3>
          <p>Try removing peaks to see how the fingerprint prediction changes.</p>
        </div>
        <div className="ablation-actions">
          <button
            type="button"
            onClick={handleResetToOriginal}
            disabled={!originalSpectralData}
          >
            Reset to Original Inputs
          </button>
        </div>
      </div>

      <div className="ablation-inputs">
        <UnifiedSpreadsheetTable
          hsqc={propsToUse.hsqc}
          h_nmr={propsToUse.h_nmr}
          c_nmr={propsToUse.c_nmr}
          mass_spec={propsToUse.mass_spec}
          onHSQCChange={(data) => handleSpectralUpdate({ hsqc: data })}
          onHNMRChange={(data) => handleSpectralUpdate({ h_nmr: data })}
          onCNMRChange={(data) => handleSpectralUpdate({ c_nmr: data })}
          onMassSpecChange={(data) => handleSpectralUpdate({ mass_spec: data })}
          onValidationChange={(summary) => setValidationSummary({ anyInvalid: summary.anyInvalid })}
        />
        <div className="ablation-mw-input">
          <label>
            Molecular Weight (Da)
            <input
              type="number"
              step="0.01"
              value={mwDisplay}
              onChange={(event) => {
                const value = event.target.value
                handleSpectralUpdate({ mw: value ? parseFloat(value) : null })
              }}
              placeholder="Enter molecular weight"
            />
          </label>
        </div>
      </div>

      <div className="ablation-predict-section">
        {validationSummary.anyInvalid && (
          <div className="validation-warning">Incomplete data</div>
        )}
        <button
          type="button"
          onClick={handleRunAblation}
          disabled={!selectedMolecule || ablationMutation.isPending || validationSummary.anyInvalid}
        >
          {ablationMutation.isPending ? 'Running…' : 'Run Prediction'}
        </button>
      </div>

      {ablationMutation.isError && (
        <div className="ablation-error">
          {ablationMutation.error?.message || 'Ablation prediction failed'}
        </div>
      )}

      {hasResults ? (
        <div className="ablation-results">
          <BitComparison
            originalFp={originalPredictedFp}
            newFp={ablationPredictedFp}
            threshold={DEFAULT_THRESHOLD}
          />

          <div className="similarity-map-grid">
            <SimilarityMap
              title="Original Prediction"
              image={originalSimilarityMap}
              description="Similarity map using the original spectral inputs"
            />
            <SimilarityMap
              title="Ablation Prediction"
              image={ablationSimilarityMap}
              isLoading={ablationMutation.isPending}
              description="Similarity map after applying ablation edits"
            />
          </div>

          <div className="ablation-overlays">
            <h4>Bit Change Overlays</h4>
            {changeOverlaySvg && selectedMolecule ? (
              <MoleculeOverlay
                smiles={selectedMolecule.smiles}
                svg={changeOverlaySvg}
                selectedBits={changeOverlayBits}
              />
            ) : (
              <div className="change-overlay-placeholder">No change detected.</div>
            )}
          </div>
        </div>
      ) : (
        <div className="ablation-results-placeholder">
          Run an ablation prediction to visualize how the fingerprint changes.
        </div>
      )}
    </section>
  )
}

export default Ablation

