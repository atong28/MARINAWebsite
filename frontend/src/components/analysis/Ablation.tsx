import { useMemo, useCallback, useState } from 'react'
import UnifiedSpreadsheetTable from '../spreadsheet/UnifiedSpreadsheetTable'
import { useAnalysisPageStore } from '../../store/store'
import { useAblation } from '../../services/api'
import SimilarityMap from './SimilarityMap'
import BitComparison from './BitComparison'
import MoleculeOverlay from './MoleculeOverlay'
import './Ablation.css'

const DEFAULT_THRESHOLD = 0.5

function Ablation() {
  const {
    selectedMolecule,
    ablationSpectralData,
    setAblationSpectralData,
    ablationPredictedFp,
    setAblationPredictedFp,
    ablationSimilarityMap,
    setAblationSimilarityMap,
    originalPredictedFp,
    originalSimilarityMap,
    originalSpectralData,
  } = useAnalysisPageStore()

  const [changeOverlaySvg, setChangeOverlaySvg] = useState<string | null>(null)
  const hasResults = Boolean(ablationPredictedFp)

  const ablationMutation = useAblation({
    onSuccess: (data) => {
      setAblationPredictedFp(data.pred_fp)
      setAblationSimilarityMap(data.similarity_map ?? null)
      setChangeOverlaySvg(data.change_overlay_svg ?? null)
    },
  })

  const handleSpectralUpdate = useCallback(
    (partial: Partial<typeof ablationSpectralData>) => {
      setAblationSpectralData({
        ...ablationSpectralData,
        ...partial,
      })
    },
    [ablationSpectralData, setAblationSpectralData]
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

    const payload = {
      raw: {
        hsqc: ablationSpectralData.hsqc.length > 0 ? ablationSpectralData.hsqc : undefined,
        h_nmr: ablationSpectralData.h_nmr.length > 0 ? ablationSpectralData.h_nmr : undefined,
        c_nmr: ablationSpectralData.c_nmr.length > 0 ? ablationSpectralData.c_nmr : undefined,
        mass_spec: ablationSpectralData.mass_spec.length > 0 ? ablationSpectralData.mass_spec : undefined,
        mw: ablationSpectralData.mw ?? undefined,
      },
      smiles: selectedMolecule.smiles,
      bit_threshold: DEFAULT_THRESHOLD,
      reference_fp: originalPredictedFp ?? undefined,
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

  const mwDisplay = ablationSpectralData.mw ?? ''

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
          <button
            type="button"
            onClick={handleRunAblation}
            disabled={!selectedMolecule || ablationMutation.isPending}
          >
            {ablationMutation.isPending ? 'Running…' : 'Run Prediction'}
          </button>
        </div>
      </div>

      <div className="ablation-inputs">
        <UnifiedSpreadsheetTable
          hsqc={ablationSpectralData.hsqc}
          h_nmr={ablationSpectralData.h_nmr}
          c_nmr={ablationSpectralData.c_nmr}
          mass_spec={ablationSpectralData.mass_spec}
          onHSQCChange={(data) => handleSpectralUpdate({ hsqc: data })}
          onHNMRChange={(data) => handleSpectralUpdate({ h_nmr: data })}
          onCNMRChange={(data) => handleSpectralUpdate({ c_nmr: data })}
          onMassSpecChange={(data) => handleSpectralUpdate({ mass_spec: data })}
        />
        <div className="ablation-mw-input">
          <label>
            Molecular Weight (g/mol)
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

