import { useEffect, useRef } from 'react'
import { useNavigate, useLocation } from 'react-router-dom'
import { useAnalysisPageStore, useMainPageStore } from '../store/store'
import { api, useAnalyze } from '../services/api'
import FingerprintIndices from '../components/analysis/FingerprintIndices'
import MoleculeOverlay from '../components/analysis/MoleculeOverlay'
import SecondaryRetrieval from '../components/analysis/SecondaryRetrieval'
import Ablation from '../components/analysis/Ablation'
import './AnalysisPage.css'

const ABALATION_BIT_THRESHOLD = 0.5

function AnalysisPage() {
  const navigate = useNavigate()
  const location = useLocation()
  const { result } = location.state || {}
  
  const {
    selectedMolecule,
    selectedBits,
    retrievedFpIndices,
    bitEnvironments,
    moleculeSvgWithOverlays,
    originalSpectralData,
    originalSimilarityMap,
    setSelectedMolecule,
    setRetrievedFpIndices,
    setBitEnvironments,
    setMoleculeSvgWithOverlays,
    setOriginalSpectralData,
    setOriginalPredictedFp,
    setOriginalSimilarityMap,
    setAblationSpectralData,
    setAblationPredictedFp,
    setAblationSimilarityMap,
    ablationSeedKey,
    setAblationSeedKey,
    clearAnalysis,
  } = useAnalysisPageStore()
  
  const {
    hsqc,
    h_nmr,
    c_nmr,
    mass_spec,
    mw,
    predictedFp,
  } = useMainPageStore()
  
  // Track which result we've already processed to prevent duplicate calls
  const processedResultRef = useRef<string | null>(null)
  
  const analyzeMutation = useAnalyze({
    onSuccess: (data) => {
        setMoleculeSvgWithOverlays(data.molecule_svg)
        setRetrievedFpIndices(data.retrieved_molecule_fp_indices)
    },
    onError: (error) => {
      console.error('[AnalysisPage] Analyze endpoint error:', {
        message: error.message,
        error: error,
        stack: error.stack,
      })
    },
  })
  
  useEffect(() => {
    if (!selectedMolecule) {
      return
    }

    const seedKey = `${selectedMolecule.smiles}-${selectedMolecule.index}`

    if (ablationSeedKey === seedKey) {
      return
    }

    const snapshot = {
      hsqc: Array.from(hsqc),
      h_nmr: Array.from(h_nmr),
      c_nmr: Array.from(c_nmr),
      mass_spec: Array.from(mass_spec),
      mw: mw ?? null,
    }

    setOriginalSpectralData(snapshot)
    setAblationSpectralData(snapshot)
    setOriginalPredictedFp(predictedFp ? Array.from(predictedFp) : null)
    setOriginalSimilarityMap(null)
    setAblationPredictedFp(null)
    setAblationSimilarityMap(null)
    setAblationSeedKey(seedKey)
  }, [
    selectedMolecule,
    ablationSeedKey,
    hsqc,
    h_nmr,
    c_nmr,
    mass_spec,
    mw,
    predictedFp,
    setOriginalSpectralData,
    setAblationSpectralData,
    setOriginalPredictedFp,
    setOriginalSimilarityMap,
    setAblationPredictedFp,
    setAblationSimilarityMap,
    setAblationSeedKey,
  ])

  useEffect(() => {
    let isMounted = true

    const fetchOriginalSimilarity = async () => {
      if (!selectedMolecule || !originalSpectralData || !predictedFp) {
        return
      }

      const seedKey = `${selectedMolecule.smiles}-${selectedMolecule.index}`
      if (ablationSeedKey !== seedKey) {
        return
      }

      if (originalSimilarityMap !== null) {
        return
      }

      try {
        const response = await api.ablation({
          raw: {
            hsqc: originalSpectralData.hsqc.length ? originalSpectralData.hsqc : undefined,
            h_nmr: originalSpectralData.h_nmr.length ? originalSpectralData.h_nmr : undefined,
            c_nmr: originalSpectralData.c_nmr.length ? originalSpectralData.c_nmr : undefined,
            mass_spec: originalSpectralData.mass_spec.length ? originalSpectralData.mass_spec : undefined,
            mw: originalSpectralData.mw ?? undefined,
          },
          smiles: selectedMolecule.smiles,
          bit_threshold: ABALATION_BIT_THRESHOLD,
          reference_fp: predictedFp,
        })

        if (isMounted) {
          setOriginalSimilarityMap(response.similarity_map ?? null)
        }
      } catch (error) {
        if (isMounted) {
          console.error('[AnalysisPage] Failed to fetch original similarity map:', error)
          setOriginalSimilarityMap(null)
        }
      }
    }

    fetchOriginalSimilarity()

    return () => {
      isMounted = false
    }
  }, [
    selectedMolecule,
    originalSpectralData,
    predictedFp,
    originalSimilarityMap,
    ablationSeedKey,
    setOriginalSimilarityMap,
  ])

  useEffect(() => {
    if (result) {
      // Create a unique key for this result to track if we've processed it
      const resultKey = `${result.smiles}-${result.index}`
      
      // Skip if we've already processed this result
      if (processedResultRef.current === resultKey) {
        return
      }
      
      // Mark this result as processed
      processedResultRef.current = resultKey
      
      setSelectedMolecule(result)
      
      // Extract fingerprint indices from result
      if (result.retrieved_molecule_fp_indices) {
        setRetrievedFpIndices(result.retrieved_molecule_fp_indices)
      }
      if (result.bit_environments) {
        setBitEnvironments(result.bit_environments)
      }
      
      // Call analyze endpoint to get SVG with overlays
      // We need to construct retrieved_fp from indices (create sparse vector)
      // Only construct if we have valid indices
      let retrievedFp: number[] | undefined = undefined
      if (result.retrieved_molecule_fp_indices && result.retrieved_molecule_fp_indices.length > 0) {
        try {
          // Create a sparse vector with 1.0 at active indices
          const maxIndex = Math.max(...result.retrieved_molecule_fp_indices)
          if (maxIndex >= 0) {
            retrievedFp = new Array(maxIndex + 1).fill(0)
            result.retrieved_molecule_fp_indices.forEach((idx: number) => {
              if (idx >= 0 && idx <= maxIndex) {
                retrievedFp![idx] = 1.0
              }
            })
          }
        } catch (error) {
          console.error('[AnalysisPage] Error constructing retrieved_fp:', error)
          retrievedFp = undefined
        }
      }
      
      if (!retrievedFp || retrievedFp.length === 0) {
        console.warn('[AnalysisPage] No retrieved_fp available, cannot call analyze endpoint')
        return
      }
      
      // Use mutate directly - it's stable and doesn't need to be in dependencies
      analyzeMutation.mutate({
        smiles: result.smiles,
        retrieved_fp: retrievedFp,
      })
    }
    
    return () => {
      // Reset processed result ref when component unmounts or result changes
      processedResultRef.current = null
      clearAnalysis()
    }
  }, [result, setSelectedMolecule, setRetrievedFpIndices, setBitEnvironments, clearAnalysis])
  
  const handleBack = () => {
    navigate('/')
  }
  
  if (!selectedMolecule) {
    return (
      <div className="analysis-page">
        <button onClick={handleBack}>Back to Main</button>
        <p>No molecule selected for analysis.</p>
      </div>
    )
  }
  
  return (
    <div className="analysis-page">
      <header className="analysis-header">
        <button onClick={handleBack}>← Back to Main</button>
        <h2>Analysis: {selectedMolecule.name || selectedMolecule.smiles || `Molecule #${selectedMolecule.index}`}</h2>
      </header>
      
      <div className="molecule-info-section">
        <div className="smiles-display">
          <strong>SMILES:</strong> <code>{selectedMolecule.smiles}</code>
        </div>
        {selectedMolecule.database_links && (
          (selectedMolecule.database_links.coconut || selectedMolecule.database_links.lotus || selectedMolecule.database_links.npmrd) && (
            <div className="database-buttons">
              {selectedMolecule.database_links.coconut && (
                <a
                  href={selectedMolecule.database_links.coconut}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="database-btn coconut"
                >
                  View on COCONUT
                </a>
              )}
              {selectedMolecule.database_links.lotus && (
                <a
                  href={selectedMolecule.database_links.lotus}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="database-btn lotus"
                >
                  View on LOTUS
                </a>
              )}
              {selectedMolecule.database_links.npmrd && (
                <a
                  href={selectedMolecule.database_links.npmrd}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="database-btn npmrd"
                >
                  View on NPMRD
                </a>
              )}
            </div>
          )
        )}
      </div>
      
      <div className="analysis-layout">
        <div className="molecule-section">
          <MoleculeOverlay
            smiles={selectedMolecule.smiles}
            svg={moleculeSvgWithOverlays || selectedMolecule.svg}
            selectedBits={selectedBits}
          />
        </div>
        
        <div className="fingerprint-section">
          <FingerprintIndices
            retrievedIndices={retrievedFpIndices}
            bitEnvironments={bitEnvironments}
            predictedFp={predictedFp}
          />
        </div>
      </div>
      
      {predictedFp && (
        <SecondaryRetrieval
          predictedFp={predictedFp}
          retrievedFpIndices={retrievedFpIndices}
          k={10}
        />
      )}

      <Ablation />
      
      {analyzeMutation.isError && (
        <div className="error-message" style={{ margin: '20px', padding: '15px', backgroundColor: '#fee', border: '1px solid #fcc', borderRadius: '4px' }}>
          <h3>Analysis Error</h3>
          <p>{analyzeMutation.error?.message || 'Unknown error occurred'}</p>
          <button onClick={() => analyzeMutation.reset()}>Dismiss</button>
        </div>
      )}
      
      {analyzeMutation.isPending && (
        <div className="loading-message" style={{ margin: '20px', padding: '15px' }}>
          <p>Analyzing molecule...</p>
        </div>
      )}
    </div>
  )
}

export default AnalysisPage

