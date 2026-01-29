import { useState, useEffect, useRef, useCallback } from 'react'
import { useNavigate } from 'react-router-dom'
import { usePredict, useSmilesSearch, useHealth, api, cancelInFlight } from '../services/api'
import { useMainPageStore } from '../store/store'
import { getAvailableExamples, loadExample, type ExampleMetadata } from '../services/exampleLoader'
import { ROUTES } from '../routes'
import SpectralInputTabs from '../components/spectral/SpectralInputTabs'
import SmilesSearchBox from '../components/spectral/SmilesSearchBox'
import CustomSmilesBox from '../components/spectral/CustomSmilesBox'
import ResultsGrid from '../components/results/ResultsGrid'
import CustomResultsGrid from '../components/results/CustomResultsGrid'
import StatusIndicator from '../components/common/StatusIndicator'
import ModelSelector from '../components/common/ModelSelector'
import './MainPage.css'

const API_DOCS_URL = '/docs'

function MainPage() {
  const navigate = useNavigate()
  const [k, setK] = useState(10)
  const [activeTab, setActiveTab] = useState<'spectral' | 'smiles'>('spectral')
  const [availableExamples, setAvailableExamples] = useState<ExampleMetadata[]>([])
  const [selectedExample, setSelectedExample] = useState<string>('')
  const [isLoadingExamples, setIsLoadingExamples] = useState(false)
  const [hasInvalidSpreadsheet, setHasInvalidSpreadsheet] = useState(false)
  
  const {
    selectedModelId,
    hsqc,
    h_nmr,
    c_nmr,
    mass_spec,
    mw,
    smilesInput,
    results,
    setResults,
    setAnalysisSource,
    setSmilesInput,
    setHSQC,
    setHNMR,
    setCNMR,
    setMassSpec,
    setMW,
    retrievalMwMin,
    retrievalMwMax,
    customResults,
    addCustomResult,
    removeCustomResult,
  } = useMainPageStore()
  const [customSmiles, setCustomSmiles] = useState('')
  const [customResultError, setCustomResultError] = useState<string | null>(null)
  const [isCustomLoading, setIsCustomLoading] = useState(false)

  // Monotonic request sequencing to ensure only latest request updates UI state.
  const requestSeqRef = useRef(0)
  const latestPredictSeqRef = useRef<number | null>(null)
  const latestSmilesSeqRef = useRef<number | null>(null)
  
  const { data: health } = useHealth()
  const predictMutation = usePredict({
    onSuccess: (data) => {
      if (latestPredictSeqRef.current === null) return
      if (latestPredictSeqRef.current !== requestSeqRef.current) return
      setResults(data.results, data.pred_fp || null, null)
      setAnalysisSource('prediction')
    },
  })
  const smilesSearchMutation = useSmilesSearch({
    onSuccess: (data) => {
      if (latestSmilesSeqRef.current === null) return
      if (latestSmilesSeqRef.current !== requestSeqRef.current) return
      setResults(data.results, null, data.query_fp || null)
      setAnalysisSource('smiles-search')
    },
  })

  // Load available examples on mount
  useEffect(() => {
    const loadExamples = async () => {
      try {
        const examples = await getAvailableExamples()
        setAvailableExamples(examples)
      } catch (error) {
        console.error('Failed to load examples:', error)
      }
    }
    loadExamples()
  }, [])
  
  // Helper function to filter out NaN and null values from arrays
  const filterNaNValues = (arr: number[]): number[] => {
    return arr.filter((val) => {
      // Filter out NaN, null, and undefined values
      return val !== null && val !== undefined && !Number.isNaN(val) && isFinite(val)
    })
  }

  const mwInvalid =
    (retrievalMwMin !== null && retrievalMwMin !== undefined && retrievalMwMin < 0) ||
    (retrievalMwMax !== null && retrievalMwMax !== undefined && retrievalMwMax < 0) ||
    (retrievalMwMin !== null &&
      retrievalMwMax !== null &&
      retrievalMwMin !== undefined &&
      retrievalMwMax !== undefined &&
      retrievalMwMin > retrievalMwMax)
  
  const handlePredict = useCallback(() => {
    // New request starts: cancel other in-flight search to avoid stale writes/races.
    cancelInFlight('smilesSearch')
    cancelInFlight('predict')
    requestSeqRef.current += 1
    latestPredictSeqRef.current = requestSeqRef.current

    setCustomResultError(null)

    // Sanitize arrays by filtering out NaN values
    const sanitizedHSQC = filterNaNValues(hsqc)
    const sanitizedHNMR = filterNaNValues(h_nmr)
    const sanitizedCNMR = filterNaNValues(c_nmr)
    const sanitizedMassSpec = filterNaNValues(mass_spec)
    
    // Build request payload - only include arrays with at least one valid value
    const raw: any = {}
    if (sanitizedHSQC.length > 0) raw.hsqc = sanitizedHSQC
    if (sanitizedHNMR.length > 0) raw.h_nmr = sanitizedHNMR
    if (sanitizedCNMR.length > 0) raw.c_nmr = sanitizedCNMR
    if (sanitizedMassSpec.length > 0) raw.mass_spec = sanitizedMassSpec
    if (mw !== null && mw !== undefined && !Number.isNaN(mw) && isFinite(mw)) {
      raw.mw = mw
    }
    
    const payload: any = { raw, k }
    if (selectedModelId) payload.model_id = selectedModelId
    if (retrievalMwMin !== null && retrievalMwMin !== undefined && !Number.isNaN(retrievalMwMin) && isFinite(retrievalMwMin)) {
      payload.mw_min = retrievalMwMin
    }
    if (retrievalMwMax !== null && retrievalMwMax !== undefined && !Number.isNaN(retrievalMwMax) && isFinite(retrievalMwMax)) {
      payload.mw_max = retrievalMwMax
    }
    predictMutation.mutate(payload)
  }, [
    c_nmr,
    h_nmr,
    hsqc,
    mass_spec,
    mw,
    k,
    retrievalMwMin,
    retrievalMwMax,
    selectedModelId,
    predictMutation,
    setCustomResultError,
  ])
  
  const handleSmilesSearch = useCallback(() => {
    if (!smilesInput.trim()) return
    cancelInFlight('predict')
    cancelInFlight('smilesSearch')
    requestSeqRef.current += 1
    latestSmilesSeqRef.current = requestSeqRef.current

    setCustomResultError(null)

    const payload: any = { smiles: smilesInput.trim(), k }
    if (selectedModelId) payload.model_id = selectedModelId
    if (retrievalMwMin !== null && retrievalMwMin !== undefined && !Number.isNaN(retrievalMwMin) && isFinite(retrievalMwMin)) {
      payload.mw_min = retrievalMwMin
    }
    if (retrievalMwMax !== null && retrievalMwMax !== undefined && !Number.isNaN(retrievalMwMax) && isFinite(retrievalMwMax)) {
      payload.mw_max = retrievalMwMax
    }
    smilesSearchMutation.mutate(payload)
  }, [
    smilesInput,
    k,
    retrievalMwMin,
    retrievalMwMax,
    selectedModelId,
    smilesSearchMutation,
    setCustomResultError,
  ])
  
  const handleAnalyze = useCallback((index: number) => {
    const result = results[index]
    if (!result) return
    navigate(ROUTES.ANALYSIS, { state: { result, index } })
  }, [navigate, results])

  const handleCustomAnalyze = useCallback((index: number) => {
    const result = customResults[index]
    if (!result) return
    // Index is not meaningful for custom results, but AnalysisPage only needs the result itself.
    navigate(ROUTES.ANALYSIS, { state: { result, index: 0 } })
  }, [customResults, navigate])

  const handleLoadExample = async () => {
    if (!selectedExample) return
    
    setIsLoadingExamples(true)
    try {
      const exampleData = await loadExample(selectedExample)
      setHSQC(exampleData.hsqc || [])
      setHNMR(exampleData.h_nmr || [])
      setCNMR(exampleData.c_nmr || [])
      setMassSpec(exampleData.mass_spec || [])
      setMW(exampleData.mw || null)
    } catch (error) {
      console.error('Failed to load example:', error)
      alert(`Failed to load example: ${error instanceof Error ? error.message : 'Unknown error'}`)
    } finally {
      setIsLoadingExamples(false)
    }
  }
  
  return (
    <div className="main-page">
      <header className="header">
        <div className="header-top">
          <a
            className="api-docs-btn"
            href={API_DOCS_URL}
            target="_blank"
            rel="noopener noreferrer"
          >
            API Docs
          </a>
          <ModelSelector />
          <StatusIndicator health={health} />
        </div>
        <h1>
          <i className="fas fa-flask"></i> MARINA: Molecular Structure Annotator
        </h1>
        <p className="subtitle">Predict molecular structures from spectral data using deep learning</p>
      </header>
      
      <div className="tabs">
        <button
          className={activeTab === 'spectral' ? 'active' : ''}
          onClick={() => setActiveTab('spectral')}
        >
          Spectral Data
        </button>
        <button
          className={activeTab === 'smiles' ? 'active' : ''}
          onClick={() => setActiveTab('smiles')}
        >
          SMILES Search
        </button>
      </div>
      
      {activeTab === 'spectral' && (
        <div className="spectral-input-section">
          <SpectralInputTabs onValidationChange={(s) => setHasInvalidSpreadsheet(s.anyInvalid)} />
          <div className="controls">
            <label>
              Number of results:
              <input
                type="number"
                min="1"
                max="50"
                value={k}
                onChange={(e) => setK(parseInt(e.target.value) || 10)}
              />
            </label>
            <button
              onClick={handlePredict}
              disabled={predictMutation.isPending || !health?.model_loaded || hasInvalidSpreadsheet || mwInvalid}
              title={
                hasInvalidSpreadsheet
                  ? 'You have incomplete data in the spreadsheet.'
                  : mwInvalid
                  ? 'Invalid MW filter: ensure 0 ≤ min ≤ max.'
                  : ''
              }
            >
              {predictMutation.isPending ? 'Predicting...' : 'Predict Structure'}
            </button>
            {hasInvalidSpreadsheet && (
              <span style={{ color: '#c00', fontSize: 12 }}>Incomplete data detected. Please fix or use Condense rows.</span>
            )}
            <div className="example-loader">
              <select
                value={selectedExample}
                onChange={(e) => setSelectedExample(e.target.value)}
                disabled={isLoadingExamples || availableExamples.length === 0}
              >
                <option value="">Select Example...</option>
                {availableExamples.map((example) => (
                  <option key={example.filename} value={example.filename}>
                    {example.name}
                  </option>
                ))}
              </select>
              <button
                onClick={handleLoadExample}
                disabled={!selectedExample || isLoadingExamples}
              >
                {isLoadingExamples ? 'Loading...' : 'Load Example'}
              </button>
            </div>
          </div>
        </div>
      )}
      
      {activeTab === 'smiles' && (
        <div className="smiles-input-section">
          <SmilesSearchBox
            value={smilesInput}
            onChange={setSmilesInput}
          />
          <div className="controls">
            <label>
              Number of results:
              <input
                type="number"
                min="1"
                max="50"
                value={k}
                onChange={(e) => setK(parseInt(e.target.value) || 10)}
              />
            </label>
            <button
              onClick={handleSmilesSearch}
              disabled={smilesSearchMutation.isPending || !health?.model_loaded || mwInvalid}
            >
              {smilesSearchMutation.isPending ? 'Searching...' : 'Search'}
            </button>
          </div>
        </div>
      )}
      
      {results.length > 0 && (
        <>
          <div className="smiles-input-section custom-smiles-section">
            <CustomSmilesBox value={customSmiles} onChange={setCustomSmiles} />
            <div className="controls">
              <button
                onClick={async () => {
                  setCustomResultError(null)
                  const trimmed = customSmiles.trim()
                  if (!trimmed) return
                  try {
                    setIsCustomLoading(true)
                    const state = useMainPageStore.getState()
                    const refFp =
                      state.analysisSource === 'prediction'
                        ? state.predictedFp
                        : state.queryFp
                    if (!refFp || refFp.length === 0) {
                      setCustomResultError('No reference fingerprint available from the current session.')
                      return
                    }
                    const response = await api.customSmilesCard({
                      smiles: trimmed,
                      reference_fp: refFp,
                      model_id: state.selectedModelId ?? undefined,
                    })
                    if (response.result) {
                      addCustomResult(response.result)
                    }
                  } catch (error) {
                    const message = error instanceof Error ? error.message : 'Failed to fetch custom SMILES result'
                    setCustomResultError(message)
                  } finally {
                    setIsCustomLoading(false)
                  }
                }}
                disabled={isCustomLoading}
              >
                {isCustomLoading ? 'Adding...' : 'Add Custom Result'}
              </button>
            </div>
            {customResultError && (
              <div className="error-message" style={{ marginTop: 10 }}>
                {customResultError}
              </div>
            )}
          </div>

          <CustomResultsGrid
            results={customResults}
            onRemove={(index) => removeCustomResult(index)}
            onAnalyze={handleCustomAnalyze}
          />

          <ResultsGrid results={results} onAnalyze={handleAnalyze} />
        </>
      )}
      
      {(predictMutation.error || smilesSearchMutation.error) && (
        <div className="error-message">
          Error: {predictMutation.error?.message || smilesSearchMutation.error?.message}
        </div>
      )}
    </div>
  )
}

export default MainPage

