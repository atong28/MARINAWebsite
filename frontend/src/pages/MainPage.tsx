import { useState, useEffect } from 'react'
import { useNavigate } from 'react-router-dom'
import { usePredict, useSmilesSearch, useHealth } from '../services/api'
import { useMainPageStore } from '../store/store'
import { getAvailableExamples, loadExample, type ExampleMetadata } from '../services/exampleLoader'
import SpectralInputTabs from '../components/spectral/SpectralInputTabs'
import ResultsGrid from '../components/results/ResultsGrid'
import StatusIndicator from '../components/common/StatusIndicator'
import './MainPage.css'

const API_BASE = (import.meta.env.VITE_API_BASE as string | undefined) || 'http://localhost:5000'
const API_DOCS_URL = `${API_BASE.replace(/\/$/, '')}/docs`

function MainPage() {
  const navigate = useNavigate()
  const [k, setK] = useState(10)
  const [activeTab, setActiveTab] = useState<'spectral' | 'smiles'>('spectral')
  const [availableExamples, setAvailableExamples] = useState<ExampleMetadata[]>([])
  const [selectedExample, setSelectedExample] = useState<string>('')
  const [isLoadingExamples, setIsLoadingExamples] = useState(false)
  
  const {
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
  } = useMainPageStore()
  
  const { data: health } = useHealth()
  const predictMutation = usePredict({
    onSuccess: (data) => {
      setResults(data.results, data.pred_fp || null)
      setAnalysisSource('prediction')
    },
  })
  const smilesSearchMutation = useSmilesSearch({
    onSuccess: (data) => {
      setResults(data.results)
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
  
  const handlePredict = () => {
    const raw: any = {}
    if (hsqc.length > 0) raw.hsqc = hsqc
    if (h_nmr.length > 0) raw.h_nmr = h_nmr
    if (c_nmr.length > 0) raw.c_nmr = c_nmr
    if (mass_spec.length > 0) raw.mass_spec = mass_spec
    if (mw !== null) raw.mw = mw
    
    predictMutation.mutate({ raw, k })
  }
  
  const handleSmilesSearch = () => {
    if (!smilesInput.trim()) return
    smilesSearchMutation.mutate({ smiles: smilesInput.trim(), k })
  }
  
  const handleAnalyze = (index: number) => {
    const result = results[index]
    if (!result) return
    navigate('/analysis', { state: { result, index } })
  }

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
          <SpectralInputTabs />
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
              disabled={predictMutation.isPending || !health?.model_loaded}
            >
              {predictMutation.isPending ? 'Predicting...' : 'Predict Structure'}
            </button>
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
          <label>
            SMILES String:
            <input
              type="text"
              value={smilesInput}
              onChange={(e) => setSmilesInput(e.target.value)}
              placeholder="Enter SMILES string..."
            />
          </label>
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
              disabled={smilesSearchMutation.isPending || !health?.model_loaded}
            >
              {smilesSearchMutation.isPending ? 'Searching...' : 'Search'}
            </button>
          </div>
        </div>
      )}
      
      {results.length > 0 && (
        <ResultsGrid results={results} onAnalyze={handleAnalyze} />
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

