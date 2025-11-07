import { useState } from 'react'
import { useNavigate } from 'react-router-dom'
import { usePredict, useSmilesSearch, useHealth } from '../services/api'
import { useMainPageStore } from '../store/store'
import SpectralInputTabs from '../components/spectral/SpectralInputTabs'
import ResultsGrid from '../components/results/ResultsGrid'
import StatusIndicator from '../components/common/StatusIndicator'
import './MainPage.css'

function MainPage() {
  const navigate = useNavigate()
  const [k, setK] = useState(10)
  const [activeTab, setActiveTab] = useState<'spectral' | 'smiles'>('spectral')
  
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
  
  return (
    <div className="main-page">
      <header className="header">
        <h1>
          <i className="fas fa-flask"></i> MARINA: Molecular Structure Annotator
        </h1>
        <p className="subtitle">Predict molecular structures from spectral data using deep learning</p>
        <StatusIndicator health={health} />
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

