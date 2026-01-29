import { useMemo, useState } from 'react'
import { useSecondaryRetrieval } from '../../services/api'
import ResultsGrid from '../results/ResultsGrid'
import { useMainPageStore } from '../../store/store'
import './SecondaryRetrieval.css'

interface SecondaryRetrievalProps {
  predictedFp: number[] | null
  retrievedFpIndices: number[]
  k?: number
}

function SecondaryRetrieval({ predictedFp, retrievedFpIndices, k = 10 }: SecondaryRetrievalProps) {
  const [showResults, setShowResults] = useState(false)
  const selectedModelId = useMainPageStore((s) => s.selectedModelId)
  
  // Calculate difference bits: bits that are 1 in predicted_fp but 0 in retrieved_fp
  const differenceBits = useMemo(() => {
    if (!predictedFp || predictedFp.length === 0) {
      return []
    }
    
    // Create a set of retrieved indices for fast lookup
    const retrievedSet = new Set(retrievedFpIndices)
    
    // Find bits that are 1 in predicted_fp but not in retrieved_fp
    const diff: number[] = []
    for (let i = 0; i < predictedFp.length; i++) {
      const bitValue = predictedFp[i]
      if (bitValue !== undefined && bitValue > 0.5 && !retrievedSet.has(i)) {
        diff.push(i)
      }
    }
    
    return diff
  }, [predictedFp, retrievedFpIndices])
  
  // Construct retrieved_fp from indices (sparse vector)
  const retrievedFp = useMemo(() => {
    if (retrievedFpIndices.length === 0) {
      return []
    }
    
    const maxIndex = Math.max(...retrievedFpIndices)
    if (maxIndex < 0 || !isFinite(maxIndex)) {
      return []
    }
    
    const fp = new Array(maxIndex + 1).fill(0)
    retrievedFpIndices.forEach((idx) => {
      if (idx >= 0 && idx <= maxIndex) {
        fp[idx] = 1.0
      }
    })
    
    return fp
  }, [retrievedFpIndices])
  
  const secondaryRetrievalMutation = useSecondaryRetrieval({
    onSuccess: () => {
      setShowResults(true)
    },
  })
  
  const handleRetrieve = () => {
    if (!predictedFp || predictedFp.length === 0 || retrievedFp.length === 0) {
      return
    }
    
    // Ensure both fingerprints have the same length
    const maxLen = Math.max(predictedFp.length, retrievedFp.length)
    const paddedPredicted = [...predictedFp]
    const paddedRetrieved = [...retrievedFp]
    
    while (paddedPredicted.length < maxLen) {
      paddedPredicted.push(0)
    }
    while (paddedRetrieved.length < maxLen) {
      paddedRetrieved.push(0)
    }
    
    secondaryRetrievalMutation.mutate({
      predicted_fp: paddedPredicted,
      retrieved_fp: paddedRetrieved,
      k,
      model_id: selectedModelId ?? undefined,
    })
  }
  
  if (!predictedFp || predictedFp.length === 0) {
    return null
  }
  
  if (differenceBits.length === 0) {
    return (
      <div className="secondary-retrieval">
        <h3>Secondary Retrieval</h3>
        <p>No difference bits found. All predicted bits are present in the retrieved molecule.</p>
      </div>
    )
  }
  
  return (
    <div className="secondary-retrieval">
      <h3>Secondary Retrieval</h3>
      <p className="difference-bits-info">
        Difference Bits ({differenceBits.length}): Bits present in predicted fingerprint but not in retrieved molecule
      </p>
      
      <div className="bit-badges">
        {differenceBits.map((bit) => (
          <span key={`diff-${bit}`} className="bit-badge">
            {bit}
          </span>
        ))}
      </div>
      
      <button
        className="retrieve-button"
        onClick={handleRetrieve}
        disabled={secondaryRetrievalMutation.isPending}
      >
        {secondaryRetrievalMutation.isPending ? 'Retrieving...' : `Retrieve Top ${k} Compounds`}
      </button>
      
      {secondaryRetrievalMutation.isError && (
        <div className="error-message">
          <p>Error: {secondaryRetrievalMutation.error?.message || 'Unknown error occurred'}</p>
          <button onClick={() => secondaryRetrievalMutation.reset()}>Dismiss</button>
        </div>
      )}
      
      {showResults && secondaryRetrievalMutation.data && (
        <div className="secondary-results">
          <h4>Secondary Retrieval Results ({secondaryRetrievalMutation.data.total_count})</h4>
          <ResultsGrid
            results={secondaryRetrievalMutation.data.results}
            onAnalyze={() => {}}
            showAnalyzeButton={false}
          />
        </div>
      )}
    </div>
  )
}

export default SecondaryRetrieval

