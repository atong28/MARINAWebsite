import { useCallback, useMemo } from 'react'
import { useAnalysisPageStore } from '../../store/store'
import './FingerprintIndices.css'

interface FingerprintIndicesProps {
  retrievedIndices: number[]
  bitEnvironments: Record<number, any>
  predictedFp?: number[] | null
}

function FingerprintIndices({ retrievedIndices, bitEnvironments, predictedFp }: FingerprintIndicesProps) {
  const { selectedBits, setSelectedBits } = useAnalysisPageStore()
  
  const predictedOverlapSet = useMemo(() => {
    if (!predictedFp) return new Set<number>()
    const set = new Set<number>()
    for (let i = 0; i < predictedFp.length; i++) {
      const val = predictedFp[i]
      if (typeof val === 'number' && val > 0.5) set.add(i)
    }
    return set
  }, [predictedFp])

  const unavailableSet = useMemo(() => {
    // Only need to compute for bits we display.
    const set = new Set<number>()
    for (const bit of retrievedIndices) {
      if (!bitEnvironments[bit]) set.add(bit)
    }
    return set
  }, [bitEnvironments, retrievedIndices])

  const handleBitClick = useCallback(
    (bit: number) => {
      const newSelected = new Set(selectedBits)
      if (newSelected.has(bit)) newSelected.delete(bit)
      else newSelected.add(bit)
      setSelectedBits(newSelected)
    },
    [selectedBits, setSelectedBits],
  )
  
  return (
    <div className="fingerprint-indices">
      <h3>Retrieved Bits ({retrievedIndices.length})</h3>
      
        <div className="bit-badges">
          {retrievedIndices.map((bit) => (
            <span
              key={`ret-${bit}`}
              className={`bit-badge ${selectedBits.has(bit) ? 'selected' : ''} ${unavailableSet.has(bit) ? 'unavailable' : ''} ${predictedOverlapSet.has(bit) ? 'overlap' : ''}`}
              onClick={() => !unavailableSet.has(bit) && handleBitClick(bit)}
              title={unavailableSet.has(bit) ? 'Bit environment not available' : `Bit ${bit}`}
            >
              {bit}
            </span>
          ))}
      </div>
    </div>
  )
}

export default FingerprintIndices

