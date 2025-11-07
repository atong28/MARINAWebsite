import { useAnalysisPageStore } from '../../store/store'
import './FingerprintIndices.css'

interface FingerprintIndicesProps {
  retrievedIndices: number[]
  bitEnvironments: Record<number, any>
}

function FingerprintIndices({ retrievedIndices, bitEnvironments }: FingerprintIndicesProps) {
  const { selectedBits, setSelectedBits } = useAnalysisPageStore()
  
  const handleBitClick = (bit: number) => {
    const newSelected = new Set(selectedBits)
    const wasSelected = newSelected.has(bit)
    
    if (wasSelected) {
      newSelected.delete(bit)
    } else {
      newSelected.add(bit)
    }
    setSelectedBits(newSelected)
  }
  
  const isUnavailable = (bit: number) => !bitEnvironments[bit]
  const isSelected = (bit: number) => selectedBits.has(bit)
  
  return (
    <div className="fingerprint-indices">
      <h3>Retrieved Bits ({retrievedIndices.length})</h3>
      
        <div className="bit-badges">
          {retrievedIndices.map((bit) => (
            <span
              key={`ret-${bit}`}
              className={`bit-badge ${isSelected(bit) ? 'selected' : ''} ${isUnavailable(bit) ? 'unavailable' : ''}`}
              onClick={() => !isUnavailable(bit) && handleBitClick(bit)}
              title={isUnavailable(bit) ? 'Bit environment not available' : `Bit ${bit}`}
            >
              {bit}
            </span>
          ))}
      </div>
    </div>
  )
}

export default FingerprintIndices

