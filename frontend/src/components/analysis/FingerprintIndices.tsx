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
      console.log('[FingerprintIndices] Deselecting bit:', bit, {
        previousSelected: Array.from(selectedBits),
        newSelectedCount: newSelected.size,
        newSelectedBits: Array.from(newSelected),
      })
    } else {
      newSelected.add(bit)
      console.log('[FingerprintIndices] Selecting bit:', bit, {
        previousSelected: Array.from(selectedBits),
        newSelectedCount: newSelected.size,
        newSelectedBits: Array.from(newSelected),
        isUnavailable: isUnavailable(bit),
      })
    }
    
    console.log('[FingerprintIndices] Calling setSelectedBits with:', Array.from(newSelected))
    setSelectedBits(newSelected)
    
    // Verify the state was updated (this will log on next render)
    setTimeout(() => {
      console.log('[FingerprintIndices] State after update - selectedBits:', Array.from(selectedBits))
    }, 0)
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

