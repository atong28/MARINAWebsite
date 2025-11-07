import { useMemo } from 'react'
import './Ablation.css'

type BitStatus = 'unchanged' | 'gained' | 'lost'

interface BitComparisonProps {
  originalFp?: number[] | null
  newFp?: number[] | null
  threshold?: number
  onBitClick?: (bit: number, status: BitStatus) => void
}

interface BitGroups {
  unchanged: number[]
  gained: number[]
  lost: number[]
}

function BitComparison({ originalFp, newFp, threshold = 0.5, onBitClick }: BitComparisonProps) {
  const groups = useMemo<BitGroups>(() => {
    const base: BitGroups = {
      unchanged: [],
      gained: [],
      lost: [],
    }

    if (!originalFp || !newFp || originalFp.length === 0 || newFp.length === 0) {
      return base
    }

    const originalBits = new Set<number>()
    const newBits = new Set<number>()

    originalFp.forEach((value, index) => {
      if (value !== undefined && value > threshold) {
        originalBits.add(index)
      }
    })

    newFp.forEach((value, index) => {
      if (value !== undefined && value > threshold) {
        newBits.add(index)
      }
    })

    originalBits.forEach((bit) => {
      if (newBits.has(bit)) {
        base.unchanged.push(bit)
      } else {
        base.lost.push(bit)
      }
    })

    newBits.forEach((bit) => {
      if (!originalBits.has(bit)) {
        base.gained.push(bit)
      }
    })

    base.unchanged.sort((a, b) => a - b)
    base.gained.sort((a, b) => a - b)
    base.lost.sort((a, b) => a - b)

    return base
  }, [originalFp, newFp, threshold])

  const renderBitGroup = (label: string, bits: number[], status: BitStatus, className: string) => (
    <div className={`bit-comparison-group ${className}`}>
      <h4>{label}</h4>
      <div className="bit-badges">
        {bits.length > 0 ? (
          bits.map((bit) => (
            <button
              type="button"
              key={`${status}-${bit}`}
              className={`bit-badge ${status}`}
              onClick={() => onBitClick?.(bit, status)}
            >
              {bit}
            </button>
          ))
        ) : (
          <span className="bit-badge empty">None</span>
        )}
      </div>
    </div>
  )

  return (
    <div className="bit-comparison">
      {renderBitGroup('Unchanged Bits', groups.unchanged, 'unchanged', 'unchanged')}
      {renderBitGroup('Gained Bits', groups.gained, 'gained', 'gained')}
      {renderBitGroup('Lost Bits', groups.lost, 'lost', 'lost')}
    </div>
  )
}

export default BitComparison

