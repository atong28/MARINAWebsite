import './RetrievalMWFilter.css'

interface RetrievalMWFilterProps {
  minMw: number | null
  maxMw: number | null
  onChange: (min: number | null, max: number | null) => void
}

function parseMw(value: string): number | null {
  if (!value.trim()) return null
  const parsed = parseFloat(value)
  if (!Number.isFinite(parsed) || parsed <= 0) return null
  return parsed
}

function RetrievalMWFilter({ minMw, maxMw, onChange }: RetrievalMWFilterProps) {
  const handleMinChange = (value: string) => {
    onChange(parseMw(value), maxMw)
  }

  const handleMaxChange = (value: string) => {
    onChange(minMw, parseMw(value))
  }

  return (
    <div className="retrieval-mw-filter">
      <span className="retrieval-mw-filter__label">Retrieval MW Filter (Da):</span>
      <div className="retrieval-mw-filter__inputs">
        <input
          type="number"
          step="0.01"
          placeholder="Min"
          value={minMw ?? ''}
          onChange={(e) => handleMinChange(e.target.value)}
        />
        <span className="retrieval-mw-filter__dash">–</span>
        <input
          type="number"
          step="0.01"
          placeholder="Max"
          value={maxMw ?? ''}
          onChange={(e) => handleMaxChange(e.target.value)}
        />
      </div>
      <span className="retrieval-mw-filter__hint">
        Leave blank for no bound. One-sided ranges are allowed.
      </span>
    </div>
  )
}

export default RetrievalMWFilter


