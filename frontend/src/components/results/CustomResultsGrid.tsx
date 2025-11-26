import { ResultCard } from '../../services/api'
import ResultCardComponent from './ResultCard'
import './ResultsGrid.css'

interface CustomResultsGridProps {
  results: ResultCard[]
  onRemove: (index: number) => void
  onAnalyze: (index: number) => void
}

function CustomResultsGrid({ results, onRemove, onAnalyze }: CustomResultsGridProps) {
  if (results.length === 0) return null

  return (
    <div className="results-grid">
      <h2>Custom Results ({results.length})</h2>
      <div className="results-list">
        {results.map((result, index) => (
          <ResultCardComponent
            key={`${result.index}-${index}-${result.smiles}`}
            result={result}
            position={0}
            onAnalyze={() => onAnalyze(index)}
            showAnalyzeButton
            isCustom
            onRemove={() => onRemove(index)}
          />
        ))}
      </div>
    </div>
  )
}

export default CustomResultsGrid


