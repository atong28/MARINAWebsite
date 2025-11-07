import { ResultCard } from '../../services/api'
import ResultCardComponent from './ResultCard'
import './ResultsGrid.css'

interface ResultsGridProps {
  results: ResultCard[]
  onAnalyze: (index: number) => void
  showAnalyzeButton?: boolean
}

function ResultsGrid({ results, onAnalyze, showAnalyzeButton = true }: ResultsGridProps) {
  return (
    <div className="results-grid">
      <h2>Results ({results.length})</h2>
      <div className="results-list">
        {results.map((result, index) => (
          <ResultCardComponent
            key={result.index}
            result={result}
            position={index + 1}
            onAnalyze={() => onAnalyze(index)}
            showAnalyzeButton={showAnalyzeButton}
          />
        ))}
      </div>
    </div>
  )
}

export default ResultsGrid

