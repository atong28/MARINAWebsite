import { ResultCard } from '../../services/api'
import ResultCardComponent from './ResultCard'
import './ResultsGrid.css'

interface ResultsGridProps {
  results: ResultCard[]
  onAnalyze: (index: number) => void
}

function ResultsGrid({ results, onAnalyze }: ResultsGridProps) {
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
          />
        ))}
      </div>
    </div>
  )
}

export default ResultsGrid

