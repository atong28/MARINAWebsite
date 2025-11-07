import { ResultCard as ResultCardType } from '../../services/api'
import MoleculeViewer from './MoleculeViewer'
import './ResultCard.css'

interface ResultCardProps {
  result: ResultCardType
  position: number
  onAnalyze: () => void
}

function ResultCard({ result, position, onAnalyze }: ResultCardProps) {
  return (
    <div className="result-card">
      <div className="result-header">
        <span className="result-position">#{position}</span>
        {result.name && <h3>{result.name}</h3>}
        <span className="similarity">Similarity: {(result.similarity * 100).toFixed(1)}%</span>
      </div>
      
      <MoleculeViewer svg={result.svg || result.plain_svg} smiles={result.smiles} />
      
      <div className="result-footer">
        {result.primary_link && (
          <a href={result.primary_link} target="_blank" rel="noopener noreferrer">
            View Details
          </a>
        )}
        <button onClick={onAnalyze} className="analyze-btn">
          Analyze
        </button>
      </div>
    </div>
  )
}

export default ResultCard

