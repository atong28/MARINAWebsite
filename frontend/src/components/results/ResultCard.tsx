import { ResultCard as ResultCardType } from '../../services/api'
import MoleculeViewer from './MoleculeViewer'
import './ResultCard.css'

interface ResultCardProps {
  result: ResultCardType
  position: number
  onAnalyze: () => void
  showAnalyzeButton?: boolean
}

function ResultCard({ result, position, onAnalyze, showAnalyzeButton = true }: ResultCardProps) {
  const { database_links } = result
  const hasDatabaseLinks = database_links && (
    database_links.coconut || database_links.lotus || database_links.npmrd
  )
  
  return (
    <div className="result-card">
      <div className="result-top-line">
        <span className="result-position">#{position}</span>
        <span className="similarity-badge">{(result.similarity * 100).toFixed(1)}%</span>
      </div>
      
      <div className="result-name-section">
        <h3 className="result-name">{result.name || result.smiles}</h3>
      </div>
      
      <MoleculeViewer svg={result.svg || result.plain_svg} smiles={result.smiles} />
      
      <div className="result-footer">
        <div className="database-links">
          {hasDatabaseLinks ? (
            <>
              {database_links.coconut && (
                <a
                  href={database_links.coconut}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="database-btn coconut"
                >
                  View on COCONUT
                </a>
              )}
              {database_links.lotus && (
                <a
                  href={database_links.lotus}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="database-btn lotus"
                >
                  View on LOTUS
                </a>
              )}
              {database_links.npmrd && (
                <a
                  href={database_links.npmrd}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="database-btn npmrd"
                >
                  View on NPMRD
                </a>
              )}
            </>
          ) : (
            <span className="local-dataset-note">Local dataset only</span>
          )}
        </div>
        {showAnalyzeButton && (
          <button onClick={onAnalyze} className="analyze-btn">
            Analyze
          </button>
        )}
      </div>
    </div>
  )
}

export default ResultCard

