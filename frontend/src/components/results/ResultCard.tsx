import { memo, useMemo } from 'react'
import { ResultCard as ResultCardType } from '../../services/api'
import MoleculeViewer from './MoleculeViewer'
import './ResultCard.css'

interface ResultCardProps {
  result: ResultCardType
  position: number
  onAnalyze: () => void
  showAnalyzeButton?: boolean
  isCustom?: boolean
  onRemove?: () => void
}

function ResultCard({
  result,
  position,
  onAnalyze,
  showAnalyzeButton = true,
  isCustom = false,
  onRemove,
}: ResultCardProps) {
  const { database_links } = result
  const hasDatabaseLinks = useMemo(() => {
    return Boolean(
      database_links &&
        (database_links.coconut || database_links.lotus || database_links.npmrd),
    )
  }, [database_links])

  const cosine = useMemo(() => {
    return typeof result.cosine_similarity === 'number' ? result.cosine_similarity : result.similarity
  }, [result.cosine_similarity, result.similarity])

  const tanimoto = useMemo(() => {
    return typeof result.tanimoto_similarity === 'number' ? result.tanimoto_similarity : undefined
  }, [result.tanimoto_similarity])
  const positionLabel = position > 0 ? `#${position}` : 'Custom'
  
  return (
    <div className="result-card">
      <div className="result-top-line">
        <span className="result-position">{positionLabel}</span>
        <div className="similarity-badges">
          <span className="similarity-badge cosine-badge">
            C:{cosine.toFixed(3)}
          </span>
          {typeof tanimoto === 'number' && (
            <span className="similarity-badge tanimoto-badge">
              T:{tanimoto.toFixed(3)}
            </span>
          )}
        </div>
      </div>
      
      <div className="result-name-section">
        <h3 className="result-name">{result.name || result.smiles}</h3>
        {result.exact_mass && (
          <div className="result-exact-mass">
            Exact Mass: {result.exact_mass.toFixed(4)} Da
          </div>
        )}
      </div>
      
      <MoleculeViewer svg={result.svg || result.plain_svg} smiles={result.smiles} />
      
      <div className="result-footer">
        <div className="database-links">
          {isCustom ? (
            onRemove && (
              <button className="custom-remove-btn" onClick={onRemove}>
                ✕ Remove
              </button>
            )
          ) : hasDatabaseLinks ? (
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

export default memo(ResultCard)

