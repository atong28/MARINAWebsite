import './Ablation.css'

interface SimilarityMapProps {
  title: string
  image?: string | null
  isLoading?: boolean
  description?: string
}

function SimilarityMap({ title, image, isLoading = false, description }: SimilarityMapProps) {
  return (
    <div className="similarity-map-card">
      <div className="similarity-map-header">
        <h4>{title}</h4>
        {description ? <p className="similarity-map-description">{description}</p> : null}
      </div>
      <div className="similarity-map-body">
        {isLoading ? (
          <div className="similarity-map-placeholder">Generating similarity map…</div>
        ) : image ? (
          <img src={image} alt={title} className="similarity-map-image" />
        ) : (
          <div className="similarity-map-placeholder">No similarity map available</div>
        )}
      </div>
    </div>
  )
}

export default SimilarityMap

