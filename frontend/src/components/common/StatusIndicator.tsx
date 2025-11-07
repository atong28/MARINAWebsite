import { HealthResponse } from '../../services/api'
import './StatusIndicator.css'

interface StatusIndicatorProps {
  health?: HealthResponse
}

function StatusIndicator({ health }: StatusIndicatorProps) {
  const status = health?.status || 'unknown'
  const isReady = health?.model_loaded || false
  
  return (
    <div className="status-indicator">
      <i className={`fas fa-circle status-${status}`}></i>
      <span className="status-text">
        {isReady ? 'Model Ready' : 'Model Loading...'}
      </span>
    </div>
  )
}

export default StatusIndicator

