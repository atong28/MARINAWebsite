import { useMainPageStore } from '../../store/store'
import './ModelSelector.css'

interface ModelSelectorProps {
  compact?: boolean
}

function ModelSelector({ compact = true }: ModelSelectorProps) {
  const { availableModels, selectedModelId, defaultModelId, setSelectedModelId } = useMainPageStore()

  if (!availableModels || availableModels.length === 0) {
    return (
      <div className={`model-selector ${compact ? 'compact' : ''}`}>
        <span className="model-selector-label">Model</span>
        <span className="model-selector-loading">Loading…</span>
      </div>
    )
  }

  const selected = selectedModelId ?? defaultModelId ?? availableModels[0]?.id ?? ''

  return (
    <div className={`model-selector ${compact ? 'compact' : ''}`}>
      <label className="model-selector-label" htmlFor="model-selector">
        Model
      </label>
      <select
        id="model-selector"
        className="model-selector-select"
        value={selected}
        onChange={(e) => setSelectedModelId(e.target.value)}
      >
        {availableModels.map((m) => {
          const suffixParts: string[] = []
          if (m.default) suffixParts.push('default')
          if (!m.loaded) suffixParts.push('not loaded')
          if (m.type) suffixParts.push(m.type)
          const suffix = suffixParts.length ? ` (${suffixParts.join(', ')})` : ''
          return (
            <option key={m.id} value={m.id}>
              {m.id}
              {suffix}
            </option>
          )
        })}
      </select>
    </div>
  )
}

export default ModelSelector

