import { useState } from 'react'
import { useQueryClient } from '@tanstack/react-query'
import { useMainPageStore } from '../../store/store'
import { useLoadModel } from '../../services/api'
import './ModelSelector.css'

interface ModelSelectorProps {
  compact?: boolean
}

function ModelSelector({ compact = true }: ModelSelectorProps) {
  const { availableModels, selectedModelId, defaultModelId, setSelectedModelId } = useMainPageStore()
  const queryClient = useQueryClient()
  const [loadError, setLoadError] = useState<string | null>(null)

  const selected = selectedModelId ?? defaultModelId ?? availableModels?.[0]?.id ?? ''
  const selectedModel = availableModels?.find(m => m.id === selected)
  const needsLoading = selectedModel && !selectedModel.loaded

  const loadModelMutation = useLoadModel({
    onSuccess: () => {
      setLoadError(null)
      // Invalidate models query to refresh loaded status
      queryClient.invalidateQueries({ queryKey: ['models'] })
    },
    onError: (error) => {
      setLoadError(error.message || 'Failed to load model')
    },
  })

  const handleLoad = () => {
    if (!selected) return
    setLoadError(null)
    loadModelMutation.mutate(selected)
  }

  if (!availableModels || availableModels.length === 0) {
    return (
      <div className={`model-selector ${compact ? 'compact' : ''}`}>
        <span className="model-selector-label">Model</span>
        <span className="model-selector-loading">Loading…</span>
      </div>
    )
  }

  return (
    <div className={`model-selector ${compact ? 'compact' : ''}`}>
      <label className="model-selector-label" htmlFor="model-selector">
        Model
      </label>
      <select
        id="model-selector"
        className="model-selector-select"
        value={selected}
        onChange={(e) => {
          setSelectedModelId(e.target.value)
          setLoadError(null)
        }}
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
      {needsLoading && (
        <div className="model-selector-load-container">
          <button
            className="model-selector-load-button"
            onClick={handleLoad}
            disabled={loadModelMutation.isPending}
            title="Load model and preload resources for faster predictions"
          >
            {loadModelMutation.isPending ? 'Loading...' : 'Load'}
          </button>
          {loadError && (
            <span className="model-selector-load-error" title={loadError}>
              ⚠
            </span>
          )}
        </div>
      )}
    </div>
  )
}

export default ModelSelector

