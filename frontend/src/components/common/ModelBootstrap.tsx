import { useEffect } from 'react'
import { useModels } from '../../services/api'
import { useMainPageStore } from '../../store/store'

/**
 * Fetches available backend models once and initializes global model selection.
 * Render this near the top of the app.
 */
function ModelBootstrap() {
  const { data, isError, error } = useModels()
  const initializeModelSelection = useMainPageStore((s) => s.initializeModelSelection)
  const setAvailableModels = useMainPageStore((s) => s.setAvailableModels)

  useEffect(() => {
    if (!data) return
    initializeModelSelection(data.models, data.default_model_id)
  }, [data, initializeModelSelection])

  useEffect(() => {
    if (!isError) return
    // If model discovery fails, keep app usable; we just won't have selector options.
    // Still, record that models are unavailable so UI can show a consistent state.
    setAvailableModels([], '')
    console.error('[ModelBootstrap] Failed to load models:', error)
  }, [isError, error, setAvailableModels])

  return null
}

export default ModelBootstrap

