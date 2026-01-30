/**
 * Typed API client for MARINA backend.
 */
import { useQuery, useMutation, UseQueryOptions, UseMutationOptions } from '@tanstack/react-query'

// Get API base URLs from environment variables
// In production, use relative paths so it works with nginx reverse proxy
// For development, these can be full URLs (e.g., http://localhost:5000/api/marina)
const API_BASE = import.meta.env.VITE_API_BASE || '/api'
const API_BASE_MARINA = import.meta.env.VITE_API_BASE_MARINA || `${API_BASE}`
const API_BASE_SPECTRE = import.meta.env.VITE_API_BASE_SPECTRE || `${API_BASE}`

type AbortableRequestInit = RequestInit & {
  signal?: AbortSignal
}

const modelTypeById = new Map<string, string>()

function setModelTypeCache(models: ModelInfo[]) {
  modelTypeById.clear()
  models.forEach((m) => {
    if (m.type) {
      modelTypeById.set(m.id, m.type)
    }
  })
}

function resolveApiBase(modelId?: string) {
  const modelType = modelId ? modelTypeById.get(modelId) : null
  if (modelType === 'spectre') return API_BASE_SPECTRE
  return API_BASE_MARINA
}

// Types (will be generated from OpenAPI schema in production)
export interface SpectralDataInput {
  hsqc?: number[]
  h_nmr?: number[]
  c_nmr?: number[]
  mass_spec?: number[]
  mw?: number
}

export interface PredictRequest {
  raw: SpectralDataInput
  k?: number
  mw_min?: number
  mw_max?: number
  model_id?: string
}

export interface ResultCard {
  index: number
  smiles: string
  similarity: number
  cosine_similarity?: number
  tanimoto_similarity?: number
  svg?: string
  plain_svg?: string
  name?: string
  primary_link?: string
  database_links: {
    coconut?: string
    lotus?: string
    npmrd?: string
  }
  retrieved_molecule_fp_indices: number[]
  bit_environments: Record<number, any>
  exact_mass?: number
}

export interface PredictResponse {
  results: ResultCard[]
  total_count: number
  offset: number
  limit: number
  pred_fp?: number[]
}

export interface SmilesSearchRequest {
  smiles: string
  k?: number
  mw_min?: number
  mw_max?: number
  model_id?: string
}

export interface SmilesSearchResponse {
  results: ResultCard[]
  total_count: number
  offset: number
  limit: number
  query_smiles: string
  query_fp?: number[]
}

export interface AnalysisRequest {
  smiles: string
  predicted_fp?: number[]
  retrieved_fp?: number[]
  model_id?: string
}

export interface AnalysisResponse {
  retrieved_molecule_fp_indices: number[]
  molecule_svg: string  // SVG with embedded bit environment overlays
}

export interface HealthResponse {
  status: string
  model_loaded: boolean
  message: string
  uptime_seconds: number
}

export interface ModelInfo {
  id: string
  root: string
  type: string
  default: boolean
  loaded: boolean
  display_name?: string
}

export interface ModelsResponse {
  models: ModelInfo[]
  default_model_id: string
}

export interface SecondaryRetrievalRequest {
  predicted_fp: number[]
  retrieved_fp: number[]
  k?: number
  model_id?: string
}

export interface SecondaryRetrievalResponse {
  results: ResultCard[]
  total_count: number
  difference_fp?: number[]
}

export interface CustomSmilesCardRequest {
  smiles: string
  reference_fp: number[]
  model_id?: string
}

export interface CustomSmilesCardResponse {
  result: ResultCard
}

export interface AblationRequest {
  raw: SpectralDataInput
  smiles: string
  bit_threshold?: number
  max_bits?: number
  reference_fp?: number[]
  model_id?: string
}

export interface AblationResponse {
  pred_fp: number[]
  active_bit_indices: number[]
  similarity_map?: string | null
  bit_environments: Record<number, any>
  change_overlay_svg?: string | null
}

// API functions
async function fetchJson<T>(
  endpoint: string,
  options?: AbortableRequestInit,
  baseOverride?: string
): Promise<T> {
  try {
    const base = baseOverride ?? API_BASE
    const response = await fetch(`${base}${endpoint}`, {
      ...options,
      signal: options?.signal,
      headers: {
        'Content-Type': 'application/json',
        ...options?.headers,
      },
    })

    if (!response.ok) {
      let errorMessage = `HTTP ${response.status}`
      try {
        const errorData = await response.json()
        errorMessage = errorData.detail || errorData.error || errorMessage
      } catch {
        // If response is not JSON, use status text
        errorMessage = response.statusText || errorMessage
      }
      throw new Error(errorMessage)
    }

    return response.json()
  } catch (error) {
    // Normalize AbortError into a standard Error so downstream code can ignore it if desired.
    if (error instanceof DOMException && error.name === 'AbortError') {
      throw new Error('Request aborted')
    }
    if (error instanceof Error) {
      throw error
    }
    throw new Error('Network error: Failed to fetch')
  }
}

type AbortKey =
  | 'predict'
  | 'smilesSearch'
  | 'analyze'
  | 'secondaryRetrieval'
  | 'customSmilesCard'
  | 'ablation'

const abortControllers = new Map<AbortKey, AbortController>()

function abortPreviousAndCreate(key: AbortKey): AbortController {
  const prev = abortControllers.get(key)
  if (prev) prev.abort()
  const next = new AbortController()
  abortControllers.set(key, next)
  return next
}

export function cancelInFlight(key: AbortKey) {
  const controller = abortControllers.get(key)
  if (controller) controller.abort()
}

export const api = {
  health: (modelId?: string) => fetchJson<HealthResponse>('/health', undefined, resolveApiBase(modelId)),

  models: async () => {
    const [marinaRes, spectreRes] = await Promise.allSettled([
      fetchJson<ModelsResponse>('/models', undefined, API_BASE_MARINA),
      fetchJson<ModelsResponse>('/models', undefined, API_BASE_SPECTRE),
    ])

    const models: ModelInfo[] = []
    let defaultModelId = ''

    if (marinaRes.status === 'fulfilled') {
      models.push(...marinaRes.value.models)
      defaultModelId = marinaRes.value.default_model_id || defaultModelId
    }
    if (spectreRes.status === 'fulfilled') {
      models.push(...spectreRes.value.models)
      if (!defaultModelId) defaultModelId = spectreRes.value.default_model_id
    }

    setModelTypeCache(models)
    if (!defaultModelId && models.length > 0) {
      defaultModelId = models[0]?.id ?? ''
    }

    return { models, default_model_id: defaultModelId }
  },
  
  predict: (data: PredictRequest) => {
    const controller = abortPreviousAndCreate('predict')
    return fetchJson<PredictResponse>('/predict', {
      method: 'POST',
      body: JSON.stringify(data),
      signal: controller.signal,
    }, resolveApiBase(data.model_id))
  },
  
  smilesSearch: (data: SmilesSearchRequest) => {
    const controller = abortPreviousAndCreate('smilesSearch')
    return fetchJson<SmilesSearchResponse>('/smiles-search', {
      method: 'POST',
      body: JSON.stringify(data),
      signal: controller.signal,
    }, resolveApiBase(data.model_id))
  },
  
  analyze: (data: AnalysisRequest) => {
    const controller = abortPreviousAndCreate('analyze')
    return fetchJson<AnalysisResponse>('/analyze', {
      method: 'POST',
      body: JSON.stringify(data),
      signal: controller.signal,
    }, resolveApiBase(data.model_id))
  },
  
  secondaryRetrieval: (data: SecondaryRetrievalRequest) => {
    const controller = abortPreviousAndCreate('secondaryRetrieval')
    return fetchJson<SecondaryRetrievalResponse>('/secondary-retrieval', {
      method: 'POST',
      body: JSON.stringify(data),
      signal: controller.signal,
    }, resolveApiBase(data.model_id))
  },

  customSmilesCard: (data: CustomSmilesCardRequest) => {
    const controller = abortPreviousAndCreate('customSmilesCard')
    return fetchJson<CustomSmilesCardResponse>('/custom-smiles-card', {
      method: 'POST',
      body: JSON.stringify(data),
      signal: controller.signal,
    }, resolveApiBase(data.model_id))
  },

  ablation: (data: AblationRequest) => {
    const controller = abortPreviousAndCreate('ablation')
    return fetchJson<AblationResponse>('/ablation', {
      method: 'POST',
      body: JSON.stringify(data),
      signal: controller.signal,
    }, resolveApiBase(data.model_id))
  },
}

// React Query hooks
export function useHealth(modelId?: string, options?: UseQueryOptions<HealthResponse>) {
  return useQuery({
    queryKey: ['health', modelId ?? 'default'],
    queryFn: () => api.health(modelId),
    refetchInterval: 5000, // Check every 2 seconds for responsive loading state
    refetchOnMount: true, // Always refetch on component mount (page load/refresh)
    refetchOnWindowFocus: true, // Refetch when window regains focus (overrides global setting)
    retry: 1, // Retry once on failure
    retryDelay: 1000,
    ...options,
  })
}

export function useModels(options?: UseQueryOptions<ModelsResponse>) {
  return useQuery({
    queryKey: ['models'],
    queryFn: api.models,
    staleTime: 5 * 60 * 1000,
    ...options,
  })
}

export function usePredict(options?: UseMutationOptions<PredictResponse, Error, PredictRequest>) {
  return useMutation({
    mutationFn: api.predict,
    ...options,
  })
}

export function useSmilesSearch(options?: UseMutationOptions<SmilesSearchResponse, Error, SmilesSearchRequest>) {
  return useMutation({
    mutationFn: api.smilesSearch,
    ...options,
  })
}

export function useAnalyze(options?: UseMutationOptions<AnalysisResponse, Error, AnalysisRequest>) {
  return useMutation({
    mutationFn: api.analyze,
    ...options,
  })
}

export function useSecondaryRetrieval(options?: UseMutationOptions<SecondaryRetrievalResponse, Error, SecondaryRetrievalRequest>) {
  return useMutation({
    mutationFn: api.secondaryRetrieval,
    ...options,
  })
}

export function useAblation(options?: UseMutationOptions<AblationResponse, Error, AblationRequest>) {
  return useMutation({
    mutationFn: api.ablation,
    ...options,
  })
}

