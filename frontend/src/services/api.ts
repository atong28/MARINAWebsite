/**
 * Typed API client for MARINA backend.
 */
import { useQuery, useMutation, UseQueryOptions, UseMutationOptions } from '@tanstack/react-query'

// Get API base URL from environment variable
// In production, use relative path '/api' so it works with nginx reverse proxy
// For development, VITE_API_BASE can be set to full URL (e.g., http://localhost:5000/api)
const API_BASE = import.meta.env.VITE_API_BASE || '/api'

type AbortableRequestInit = RequestInit & { signal?: AbortSignal }

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
}

export interface ModelsResponse {
  models: ModelInfo[]
  default_model_id: string
}

export interface LoadModelResponse {
  model_id: string
  status: 'loaded' | 'already_loaded' | 'error'
  message: string
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
async function fetchJson<T>(endpoint: string, options?: AbortableRequestInit): Promise<T> {
  try {
    const response = await fetch(`${API_BASE}${endpoint}`, {
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
  health: () => fetchJson<HealthResponse>('/health'),

  models: () => fetchJson<ModelsResponse>('/models'),

  loadModel: (modelId: string) => fetchJson<LoadModelResponse>(`/models/${modelId}/load`, {
    method: 'POST',
  }),
  
  predict: (data: PredictRequest) => {
    const controller = abortPreviousAndCreate('predict')
    return fetchJson<PredictResponse>('/predict', {
      method: 'POST',
      body: JSON.stringify(data),
      signal: controller.signal,
    })
  },
  
  smilesSearch: (data: SmilesSearchRequest) => {
    const controller = abortPreviousAndCreate('smilesSearch')
    return fetchJson<SmilesSearchResponse>('/smiles-search', {
      method: 'POST',
      body: JSON.stringify(data),
      signal: controller.signal,
    })
  },
  
  analyze: (data: AnalysisRequest) => {
    const controller = abortPreviousAndCreate('analyze')
    return fetchJson<AnalysisResponse>('/analyze', {
      method: 'POST',
      body: JSON.stringify(data),
      signal: controller.signal,
    })
  },
  
  secondaryRetrieval: (data: SecondaryRetrievalRequest) => {
    const controller = abortPreviousAndCreate('secondaryRetrieval')
    return fetchJson<SecondaryRetrievalResponse>('/secondary-retrieval', {
      method: 'POST',
      body: JSON.stringify(data),
      signal: controller.signal,
    })
  },

  customSmilesCard: (data: CustomSmilesCardRequest) => {
    const controller = abortPreviousAndCreate('customSmilesCard')
    return fetchJson<CustomSmilesCardResponse>('/custom-smiles-card', {
      method: 'POST',
      body: JSON.stringify(data),
      signal: controller.signal,
    })
  },

  ablation: (data: AblationRequest) => {
    const controller = abortPreviousAndCreate('ablation')
    return fetchJson<AblationResponse>('/ablation', {
      method: 'POST',
      body: JSON.stringify(data),
      signal: controller.signal,
    })
  },
}

// React Query hooks
export function useHealth(options?: UseQueryOptions<HealthResponse>) {
  return useQuery({
    queryKey: ['health'],
    queryFn: api.health,
    refetchInterval: 30000, // Check every 30 seconds
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

export function useLoadModel(options?: UseMutationOptions<LoadModelResponse, Error, string>) {
  return useMutation({
    mutationFn: api.loadModel,
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

