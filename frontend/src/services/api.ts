/**
 * Typed API client for MARINA backend.
 */
import { useQuery, useMutation, UseQueryOptions, UseMutationOptions } from '@tanstack/react-query'

// Get API base URL from environment variable
// VITE_API_BASE should be set in root .env file
// Falls back to localhost:5000 (default BACKEND_PORT) if not set
const API_BASE = import.meta.env.VITE_API_BASE || 'http://localhost:5000'

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
}

export interface ResultCard {
  index: number
  smiles: string
  similarity: number
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
}

export interface SmilesSearchResponse {
  results: ResultCard[]
  total_count: number
  offset: number
  limit: number
  query_smiles: string
}

export interface AnalysisRequest {
  smiles: string
  predicted_fp?: number[]
  retrieved_fp?: number[]
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

export interface SecondaryRetrievalRequest {
  predicted_fp: number[]
  retrieved_fp: number[]
  k?: number
}

export interface SecondaryRetrievalResponse {
  results: ResultCard[]
  total_count: number
  difference_fp?: number[]
}

export interface AblationRequest {
  raw: SpectralDataInput
  smiles: string
  bit_threshold?: number
  max_bits?: number
  reference_fp?: number[]
}

export interface AblationResponse {
  pred_fp: number[]
  active_bit_indices: number[]
  similarity_map?: string | null
  bit_environments: Record<number, any>
  change_overlay_svg?: string | null
}

// API functions
async function fetchJson<T>(endpoint: string, options?: RequestInit): Promise<T> {
  try {
    const response = await fetch(`${API_BASE}${endpoint}`, {
      ...options,
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
    if (error instanceof Error) {
      throw error
    }
    throw new Error('Network error: Failed to fetch')
  }
}

export const api = {
  health: () => fetchJson<HealthResponse>('/health'),
  
  predict: (data: PredictRequest) =>
    fetchJson<PredictResponse>('/predict', {
      method: 'POST',
      body: JSON.stringify(data),
    }),
  
  smilesSearch: (data: SmilesSearchRequest) =>
    fetchJson<SmilesSearchResponse>('/smiles-search', {
      method: 'POST',
      body: JSON.stringify(data),
    }),
  
  analyze: (data: AnalysisRequest) =>
    fetchJson<AnalysisResponse>('/analyze', {
      method: 'POST',
      body: JSON.stringify(data),
    }),
  
  secondaryRetrieval: (data: SecondaryRetrievalRequest) =>
    fetchJson<SecondaryRetrievalResponse>('/secondary-retrieval', {
      method: 'POST',
      body: JSON.stringify(data),
    }),

  ablation: (data: AblationRequest) =>
    fetchJson<AblationResponse>('/ablation', {
      method: 'POST',
      body: JSON.stringify(data),
    }),
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

