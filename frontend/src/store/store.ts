/**
 * Zustand store for global application state.
 */
import { create } from 'zustand'
import { ModelInfo, ResultCard } from '../services/api'

interface SpectralDataSnapshot {
  hsqc: number[]
  h_nmr: number[]
  c_nmr: number[]
  mass_spec: number[]
  mw: number | null
}

interface MainPageState {
  // Model selection
  availableModels: ModelInfo[] | null
  defaultModelId: string | null
  selectedModelId: string | null

  // Spectral input data
  hsqc: number[]
  h_nmr: number[]
  c_nmr: number[]
  mass_spec: number[]
  mw: number | null
  
  // SMILES search
  smilesInput: string
  
  // Results
  results: ResultCard[]
  predictedFp: number[] | null
  queryFp: number[] | null
  customResults: ResultCard[]

  // Retrieval MW filter
  retrievalMwMin: number | null
  retrievalMwMax: number | null
  
  // Analysis source
  analysisSource: 'prediction' | 'smiles-search' | null
  
  // Actions
  setAvailableModels: (models: ModelInfo[], defaultModelId: string) => void
  initializeModelSelection: (models: ModelInfo[], defaultModelId: string) => void
  setSelectedModelId: (modelId: string) => void

  setHSQC: (data: number[]) => void
  setHNMR: (data: number[]) => void
  setCNMR: (data: number[]) => void
  setMassSpec: (data: number[]) => void
  setMW: (mw: number | null) => void
  setSmilesInput: (smiles: string) => void
  setResults: (results: ResultCard[], predictedFp?: number[] | null, queryFp?: number[] | null) => void
  setAnalysisSource: (source: 'prediction' | 'smiles-search' | null) => void
  setRetrievalMwRange: (min: number | null, max: number | null) => void
  addCustomResult: (result: ResultCard) => void
  removeCustomResult: (index: number) => void
  clearInputs: () => void
}

interface AnalysisPageState {
  // Selected molecule
  selectedMolecule: ResultCard | null
  selectedBits: Set<number>
  retrievedFpIndices: number[]
  bitEnvironments: Record<number, any>  // Only for availability checking in FingerprintIndices
  
  // Analysis results
  moleculeSvgWithOverlays: string | null  // SVG with embedded overlays from analyze endpoint
  originalSpectralData: SpectralDataSnapshot | null
  originalPredictedFp: number[] | null
  originalSimilarityMap: string | null
  ablationSpectralData: SpectralDataSnapshot
  ablationPredictedFp: number[] | null
  ablationSimilarityMap: string | null
  ablationSeedKey: string | null
  
  // Actions
  setSelectedMolecule: (molecule: ResultCard | null) => void
  setSelectedBits: (bits: Set<number>) => void
  setRetrievedFpIndices: (indices: number[]) => void
  setBitEnvironments: (envs: Record<number, any>) => void
  setMoleculeSvgWithOverlays: (svg: string | null) => void
  setOriginalSpectralData: (data: SpectralDataSnapshot | null) => void
  setOriginalPredictedFp: (fp: number[] | null) => void
  setOriginalSimilarityMap: (map: string | null) => void
  setAblationSpectralData: (data: SpectralDataSnapshot) => void
  setAblationPredictedFp: (fp: number[] | null) => void
  setAblationSimilarityMap: (map: string | null) => void
  setAblationSeedKey: (key: string | null) => void
  resetAblation: () => void
  clearAnalysis: () => void
}

export const useMainPageStore = create<MainPageState>((set) => ({
  availableModels: null,
  defaultModelId: null,
  selectedModelId: null,

  hsqc: [],
  h_nmr: [],
  c_nmr: [],
  mass_spec: [],
  mw: null,
  smilesInput: '',
  results: [],
  predictedFp: null,
  queryFp: null,
  customResults: [],
  retrievalMwMin: null,
  retrievalMwMax: null,
  analysisSource: null,
  
  setAvailableModels: (models, defaultModelId) =>
    set({
      availableModels: models,
      defaultModelId,
    }),

  initializeModelSelection: (models, defaultModelId) => {
    const storageKey = 'marina.selectedModelId'
    const stored = (() => {
      try {
        return localStorage.getItem(storageKey)
      } catch {
        return null
      }
    })()

    const knownIds = new Set(models.map((m) => m.id))
    const selected =
      stored && knownIds.has(stored)
        ? stored
        : defaultModelId

    set({
      availableModels: models,
      defaultModelId,
      selectedModelId: selected,
    })

    try {
      localStorage.setItem(storageKey, selected)
    } catch {
      // ignore storage errors (private mode, quota, etc.)
    }
  },

  setSelectedModelId: (modelId) => {
    const storageKey = 'marina.selectedModelId'
    set({ selectedModelId: modelId })
    try {
      localStorage.setItem(storageKey, modelId)
    } catch {
      // ignore storage errors
    }
  },

  setHSQC: (data) => set({ hsqc: data }),
  setHNMR: (data) => set({ h_nmr: data }),
  setCNMR: (data) => set({ c_nmr: data }),
  setMassSpec: (data) => set({ mass_spec: data }),
  setMW: (mw) => set({ mw }),
  setSmilesInput: (smiles) => set({ smilesInput: smiles }),
  setResults: (results, predictedFp, queryFp) =>
    set({
      results,
      predictedFp: predictedFp ?? null,
      queryFp: queryFp ?? null,
    }),
  setAnalysisSource: (source) => set({ analysisSource: source }),
  setRetrievalMwRange: (min, max) => set({ retrievalMwMin: min, retrievalMwMax: max }),
  addCustomResult: (result) =>
    set((state) => ({ customResults: [...state.customResults, result] })),
  removeCustomResult: (index) =>
    set((state) => ({
      customResults: state.customResults.filter((_, i) => i !== index),
    })),
  clearInputs: () => set({
    hsqc: [],
    h_nmr: [],
    c_nmr: [],
    mass_spec: [],
    mw: null,
    smilesInput: '',
    retrievalMwMin: null,
    retrievalMwMax: null,
    customResults: [],
  }),
}))

export const useAnalysisPageStore = create<AnalysisPageState>((set) => ({
  selectedMolecule: null,
  selectedBits: new Set<number>(),
  retrievedFpIndices: [],
  bitEnvironments: {},
  moleculeSvgWithOverlays: null,
  originalSpectralData: null,
  originalPredictedFp: null,
  originalSimilarityMap: null,
  ablationSpectralData: {
    hsqc: [],
    h_nmr: [],
    c_nmr: [],
    mass_spec: [],
    mw: null,
  },
  ablationPredictedFp: null,
  ablationSimilarityMap: null,
  ablationSeedKey: null,
  
  setSelectedMolecule: (molecule) => set({ selectedMolecule: molecule }),
  setSelectedBits: (bits) => set({ selectedBits: bits }),
  setRetrievedFpIndices: (indices) => set({ retrievedFpIndices: indices }),
  setBitEnvironments: (envs) => set({ bitEnvironments: envs }),
  setMoleculeSvgWithOverlays: (svg) => set({ moleculeSvgWithOverlays: svg }),
  setOriginalSpectralData: (data) => set({ originalSpectralData: data }),
  setOriginalPredictedFp: (fp) => set({ originalPredictedFp: fp }),
  setOriginalSimilarityMap: (map) => set({ originalSimilarityMap: map }),
  setAblationSpectralData: (data) => set({ ablationSpectralData: data }),
  setAblationPredictedFp: (fp) => set({ ablationPredictedFp: fp }),
  setAblationSimilarityMap: (map) => set({ ablationSimilarityMap: map }),
  setAblationSeedKey: (key) => set({ ablationSeedKey: key }),
  resetAblation: () => set({
    ablationSpectralData: {
      hsqc: [],
      h_nmr: [],
      c_nmr: [],
      mass_spec: [],
      mw: null,
    },
    ablationPredictedFp: null,
    ablationSimilarityMap: null,
    ablationSeedKey: null,
  }),
  clearAnalysis: () => set({
    selectedMolecule: null,
    selectedBits: new Set(),
    retrievedFpIndices: [],
    bitEnvironments: {},
    moleculeSvgWithOverlays: null,
    originalSpectralData: null,
    originalPredictedFp: null,
    originalSimilarityMap: null,
    ablationSpectralData: {
      hsqc: [],
      h_nmr: [],
      c_nmr: [],
      mass_spec: [],
      mw: null,
    },
    ablationPredictedFp: null,
    ablationSimilarityMap: null,
    ablationSeedKey: null,
  }),
}))

