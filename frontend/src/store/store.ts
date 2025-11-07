/**
 * Zustand store for global application state.
 */
import { create } from 'zustand'
import { ResultCard } from '../services/api'

interface MainPageState {
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
  
  // Analysis source
  analysisSource: 'prediction' | 'smiles-search' | null
  
  // Actions
  setHSQC: (data: number[]) => void
  setHNMR: (data: number[]) => void
  setCNMR: (data: number[]) => void
  setMassSpec: (data: number[]) => void
  setMW: (mw: number | null) => void
  setSmilesInput: (smiles: string) => void
  setResults: (results: ResultCard[], predictedFp?: number[] | null) => void
  setAnalysisSource: (source: 'prediction' | 'smiles-search' | null) => void
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
  
  // Actions
  setSelectedMolecule: (molecule: ResultCard | null) => void
  setSelectedBits: (bits: Set<number>) => void
  setRetrievedFpIndices: (indices: number[]) => void
  setBitEnvironments: (envs: Record<number, any>) => void
  setMoleculeSvgWithOverlays: (svg: string | null) => void
  clearAnalysis: () => void
}

export const useMainPageStore = create<MainPageState>((set) => ({
  hsqc: [],
  h_nmr: [],
  c_nmr: [],
  mass_spec: [],
  mw: null,
  smilesInput: '',
  results: [],
  predictedFp: null,
  analysisSource: null,
  
  setHSQC: (data) => set({ hsqc: data }),
  setHNMR: (data) => set({ h_nmr: data }),
  setCNMR: (data) => set({ c_nmr: data }),
  setMassSpec: (data) => set({ mass_spec: data }),
  setMW: (mw) => set({ mw }),
  setSmilesInput: (smiles) => set({ smilesInput: smiles }),
  setResults: (results, predictedFp) => set({ results, predictedFp }),
  setAnalysisSource: (source) => set({ analysisSource: source }),
  clearInputs: () => set({
    hsqc: [],
    h_nmr: [],
    c_nmr: [],
    mass_spec: [],
    mw: null,
    smilesInput: '',
  }),
}))

export const useAnalysisPageStore = create<AnalysisPageState>((set) => ({
  selectedMolecule: null,
  selectedBits: new Set<number>(),
  retrievedFpIndices: [],
  bitEnvironments: {},
  moleculeSvgWithOverlays: null,
  
  setSelectedMolecule: (molecule) => set({ selectedMolecule: molecule }),
  setSelectedBits: (bits) => set({ selectedBits: bits }),
  setRetrievedFpIndices: (indices) => set({ retrievedFpIndices: indices }),
  setBitEnvironments: (envs) => set({ bitEnvironments: envs }),
  setMoleculeSvgWithOverlays: (svg) => set({ moleculeSvgWithOverlays: svg }),
  clearAnalysis: () => set({
    selectedMolecule: null,
    selectedBits: new Set(),
    retrievedFpIndices: [],
    bitEnvironments: {},
    moleculeSvgWithOverlays: null,
  }),
}))

