/**
 * Service for loading spectral data examples.
 */

export interface ExampleData {
  name: string
  hsqc: number[]
  h_nmr: number[]
  c_nmr: number[]
  mass_spec: number[]
  mw: number | null
}

export interface ExampleMetadata {
  name: string
  filename: string
}

/**
 * Get list of available examples using Vite's import.meta.glob
 */
export async function getAvailableExamples(): Promise<ExampleMetadata[]> {
  // Use Vite's import.meta.glob to discover all JSON files in the examples directory
  const exampleModules = import.meta.glob<{ default: ExampleData }>('../data/examples/*.json', {
    eager: false,
  })

  const examples: ExampleMetadata[] = []

  for (const path in exampleModules) {
    // Extract filename from path (e.g., '../data/examples/kavaratamide-a.json' -> 'kavaratamide-a')
    const filename = path.split('/').pop()?.replace('.json', '') || ''
    
    // Skip template.json
    if (filename === 'template') {
      continue
    }

    // Load the module to get the name
    const loadModule = exampleModules[path]
    if (!loadModule) {
      continue
    }
    const module = await loadModule()
    const data = module.default || (module as any)
    examples.push({
      name: data.name || filename,
      filename: `${filename}.json`,
    })
  }

  // Sort by name
  examples.sort((a, b) => a.name.localeCompare(b.name))

  return examples
}

/**
 * Load a specific example by filename
 */
export async function loadExample(filename: string): Promise<ExampleData> {
  // Remove .json extension if present
  const nameWithoutExt = filename.replace('.json', '')
  
  try {
    const module = await import(`../data/examples/${nameWithoutExt}.json`)
    // Vite imports JSON as default export
    return (module.default || module) as ExampleData
  } catch (error) {
    throw new Error(`Failed to load example: ${filename}`)
  }
}

