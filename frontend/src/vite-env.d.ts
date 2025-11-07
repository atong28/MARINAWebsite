/// <reference types="vite/client" />

interface ImportMeta {
  glob<T = any>(
    pattern: string,
    options?: {
      eager?: boolean
      import?: string
      query?: string | Record<string, string | boolean>
    }
  ): Record<string, () => Promise<{ default: T }>>
}

