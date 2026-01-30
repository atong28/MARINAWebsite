import { defineConfig, loadEnv } from 'vite'
import react from '@vitejs/plugin-react'

// https://vitejs.dev/config/
export default defineConfig(({ mode }) => {
  // Load env file from project root
  // Vite by default looks in the current working directory, but we want the parent (project root)
  const rootDir = process.cwd().endsWith('/frontend')
    ? process.cwd() + '/..'
    : process.cwd()
  const env = loadEnv(mode, rootDir, '')

  // Get ports from environment variables with defaults
  // These should be set in root .env file (single source of truth)
  const frontendPort = parseInt(env.FRONTEND_PORT || process.env.FRONTEND_PORT || '3000', 10)
  const backendPort = parseInt(env.BACKEND_PORT || process.env.BACKEND_PORT || '5000', 10)
  // In production, use relative path '/api' for nginx reverse proxy.
  // In development, can use full URL or relative path.
  const apiBase =
    env.VITE_API_BASE ||
    process.env.VITE_API_BASE ||
    (mode === 'production' ? '/api' : `http://localhost:${backendPort}/api`)

  // Optional per-model bases; fall back to apiBase if not provided.
  const apiBaseMarina = env.VITE_API_BASE_MARINA || process.env.VITE_API_BASE_MARINA || apiBase
  const apiBaseSpectre = env.VITE_API_BASE_SPECTRE || process.env.VITE_API_BASE_SPECTRE || apiBase
  
  return {
  plugins: [react()],
  server: {
      port: frontendPort,
    proxy: {
      '/api': {
          target: `http://localhost:${backendPort}`,
        changeOrigin: true,
        rewrite: (path) => path.replace(/^\/api/, '')
      }
    }
  },
    // Inject API base URLs into client code so they're available via import.meta.env.*
    // In production, these will typically be relative paths (e.g. '/api', '/api/marina', '/api/spectre')
    // and nginx will route to the appropriate backend.
    define: {
      'import.meta.env.VITE_API_BASE': JSON.stringify(apiBase),
      'import.meta.env.VITE_API_BASE_MARINA': JSON.stringify(apiBaseMarina),
      'import.meta.env.VITE_API_BASE_SPECTRE': JSON.stringify(apiBaseSpectre),
    },
  build: {
    outDir: 'dist',
    sourcemap: true
    }
  }
})

