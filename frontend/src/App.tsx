import { Suspense, lazy } from 'react'
import { BrowserRouter, Routes, Route } from 'react-router-dom'
import ErrorBoundary from './components/common/ErrorBoundary'
import ModelBootstrap from './components/common/ModelBootstrap'
import MainPage from './pages/MainPage'
import { ROUTES } from './routes'
import './App.css'

const AnalysisPage = lazy(() => import('./pages/AnalysisPage'))

function App() {
  return (
    <ErrorBoundary>
      <BrowserRouter>
        <ModelBootstrap />
        <Suspense fallback={<div style={{ padding: 20 }}>Loading…</div>}>
          <Routes>
            <Route path={ROUTES.MAIN} element={<MainPage />} />
            <Route path={ROUTES.ANALYSIS} element={<AnalysisPage />} />
          </Routes>
        </Suspense>
      </BrowserRouter>
    </ErrorBoundary>
  )
}

export default App

