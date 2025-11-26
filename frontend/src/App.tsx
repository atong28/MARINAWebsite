import { BrowserRouter, Routes, Route } from 'react-router-dom'
import ErrorBoundary from './components/common/ErrorBoundary'
import MainPage from './pages/MainPage'
import AnalysisPage from './pages/AnalysisPage'
import { ROUTES } from './routes'
import './App.css'

function App() {
  return (
    <ErrorBoundary>
      <BrowserRouter>
        <Routes>
          <Route path={ROUTES.MAIN} element={<MainPage />} />
          <Route path={ROUTES.ANALYSIS} element={<AnalysisPage />} />
        </Routes>
      </BrowserRouter>
    </ErrorBoundary>
  )
}

export default App

