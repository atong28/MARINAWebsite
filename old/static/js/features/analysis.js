// Analysis feature module
(function(global) {
  async function processAnalysisResults(result) {
    // Update similarity visualization
    await updateSimilarityVisualization(result.similarity_visualization);
    
    // Update change visualization
    await updateChangeVisualization(result.change_visualization);
    
    // Store analysis result separately - DO NOT overwrite currentAnalysisResult
    // The original currentAnalysisResult contains bit_environments needed for highlighting
    // Analysis result only contains visualization data and fingerprint differences
    try {
      if (window.State && State.set) {
        // Store analysis result separately for visualization functions
        State.set('analysisResult', result);
        // DO NOT overwrite currentAnalysisResult - it contains bit_environments from original predict response
      }
    } catch (e) {}

    // DO NOT call updateFingerprintIndices here
    // Fingerprint indices are set once when opening analysis (from original result)
    // Analysis results should only update visualizations and fingerprint differences
    
    // Update fingerprint differences
    if (result.fingerprint_differences) {
      updateFingerprintDifferences(result.fingerprint_differences);
    }
    
    // Secondary retrieval rendering handled by caller (processAnalysisResults wrapper in app)
  }

  global.Analysis = { processAnalysisResults };
})(window);


