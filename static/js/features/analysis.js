// Analysis feature module
(function(global) {
  async function processAnalysisResults(result) {
    // Update similarity visualization
    await updateSimilarityVisualization(result.similarity_visualization);
    
    // Update change visualization
    await updateChangeVisualization(result.change_visualization);
    
    // Update fingerprint indices
    if (typeof updateFingerprintIndices === 'function') {
      updateFingerprintIndices(result.predicted_fp_indices, result.retrieved_molecule_fp_indices);
    }
    
    // Update fingerprint differences
    if (result.fingerprint_differences) {
      updateFingerprintDifferences(result.fingerprint_differences);
    }
    
    // Secondary retrieval rendering handled by caller (processAnalysisResults wrapper in app)
  }

  global.Analysis = { processAnalysisResults };
})(window);


