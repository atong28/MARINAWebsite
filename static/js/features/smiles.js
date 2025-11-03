// SMILES search feature module
(function(global) {
  async function runSmilesSearch() {
    const smilesInput = document.getElementById('smiles-input');
    const smiles = smilesInput.value.trim();
    if (!smiles) {
      showMessage('Please enter a SMILES string', 'error');
      return;
    }

    const k = parseInt(document.getElementById('top-k').value);
    const resultsSection = document.getElementById('results-section');
    const loading = document.getElementById('loading');
    const resultsGrid = document.getElementById('results-grid');

    resultsSection.style.display = 'block';
    loading.style.display = 'flex';
    resultsGrid.innerHTML = '';

    try {
      const result = await ApiClient.postSmilesSearch({ smiles, k });
      if (result.error) {
        throw new Error(result.error);
      }
      loading.style.display = 'none';
      if (global.Results && global.Results.displayResults) {
        await global.Results.displayResults(result);
      }
    } catch (error) {
      loading.style.display = 'none';
      showMessage(error.message, 'error');
    }
  }

  global.SmilesSearch = { runSmilesSearch };
})(window);


