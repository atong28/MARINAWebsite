// Results feature module - handles displaying prediction or search results
(function(global) {
  async function displayResults(result) {
    const resultsGrid = document.getElementById('results-grid');
    if (!resultsGrid) return;

    if (result.results && Array.isArray(result.results)) {
      const results = result.results;
      if (results.length === 0) {
        resultsGrid.textContent = '';
        const p = document.createElement('p');
        p.className = 'error-message';
        p.textContent = 'No results returned from the model.';
        resultsGrid.appendChild(p);
        return;
      }
      // Store globally for analysis workflow
      if (typeof window !== 'undefined') {
        if (window.State && State.set) {
          State.set('currentResults', results);
        } else {
          window.currentResults = results;
        }
      }
      // Persist main page state after results are displayed
      try {
        if (window.State && State.persistMain) {
          const inputs = (typeof collectInputData === 'function') ? collectInputData() : {};
          const smilesInput = (document.getElementById('smiles-input') || {}).value || '';
          State.persistMain({
            inputs,
            smilesInput,
            currentResults: results
          });
        }
      } catch (e) {}
      resultsGrid.textContent = '';
      // Virtualize if many results
      if (results.length > 50) {
        mountVirtualizedResults(resultsGrid, results);
      } else {
        const maxResults = Math.min(10, results.length);
        for (let i = 0; i < maxResults; i++) {
          const resultData = results[i];
          const card = (global.ResultRenderer && global.ResultRenderer.createCard)
            ? global.ResultRenderer.createCard(i + 1, resultData, { showAnalyze: true })
            : createFallbackCard(i + 1, resultData);
          resultsGrid.appendChild(card);
        }
      }
      return;
    }

    // Fallback for old response format
    const indices = result.topk_indices || result.topk || [];
    const scores = result.topk_scores || [];
    if (!indices || indices.length === 0) {
      resultsGrid.textContent = '';
      const p = document.createElement('p');
      p.className = 'error-message';
      p.textContent = 'No results returned from the model.';
      resultsGrid.appendChild(p);
      return;
    }

    // Create placeholder cards (kept for backward compatibility)
    resultsGrid.textContent = '';
    const maxResults = Math.min(10, indices.length);
    for (let i = 0; i < maxResults; i++) {
      const card = createBlankResultCard(i + 1);
      resultsGrid.appendChild(card);
    }

    // No /meta endpoint; render minimal placeholders only
    for (let i = 0; i < maxResults; i++) {
      updateResultCard(i, null, '', scores[i] || 0.0, null);
    }
  }

  // Simple incremental rendering on scroll for large result sets
  function mountVirtualizedResults(container, results) {
    const batchSize = 20;
    let nextIndex = 0;

    function renderBatch() {
      const end = Math.min(nextIndex + batchSize, results.length);
      for (let i = nextIndex; i < end; i++) {
        const resultData = results[i];
        const card = (global.ResultRenderer && global.ResultRenderer.createCard)
          ? global.ResultRenderer.createCard(i + 1, resultData, { showAnalyze: true })
          : createFallbackCard(i + 1, resultData);
        container.appendChild(card);
      }
      nextIndex = end;
    }

    function onScroll() {
      const threshold = 300; // px from bottom
      const atBottom = container.scrollHeight - container.scrollTop - container.clientHeight < threshold;
      if (atBottom && nextIndex < results.length) {
        renderBatch();
      }
    }

    // Prepare container for scrolling if needed
    container.style.maxHeight = '1200px';
    container.style.overflowY = 'auto';
    container.addEventListener('scroll', onScroll);
    renderBatch();
  }

  // Legacy helpers extracted for backward compatibility
  function createBlankResultCard(position) {
    const card = document.createElement('div');
    card.className = 'result-card loading';

    const header = document.createElement('div');
    header.className = 'card-header';
    const title = document.createElement('div');
    title.className = 'card-title';
    title.textContent = `Result #${position}`;
    const score = document.createElement('div');
    score.className = 'similarity-score loading';
    score.textContent = 'Loading...';
    header.appendChild(title);
    header.appendChild(score);

    const body = document.createElement('div');
    body.className = 'card-body';
    const mol = document.createElement('div');
    mol.id = `molecule-${position}`;
    mol.className = 'molecule-container loading';
    const spinner = document.createElement('div');
    spinner.className = 'loading-spinner';
    mol.appendChild(spinner);
    body.appendChild(mol);

    const smilesLabel = document.createElement('div');
    const strong = document.createElement('strong');
    strong.textContent = 'SMILES:';
    smilesLabel.appendChild(strong);
    body.appendChild(smilesLabel);

    const smilesCode = document.createElement('div');
    smilesCode.className = 'smiles-code loading';
    smilesCode.textContent = 'Loading...';
    body.appendChild(smilesCode);

    const simRow = document.createElement('div');
    const simLabel = document.createElement('strong');
    simLabel.textContent = 'Similarity:';
    simRow.appendChild(simLabel);
    const simVal = document.createElement('span');
    simVal.className = 'loading';
    simVal.textContent = 'Loading...';
    simRow.appendChild(simVal);
    body.appendChild(simRow);

    card.appendChild(header);
    card.appendChild(body);
    return card;
  }

  function updateResultCard(position, idx, smiles, tanimoto, errorMessage) {
    const card = document.querySelector(`#results-grid .result-card:nth-child(${position + 1})`);
    if (!card) return;
    card.className = 'result-card';
    card.textContent = '';

    const header = document.createElement('div');
    header.className = 'card-header';
    const title = document.createElement('div');
    title.className = 'card-title';
    title.textContent = `Result #${position + 1}`;
    const score = document.createElement('div');
    score.className = errorMessage ? 'similarity-score error' : 'similarity-score';
    score.textContent = errorMessage ? 'Error' : `${(tanimoto * 100).toFixed(1)}%`;
    header.appendChild(title);
    header.appendChild(score);
    card.appendChild(header);

    const body = document.createElement('div');
    body.className = 'card-body';

    const mol = document.createElement('div');
    mol.id = `molecule-${position + 1}`;
    mol.className = errorMessage ? 'molecule-container error' : 'molecule-container';

    if (errorMessage) {
      const icon = document.createElement('i');
      icon.className = 'fas fa-exclamation-triangle';
      mol.appendChild(icon);
      body.appendChild(mol);

      const errLabel = document.createElement('div');
      const strong = document.createElement('strong');
      strong.textContent = 'Error:';
      errLabel.appendChild(strong);
      body.appendChild(errLabel);

      const errMsg = document.createElement('div');
      errMsg.className = 'error-message';
      errMsg.textContent = errorMessage;
      body.appendChild(errMsg);
    } else {
      const noMol = document.createElement('div');
      noMol.className = 'no-molecule';
      noMol.textContent = 'Loading enhanced visualization...';
      mol.appendChild(noMol);
      body.appendChild(mol);

      const smilesLabel = document.createElement('div');
      const strong = document.createElement('strong');
      strong.textContent = 'SMILES:';
      smilesLabel.appendChild(strong);
      body.appendChild(smilesLabel);

      const smilesCode = document.createElement('div');
      smilesCode.className = 'smiles-code';
      smilesCode.textContent = smiles || '';
      body.appendChild(smilesCode);

      const simRow = document.createElement('div');
      const simLabel = document.createElement('strong');
      simLabel.textContent = 'Similarity:';
      simRow.appendChild(simLabel);
      const simVal = document.createElement('span');
      simVal.textContent = ` ${(tanimoto || 0).toFixed(3)}`;
      simRow.appendChild(simVal);
      body.appendChild(simRow);
    }

    card.appendChild(body);
  }

  function createFallbackCard(position, resultData) {
    const card = document.createElement('div');
    card.className = 'result-card';
    const header = document.createElement('div');
    header.className = 'card-header';
    const title = document.createElement('div');
    title.className = 'card-title';
    title.textContent = `Result #${position}`;
    header.appendChild(title);
    card.appendChild(header);
    return card;
  }

  global.Results = { displayResults, createBlankResultCard, updateResultCard };
})(window);


