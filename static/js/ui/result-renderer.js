// Result card renderer (DOM-safe)
(function(global) {
  function createCard(position, resultData) {
    const card = document.createElement('div');
    card.className = 'result-card';

    const similarity = (resultData && typeof resultData.similarity === 'number' && !isNaN(resultData.similarity))
      ? resultData.similarity : 0.0;
    const similarityPercent = (similarity * 100).toFixed(1);

    const smiles = (resultData && resultData.smiles) ? resultData.smiles : '';
    const svg = (resultData && resultData.svg) ? resultData.svg : '';
    const name = (resultData && resultData.name) ? resultData.name : '';
    const primaryLink = (resultData && resultData.primary_link) ? resultData.primary_link : '';
    const databaseLinks = (resultData && resultData.database_links) ? resultData.database_links : {};

    // Header
    const header = document.createElement('div');
    header.className = 'card-header';
    const title = document.createElement('div');
    title.className = 'card-title';
    title.textContent = `Result #${position}`;
    const score = document.createElement('div');
    score.className = 'similarity-score';
    score.textContent = `${similarityPercent}%`;
    header.appendChild(title);
    header.appendChild(score);
    card.appendChild(header);

    // Body
    const body = document.createElement('div');
    body.className = 'card-body';

    const moleculeContainer = document.createElement('div');
    moleculeContainer.id = `molecule-${position}`;
    moleculeContainer.className = 'molecule-container';
    if (svg) {
      if (svg.startsWith('data:image/')) {
        const img = document.createElement('img');
        img.className = 'molecule-image lazy-image';
        img.setAttribute('data-src', svg);
        img.alt = 'Molecular structure';
        img.style.maxWidth = '100%';
        img.style.maxHeight = '200px';
        img.style.borderRadius = '4px';
        moleculeContainer.appendChild(img);
        if (window.Lazy && Lazy.observeImage) Lazy.observeImage(img);
      } else if (svg.startsWith('<svg') || svg.startsWith('<')) {
        // Defer SVG injection until visible
        moleculeContainer.classList.add('lazy-svg');
        moleculeContainer.setAttribute('data-svg', svg);
        if (window.Lazy && Lazy.observeSvg) Lazy.observeSvg(moleculeContainer);
      }
    } else {
      const noMol = document.createElement('div');
      noMol.className = 'no-molecule';
      noMol.textContent = 'No structure available';
      moleculeContainer.appendChild(noMol);
    }
    body.appendChild(moleculeContainer);

    if (name) {
      const nameDiv = document.createElement('div');
      nameDiv.className = 'molecule-name';
      if (primaryLink) {
        const a = document.createElement('a');
        a.href = primaryLink;
        a.target = '_blank';
        a.className = 'primary-link';
        a.textContent = name;
        nameDiv.appendChild(a);
      } else {
        nameDiv.textContent = name;
      }
      body.appendChild(nameDiv);
    }

    // Database links
    if (databaseLinks && Object.keys(databaseLinks).length > 0) {
      const linksDiv = document.createElement('div');
      linksDiv.className = 'database-links';
      if (databaseLinks.coconut) {
        const a = document.createElement('a');
        a.href = databaseLinks.coconut;
        a.textContent = 'COCONUT';
        a.target = '_blank';
        a.className = 'db-link coconut';
        linksDiv.appendChild(a);
      }
      if (databaseLinks.lotus) {
        const a = document.createElement('a');
        a.href = databaseLinks.lotus;
        a.textContent = 'LOTUS';
        a.target = '_blank';
        a.className = 'db-link lotus';
        linksDiv.appendChild(a);
      }
      if (databaseLinks.npmrd) {
        const a = document.createElement('a');
        a.href = databaseLinks.npmrd;
        a.textContent = 'NP-MRD';
        a.target = '_blank';
        a.className = 'db-link npmrd';
        linksDiv.appendChild(a);
      }
      body.appendChild(linksDiv);
    }

    // SMILES label
    const smilesLabel = document.createElement('div');
    const bold = document.createElement('strong');
    bold.textContent = 'SMILES:';
    smilesLabel.appendChild(bold);
    body.appendChild(smilesLabel);

    const smilesCode = document.createElement('div');
    smilesCode.className = 'smiles-code';
    smilesCode.textContent = smiles;
    body.appendChild(smilesCode);

    // Similarity row
    const simRow = document.createElement('div');
    const simLabel = document.createElement('strong');
    simLabel.textContent = 'Similarity:';
    simRow.appendChild(simLabel);
    const simVal = document.createElement('span');
    simVal.textContent = ` ${similarityPercent}%`;
    simRow.appendChild(simVal);
    body.appendChild(simRow);

    // Actions
    const actions = document.createElement('div');
    actions.className = 'card-actions';
    const analyzeBtn = document.createElement('button');
    analyzeBtn.className = 'btn btn-primary btn-sm';
    const icon = document.createElement('i');
    icon.className = 'fas fa-microscope';
    analyzeBtn.appendChild(icon);
    analyzeBtn.appendChild(document.createTextNode(' Analyze'));
    analyzeBtn.addEventListener('click', function() { openAnalysis(position - 1); });
    actions.appendChild(analyzeBtn);
    body.appendChild(actions);

    card.appendChild(body);
    return card;
  }

  global.ResultRenderer = { createCard };
})(window);


