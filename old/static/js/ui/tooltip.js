// Simple tooltip manager for image previews
(function(global){
  const tooltip = document.createElement('div');
  tooltip.className = 'bit-tooltip';
  tooltip.style.position = 'fixed';
  tooltip.style.display = 'none';
  tooltip.style.padding = '6px';
  tooltip.style.background = 'rgba(0,0,0,0.8)';
  tooltip.style.borderRadius = '6px';
  tooltip.style.boxShadow = '0 2px 8px rgba(0,0,0,0.3)';
  tooltip.style.zIndex = '9999';
  document.body.appendChild(tooltip);

  function showAt(x, y, imageDataUrl) {
    tooltip.textContent = '';
    const img = document.createElement('img');
    img.src = imageDataUrl;
    img.alt = 'Substructure';
    img.style.maxWidth = '220px';
    img.style.maxHeight = '220px';
    img.style.display = 'block';
    tooltip.appendChild(img);
    const pad = 12;
    tooltip.style.left = `${x + pad}px`;
    tooltip.style.top = `${y + pad}px`;
    tooltip.style.display = 'block';
  }

  function hide() {
    tooltip.style.display = 'none';
  }

  global.Tooltip = { showAt, hide };
})(window);


