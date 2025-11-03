// UI toggle helpers for collapsible sections
(function(global) {
  function toggleSection(sectionId) {
    const content = document.getElementById(`${sectionId}-content`);
    const icon = document.getElementById(`${sectionId}-icon`);
    if (content && icon) {
      const collapsed = content.classList.contains('collapsed');
      if (collapsed) {
        content.classList.remove('collapsed');
        icon.classList.remove('rotated');
      } else {
        content.classList.add('collapsed');
        icon.classList.add('rotated');
      }
    }
  }

  function toggleSectionAnalysis(sectionId) {
    const content = document.getElementById(`analysis-${sectionId}-content`);
    const icon = document.getElementById(`analysis-${sectionId}-icon`);
    if (content && icon) {
      const collapsed = content.classList.contains('collapsed');
      if (collapsed) {
        content.classList.remove('collapsed');
        icon.classList.remove('rotated');
      } else {
        content.classList.add('collapsed');
        icon.classList.add('rotated');
      }
    }
  }

  global.Toggles = { toggleSection, toggleSectionAnalysis };
})(window);


