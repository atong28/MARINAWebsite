// Lazy loader for images and SVG strings injected on visibility
(function(global) {
  const observers = {};

  function ensureImageObserver() {
    if (observers.image) return observers.image;
    observers.image = new IntersectionObserver((entries, obs) => {
      entries.forEach(entry => {
        if (entry.isIntersecting) {
          const img = entry.target;
          const src = img.getAttribute('data-src');
          if (src) {
            img.src = src;
            img.removeAttribute('data-src');
            img.classList.remove('lazy-image');
          }
          obs.unobserve(img);
        }
      });
    }, { rootMargin: '200px 0px' });
    return observers.image;
  }

  function ensureSvgObserver() {
    if (observers.svg) return observers.svg;
    observers.svg = new IntersectionObserver((entries, obs) => {
      entries.forEach(entry => {
        if (entry.isIntersecting) {
          const holder = entry.target;
          const svg = holder.getAttribute('data-svg');
          if (svg) {
            holder.innerHTML = svg;
            holder.removeAttribute('data-svg');
            holder.classList.remove('lazy-svg');
          }
          obs.unobserve(holder);
        }
      });
    }, { rootMargin: '200px 0px' });
    return observers.svg;
  }

  function observeImage(img) {
    const imageObserver = ensureImageObserver();
    imageObserver.observe(img);
  }

  function observeSvg(holder) {
    const svgObserver = ensureSvgObserver();
    svgObserver.observe(holder);
  }

  function autoRegister(container = document) {
    container.querySelectorAll('img.lazy-image[data-src]').forEach(observeImage);
    container.querySelectorAll('.lazy-svg[data-svg]').forEach(observeSvg);
  }

  document.addEventListener('DOMContentLoaded', function() {
    autoRegister(document);
  });

  global.Lazy = { observeImage, observeSvg, autoRegister };
})(window);


