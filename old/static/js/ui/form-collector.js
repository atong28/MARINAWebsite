// Form collector wrapper module
(function(global) {
  function collect() {
    if (typeof collectInputData === 'function') {
      return collectInputData();
    }
    return {};
  }

  function validate(data) {
    if (typeof validateInputData === 'function') {
      return validateInputData(data);
    }
    return [];
  }

  global.FormCollector = { collect, validate };
})(window);


