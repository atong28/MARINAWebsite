// Prediction feature wrapper module
(function(global) {
  async function run() {
    const data = (global.FormCollector && global.FormCollector.collect) ? global.FormCollector.collect() : (typeof collectInputData === 'function' ? collectInputData() : {});
    const errors = (global.FormCollector && global.FormCollector.validate) ? global.FormCollector.validate(data) : (typeof validateInputData === 'function' ? validateInputData(data) : []);
    if (errors.length > 0) {
      errors.forEach(error => showMessage(error, 'error'));
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
      const healthData = await ApiClient.getHealth();
      if (!healthData.model_loaded) {
        if (healthData.error) {
          throw new Error(`Model initialization failed: ${healthData.error}`);
        } else {
          throw new Error('Model is still initializing. Please wait a moment and try again.');
        }
      }

      const requestBody = { raw: data, k };
      const result = await ApiClient.postPredict(requestBody);
      if (result.error) throw new Error(result.error);

      if (result.pred_fp) {
        if (global.State && State.set) {
          State.set('currentPredictedFp', result.pred_fp);
        } else {
          global.currentPredictedFp = result.pred_fp;
        }
      }

      loading.style.display = 'none';
      if (global.Results && global.Results.displayResults) {
        await global.Results.displayResults(result);
      }
    } catch (error) {
      loading.style.display = 'none';
      let userMessage = error.message;
      if (error.name === 'TypeError' && error.message.includes('fetch')) {
        userMessage = 'Cannot connect to the prediction server. Please ensure the backend is running and accessible.';
      } else if (error.message.includes('timeout')) {
        userMessage = 'Request timed out. The server may be overloaded. Please try again.';
      } else if (error.message.includes('NetworkError')) {
        userMessage = 'Network error. Please check your internet connection and try again.';
      }
      showDetailedError(userMessage, error);
      console.error('Prediction error:', error);
    }
  }

  async function display(result) {
    if (global.Results && global.Results.displayResults) {
      return await global.Results.displayResults(result);
    }
  }

  global.Prediction = { run, display };
})(window);


