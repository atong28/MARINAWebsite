// Simple centralized state store with localStorage persistence
(function(global) {
  const state = Object.create(null);
  const STORAGE_KEY = 'marinaStateV1';

  function _loadFromStorage() {
    try {
      const stored = localStorage.getItem(STORAGE_KEY);
      if (stored) {
        const parsed = JSON.parse(stored);
        Object.assign(state, parsed);
      }
    } catch (e) {
      // Ignore storage errors (e.g., private browsing)
    }
  }

  function _saveToStorage() {
    try {
      localStorage.setItem(STORAGE_KEY, JSON.stringify(state));
    } catch (e) {
      // Ignore storage errors
    }
  }

  // Load persisted state on init
  _loadFromStorage();

  function set(key, value) {
    state[key] = value;
    _saveToStorage();
    return value;
  }

  function get(key, defaultValue = undefined) {
    return Object.prototype.hasOwnProperty.call(state, key) ? state[key] : defaultValue;
  }

  function remove(key) {
    if (Object.prototype.hasOwnProperty.call(state, key)) {
      delete state[key];
      _saveToStorage();
    }
  }

  function reset(keys) {
    if (!keys) {
      for (const k of Object.keys(state)) delete state[k];
      _saveToStorage();
      return;
    }
    keys.forEach(k => { delete state[k]; });
    _saveToStorage();
  }

  // Persist/restore Main page state
  function persistMain(snapshot) {
    if (!snapshot) return;
    try {
      state.main = snapshot;
      _saveToStorage();
    } catch (e) {}
  }

  function restoreMain() {
    return state.main || null;
  }

  // Persist/restore Analysis page state
  function persistAnalysis(snapshot) {
    if (!snapshot) return;
    try {
      state.analysis = snapshot;
      _saveToStorage();
    } catch (e) {}
  }

  function restoreAnalysis() {
    return state.analysis || null;
  }

  global.State = { set, get, remove, reset, persistMain, restoreMain, persistAnalysis, restoreAnalysis };
})(window);


