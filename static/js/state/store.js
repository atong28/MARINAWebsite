// Simple centralized state store
(function(global) {
  const state = Object.create(null);

  function set(key, value) {
    state[key] = value;
    return value;
  }

  function get(key, defaultValue = undefined) {
    return Object.prototype.hasOwnProperty.call(state, key) ? state[key] : defaultValue;
  }

  function reset(keys) {
    if (!keys) {
      for (const k of Object.keys(state)) delete state[k];
      return;
    }
    keys.forEach(k => { delete state[k]; });
  }

  global.State = { set, get, reset };
})(window);


