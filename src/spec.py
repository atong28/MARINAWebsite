def get_openapi_spec(base_url: str = "/") -> dict:
    return {
        "openapi": "3.0.3",
        "info": {
            "title": "MARINA Backend API",
            "version": "1.0.0"
        },
        "servers": [{"url": base_url.rstrip("/") or "/"}],
        "paths": {
            "/health": {
                "get": {
                    "summary": "Health status",
                    "responses": {"200": {"description": "OK"}}
                }
            },
            "/predict": {
                "post": {
                    "summary": "Predict and retrieve similar molecules",
                    "requestBody": {
                        "required": True,
                        "content": {
                            "application/json": {
                                "schema": {
                                    "type": "object",
                                    "properties": {
                                        "raw": {"type": "object"},
                                        "k": {"type": "integer"},
                                        "offset": {"type": "integer"},
                                        "limit": {"type": "integer"}
                                    },
                                    "required": ["raw"]
                                }
                            }
                        }
                    },
                    "responses": {"200": {"description": "Results returned"}}
                }
            },
            "/smiles-search": {
                "post": {
                    "summary": "Search by SMILES",
                    "requestBody": {
                        "required": True,
                        "content": {
                            "application/json": {
                                "schema": {
                                    "type": "object",
                                    "properties": {
                                        "smiles": {"type": "string"},
                                        "k": {"type": "integer"},
                                        "offset": {"type": "integer"},
                                        "limit": {"type": "integer"}
                                    },
                                    "required": ["smiles"]
                                }
                            }
                        }
                    },
                    "responses": {"200": {"description": "Results returned"}}
                }
            },
            "/analyze": {
                "post": {
                    "summary": "Analysis for selected result",
                    "requestBody": {
                        "required": True,
                        "content": {
                            "application/json": {
                                "schema": {
                                    "type": "object",
                                    "properties": {
                                        "current_data": {"type": "object"},
                                        "original_data": {"type": "object"},
                                        "target_smiles": {"type": "string"},
                                        "target_index": {"type": "integer"}
                                    },
                                    "required": ["current_data", "original_data", "target_smiles", "target_index"]
                                }
                            }
                        }
                    },
                    "responses": {"200": {"description": "Analysis returned"}}
                }
            },
            "/secondary-retrieval": {
                "post": {
                    "summary": "Secondary retrieval with difference fingerprint",
                    "requestBody": {
                        "required": True,
                        "content": {
                            "application/json": {
                                "schema": {
                                    "type": "object",
                                    "properties": {
                                        "predicted_fp": {"type": "array", "items": {"type": "number"}},
                                        "retrieved_fp": {"type": "array", "items": {"type": "number"}},
                                        "k": {"type": "integer"},
                                        "offset": {"type": "integer"},
                                        "limit": {"type": "integer"}
                                    },
                                    "required": ["predicted_fp", "retrieved_fp"]
                                }
                            }
                        }
                    },
                    "responses": {"200": {"description": "Results returned"}}
                }
            }
        }
    }


