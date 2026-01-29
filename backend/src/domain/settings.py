from typing import Literal, List, Dict
from pydantic.dataclasses import dataclass
from dataclasses import field

from src.marina.src.modules.spectre.args import SPECTREArgs

@dataclass
class MARINAArgs:
    # random seed
    seed: int = 0
    # path to load checkpoint from
    load_from_checkpoint: str | None = None
    # whether to do training
    train: bool = True
    # whether to do testing
    test: bool = True

    input_types: List[Literal['hsqc', 'c_nmr', 'h_nmr', 'mass_spec', 'mw']] = field(
        default_factory=lambda: ['hsqc', 'c_nmr', 'h_nmr', 'mass_spec', 'mw']
    )

    requires: List[Literal['hsqc', 'c_nmr', 'h_nmr', 'mass_spec', 'mw']] = field(
        default_factory=lambda: []
    )

    # training args
    debug: bool = False
    batch_size: int = 32
    num_workers: int = 4
    epochs: int = 750
    patience: int = 30
    persistent_workers: bool = True
    lr: float = 2e-4
    eta_min: float = 1e-5
    weight_decay: float = 0.0
    scheduler: Literal['cosine', 'none'] = 'cosine'
    freeze_weights: bool = False
    use_jaccard: bool = False
    warmup: bool = False
    accumulate_grad_batches_num: int = 4
    dropout: float = 0.1
    
    # jittering default value to wobble the spectra
    jittering: float = 1.0

    # BCE and cosine similarity loss lambda. 0 for full cosine similarity loss, 1 for full BCE loss.
    lambda_hybrid: float = 0.0
    
    # fp type for prediction and evaluation. fingerprint details should be stored in 
    #   DATASET_ROOT/RankingEntropy/
    # with the proper formatting.
    fp_type: Literal['RankingEntropy'] = 'RankingEntropy'
    
    # additional test types to be used for testing, always will test on all inputs
    additional_test_types: list[list[str]] = field(default_factory=lambda: [
        ['hsqc'], ['h_nmr'], ['c_nmr'], ['mass_spec']
    ])
    
    experiment_name: str = 'marina-development'
    project_name: str = 'MARINA'

    dim_model: int = 784
    nmr_dim_coords: List[int] = field(default_factory=lambda: [391, 391, 2])
    nmr_is_sign_encoding: List[bool] = field(default_factory=lambda: [False, False, True])
    c_nmr_dim_coords: List[int] = field(default_factory=lambda: [784])
    c_nmr_is_sign_encoding: List[bool] = field(default_factory=lambda: [False])
    h_nmr_dim_coords: List[int] = field(default_factory=lambda: [784])
    h_nmr_is_sign_encoding: List[bool] = field(default_factory=lambda: [False])
    ms_dim_coords: List[int] = field(default_factory=lambda: [392, 392])
    ms_is_sign_encoding: List[bool] = field(default_factory=lambda: [False, False])
    mw_dim_coords: List[int] = field(default_factory=lambda: [784])
    mw_is_sign_encoding: List[bool] = field(default_factory=lambda: [False])
    heads: int = 8
    layers: int = 16
    self_attn_layers: Dict[str, int] = field(default_factory=
        lambda: {'hsqc': 2, 'h_nmr': 1, 'c_nmr': 2, 'mass_spec': 1, 'mw': 1}
    )
    ff_dim: int = 3072
    out_dim: int = 16384

    c_wavelength_bounds: List[float] = field(
        default_factory=lambda: [0.01, 400.0])
    h_wavelength_bounds: List[float] = field(
        default_factory=lambda: [0.01, 20.0])
    mz_wavelength_bounds: List[float] = field(
        default_factory=lambda: [0.01, 5000.0])
    intensity_wavelength_bounds: List[float] = field(
        default_factory=lambda: [0.001, 200.0])
    mw_wavelength_bounds: List[float] = field(
        default_factory=lambda: [0.01, 7000.0])