import json
import logging
import threading
from dataclasses import dataclass, field
from typing import Optional

import torch

from src.config import CKPT_PATH, PARAMS_PATH, RETRIEVAL_PATH
from src.domain.settings import MARINAArgs
from src.domain.model import MARINA
from src.domain.fingerprint.fp_loader import FPLoader, make_fp_loader
from src.domain.ranker import RankingSet


logger = logging.getLogger(__name__)


@dataclass
class ModelSession:
    """
    Encapsulates the loaded MARINA model, fingerprint loader, and rankingset store.

    This groups together what used to be module-level globals in predictor.py into a
    single, explicitly managed object.
    """

    args: MARINAArgs
    model: MARINA
    fp_loader: FPLoader
    fp_type: str
    ckpt_path: str
    params_path: str
    retrieval_path: str

    _rankingset_store: Optional[torch.Tensor] = None
    _lock: threading.Lock = field(default_factory=threading.Lock, repr=False)

    @classmethod
    def from_paths(
        cls,
        ckpt_path: str = CKPT_PATH,
        params_path: str = PARAMS_PATH,
        retrieval_path: str = RETRIEVAL_PATH,
    ) -> "ModelSession":
        """
        Factory that loads MARINAArgs, constructs the fingerprint loader and model,
        and loads the checkpoint weights.
        """
        logger.info("Loading model params from %s", params_path)
        with open(params_path, "r") as f:
            params = json.load(f)
        args = MARINAArgs(**params)

        # create fp loader according to args
        fp_loader = make_fp_loader(
            args.fp_type,
            entropy_out_dim=args.out_dim,
            retrieval_path=retrieval_path,
        )

        logger.info("Instantiating MARINA model")
        model = MARINA(args, fp_loader)

        logger.info("Loading checkpoint from %s (cpu)", ckpt_path)
        ckpt = torch.load(ckpt_path, map_location="cpu")
        sd = ckpt.get("state_dict", ckpt)

        # Defer state-dict key normalization to MARINA / caller; we simply load here.
        res = model.load_state_dict(sd, strict=False)
        logger.info("Model load_state_dict result: %s", res)

        model.to(torch.device("cpu"))
        model.eval()
        torch.set_grad_enabled(False)

        return cls(
            args=args,
            model=model,
            fp_loader=fp_loader,
            fp_type=args.fp_type,
            ckpt_path=ckpt_path,
            params_path=params_path,
            retrieval_path=retrieval_path,
        )

    def get_rankingset(self) -> RankingSet:
        """
        Lazily load the retrieval rankingset tensor and wrap it in a RankingSet.
        """
        with self._lock:
            if self._rankingset_store is None:
                store = self.fp_loader.load_rankingset(self.fp_type)
                self._rankingset_store = store
        return RankingSet(store=self._rankingset_store, metric="cosine")


