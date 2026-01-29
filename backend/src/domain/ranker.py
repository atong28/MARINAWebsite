"""
Domain ranker module - imports RankingSet from marina and extends it.
"""
from typing import List, Tuple
import torch

# Import RankingSet directly from marina
from src.marina.src.modules.core.ranker import RankingSet as _RankingSet

# Extend RankingSet with domain-specific method
class RankingSet(_RankingSet):
    """
    Extended RankingSet with domain-specific retrieve_with_scores method.
    """
    def retrieve_with_scores(
        self, query: torch.Tensor, n: int = 50
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Public helper returning both similarity scores and indices.

        Returns:
            sims_topk: (n,)      if single query; else (n, B)
            idxs_topk: (n,) int  if single query; else (n, B)
        """
        if query.dim() == 1:
            query = query.unsqueeze(0)
        sims = self._sims(query)  # (N, B)
        k_eff = min(n, sims.size(0))
        sims_topk, idxs_topk = torch.topk(sims, k=k_eff, dim=0)
        if sims_topk.size(1) == 1:
            return sims_topk.squeeze(1), idxs_topk.squeeze(1)
        return sims_topk, idxs_topk

# Domain-specific function for filtering rankingset rows
def filter_rankingset_rows(store: torch.Tensor, kept_indices: List[int]) -> torch.Tensor:
    """
    Build a new CSR tensor containing only the rows in kept_indices, in that order.

    This is used to prefilter the retrieval rankingset (e.g., by molecular weight)
    before running similarity search. Works only for sparse_csr tensors; for dense
    tensors, simple row indexing is sufficient.
    """
    if not kept_indices:
        return torch.sparse_csr_tensor(
            torch.tensor([0, 0], dtype=torch.int64),
            torch.tensor([], dtype=torch.int64),
            torch.tensor([], dtype=torch.float32),
            size=(0, store.shape[1]),
        )

    if store.layout != torch.sparse_csr:
        # Dense fallback: simple row indexing
        return store[kept_indices]

    crow = store.crow_indices()
    col = store.col_indices()
    val = store.values()

    new_crow = [0]
    new_cols = []
    new_vals = []
    nnz_so_far = 0

    for idx in kept_indices:
        row_start = int(crow[idx])
        row_end = int(crow[idx + 1])
        if row_start >= row_end:
            new_crow.append(nnz_so_far)
            continue
        row_cols = col[row_start:row_end]
        row_vals = val[row_start:row_end]
        new_cols.append(row_cols)
        new_vals.append(row_vals)
        nnz_so_far += (row_end - row_start)
        new_crow.append(nnz_so_far)

    if nnz_so_far == 0:
        return torch.sparse_csr_tensor(
            torch.tensor([0, 0], dtype=torch.int64),
            torch.tensor([], dtype=torch.int64),
            torch.tensor([], dtype=torch.float32),
            size=(len(kept_indices), store.shape[1]),
        )

    new_crow_t = torch.tensor(new_crow, dtype=torch.int64)
    new_cols_t = torch.cat(new_cols).to(torch.int64)
    new_vals_t = torch.cat(new_vals).to(torch.float32)

    return torch.sparse_csr_tensor(
        new_crow_t,
        new_cols_t,
        new_vals_t,
        size=(len(kept_indices), store.shape[1]),
    )

__all__ = ['RankingSet', 'filter_rankingset_rows']
