"""TT-rounding — rank reduction of Tensor Trains.

Rounding reduces the TT-ranks while controlling the approximation error.
The standard algorithm (Oseledets, 2011):

1. Bring the TT to right-canonical form (left-to-right QR sweep or
   orthogonalize to site ``d-1``).
2. Right-to-left SVD sweep with truncation of singular values.

The truncation threshold at each step is
``delta = eps / sqrt(d - 1) * ||TT||_F``,
which guarantees the overall relative error is at most *eps*.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.linalg import svd

if TYPE_CHECKING:
    from collections.abc import Sequence

    from tensortrain._core import TensorTrain


def tt_round(
    tt: TensorTrain,
    eps: float | None = None,
    max_ranks: Sequence[int] | None = None,
) -> TensorTrain:
    """Round a TT by reducing ranks via left-QR / right-SVD sweeps.

    This is the standard TT-rounding algorithm (Algorithm 2 in [1]).

    Parameters
    ----------
    tt : TensorTrain
        Input tensor train.
    eps : float, optional
        Relative Frobenius-norm tolerance.
    max_ranks : sequence of int, optional
        Maximum allowed ranks ``(r_1, ..., r_{d-1})``.  Applied together
        with *eps* when both are given.

    Returns
    -------
    TensorTrain
        A new, rank-reduced TensorTrain.

    Raises
    ------
    ValueError
        If neither *eps* nor *max_ranks* is provided.

    References
    ----------
    .. [1] I. V. Oseledets, "Tensor-Train Decomposition", *SIAM J. Sci.
       Comput.*, 33(5), 2295-2317, 2011.
    """
    from tensortrain._core import TensorTrain

    if eps is None and max_ranks is None:
        raise ValueError("At least one of eps or max_ranks must be provided.")

    d = tt.ndim
    if max_ranks is not None:
        max_ranks = list(max_ranks)
        if len(max_ranks) != d - 1:
            raise ValueError(
                f"max_ranks must have {d - 1} elements, got {len(max_ranks)}."
            )

    cores = [c.copy() for c in tt._cores]

    # --- Step 1: Left-to-right QR sweep (left-orthogonalize) ---
    for k in range(d - 1):
        r_left, n_k, r_right = cores[k].shape
        Q, R = np.linalg.qr(cores[k].reshape(r_left * n_k, r_right))
        r_new = Q.shape[1]
        cores[k] = Q.reshape(r_left, n_k, r_new)
        # Absorb R into next core
        cores[k + 1] = np.einsum("ij,jkl->ikl", R, cores[k + 1])

    # --- Compute delta from the norm (now concentrated in the last core) ---
    if eps is not None:
        tt_norm = np.linalg.norm(cores[-1].ravel())
        delta = (eps / np.sqrt(max(d - 1, 1))) * tt_norm
    else:
        delta = None

    # --- Step 2: Right-to-left SVD sweep with truncation ---
    for k in range(d - 1, 0, -1):
        r_left, n_k, r_right = cores[k].shape
        # Reshape to (r_left, n_k * r_right)
        mat = cores[k].reshape(r_left, n_k * r_right)

        U, sigma, Vt = svd(mat, full_matrices=False)

        # Determine truncation rank
        r = len(sigma)

        if delta is not None:
            r = _truncation_rank(sigma, delta)

        if max_ranks is not None:
            r = min(r, max_ranks[k - 1])

        r = max(r, 1)

        # Truncate
        Vt_trunc = Vt[:r, :]
        cores[k] = Vt_trunc.reshape(r, n_k, r_right)

        # Absorb U * S into the previous core
        US = U[:, :r] * sigma[:r]  # (r_left, r)
        cores[k - 1] = np.einsum("ijk,kl->ijl", cores[k - 1], US)

    return TensorTrain(cores)


def tt_round_rl(tt: TensorTrain, eps: float) -> TensorTrain:
    """Round a TT using right-to-left SVD sweep only.

    First orthogonalizes to site ``d-1`` (left-canonical form via QR),
    then performs a right-to-left SVD sweep with epsilon truncation.
    Equivalent to ``roundTT2`` in the MATLAB code.

    Parameters
    ----------
    tt : TensorTrain
        Input tensor train.
    eps : float
        Relative Frobenius-norm tolerance.

    Returns
    -------
    TensorTrain
    """
    # This is effectively the same algorithm as tt_round with eps only
    return tt_round(tt, eps=eps)


def tt_round_to_ranks(tt: TensorTrain, ranks: Sequence[int]) -> TensorTrain:
    """Round a TT to exact target ranks.

    Parameters
    ----------
    tt : TensorTrain
        Input tensor train.
    ranks : sequence of int
        Target ranks ``(r_1, ..., r_{d-1})``.

    Returns
    -------
    TensorTrain
    """
    return tt_round(tt, max_ranks=ranks)


def _truncation_rank(sigma: np.ndarray, delta: float) -> int:
    """Find smallest rank *r* such that the tail error <= delta.

    Parameters
    ----------
    sigma : 1-D ndarray
        Singular values in descending order.
    delta : float
        Absolute error threshold.

    Returns
    -------
    int
    """
    tail_sq = np.cumsum(sigma[::-1] ** 2)[::-1]
    delta_sq = delta ** 2
    for r in range(1, len(sigma)):
        if tail_sq[r] <= delta_sq:
            return r
    return len(sigma)
