"""Decomposition of dense tensors into Tensor Train format.

The core algorithm is TT-SVD (Oseledets, 2011), which sweeps left-to-right
through the modes, performing an SVD at each step and optionally truncating
singular values to control the approximation error.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.linalg import svd

if TYPE_CHECKING:
    from collections.abc import Sequence

    from tensortrain._core import TensorTrain


def tensor_to_tt(
    tensor: np.ndarray,
    eps: float | None = None,
    max_ranks: Sequence[int] | None = None,
) -> TensorTrain:
    """Decompose a dense tensor into Tensor Train format via TT-SVD.

    Performs a left-to-right sweep of SVDs. At each step the singular values
    are truncated either by an absolute tolerance derived from *eps*, by
    explicit rank bounds, or not at all (exact decomposition).

    Parameters
    ----------
    tensor : ndarray
        A *d*-dimensional array (``d >= 2``).
    eps : float, optional
        Relative Frobenius-norm tolerance.  The resulting TT satisfies
        ``||X - X_tt||_F <= eps * ||X||_F``.  Internally the per-step
        threshold is ``delta = eps / sqrt(d - 1) * ||X||_F`` (Oseledets,
        2011, Theorem 2.2).
    max_ranks : sequence of int, optional
        Maximum allowed ranks ``(r_1, ..., r_{d-1})``.  When given together
        with *eps*, both constraints are applied.

    Returns
    -------
    TensorTrain
        The TT decomposition.

    Raises
    ------
    ValueError
        If *tensor* has fewer than 2 dimensions, or *max_ranks* has wrong
        length.

    Notes
    -----
    If neither *eps* nor *max_ranks* is given, the decomposition is exact
    (no truncation), which can produce large ranks.

    References
    ----------
    .. [1] I. V. Oseledets, "Tensor-Train Decomposition", *SIAM J. Sci.
       Comput.*, 33(5), 2295-2317, 2011.

    Examples
    --------
    >>> import numpy as np
    >>> from tensortrain import TensorTrain
    >>> X = np.random.randn(4, 5, 3)
    >>> tt = TensorTrain.from_tensor(X, eps=1e-12)
    >>> np.allclose(tt.full(), X)
    True
    """
    from tensortrain._core import TensorTrain

    tensor = np.asarray(tensor, dtype=np.float64)
    d = tensor.ndim
    if d < 2:
        raise ValueError(f"Tensor must have at least 2 dimensions, got {d}.")

    shape = tensor.shape

    # Validate max_ranks
    if max_ranks is not None:
        max_ranks = list(max_ranks)
        if len(max_ranks) != d - 1:
            raise ValueError(
                f"max_ranks must have {d - 1} elements, got {len(max_ranks)}."
            )

    # Compute the per-step truncation threshold
    if eps is not None:
        tensor_norm = np.linalg.norm(tensor.ravel())
        delta = (eps / np.sqrt(d - 1)) * tensor_norm
    else:
        delta = None

    cores: list[np.ndarray] = []
    C = tensor.copy()
    r_left = 1  # current left rank

    for k in range(d - 1):
        n_k = shape[k]
        # Reshape to 2-D: (r_left * n_k) x (remaining elements)
        C = C.reshape(r_left * n_k, -1)

        # Economy SVD
        U, sigma, Vt = svd(C, full_matrices=False)

        # Determine truncation rank
        r = len(sigma)

        if delta is not None:
            r = _truncation_rank_eps(sigma, delta)

        if max_ranks is not None:
            r = min(r, max_ranks[k])

        # Ensure at least rank 1
        r = max(r, 1)

        # Truncate
        U_trunc = U[:, :r]
        sigma_trunc = sigma[:r]
        Vt_trunc = Vt[:r, :]

        # Core k: reshape U_trunc to (r_left, n_k, r)
        cores.append(U_trunc.reshape(r_left, n_k, r))

        # Prepare residual for next iteration: S * Vt
        C = np.diag(sigma_trunc) @ Vt_trunc

        r_left = r

    # Last core: whatever remains, shape (r_left, n_{d-1}, 1)
    cores.append(C.reshape(r_left, shape[-1], 1))

    return TensorTrain(cores)


def matrix_to_ttm(
    matrix: np.ndarray,
    row_shape: Sequence[int],
    col_shape: Sequence[int],
    eps: float | None = None,
    max_ranks: Sequence[int] | None = None,
) -> TTMatrix:
    """Decompose a dense matrix into TT-matrix format.

    The matrix is reshaped into a tensor with alternating row/col modes,
    decomposed via TT-SVD, and then the TT cores are reshaped into
    4-D TTm cores.

    Parameters
    ----------
    matrix : ndarray
        A 2-D array of shape ``(prod(row_shape), prod(col_shape))``.
    row_shape : sequence of int
        Row mode sizes ``(m_1, ..., m_d)``.
    col_shape : sequence of int
        Column mode sizes ``(n_1, ..., n_d)``.
    eps : float, optional
        Relative Frobenius-norm tolerance.
    max_ranks : sequence of int, optional
        Maximum TT-ranks.

    Returns
    -------
    TTMatrix
    """
    from tensortrain.convert import tt_to_ttm

    matrix = np.asarray(matrix, dtype=np.float64)
    row_shape = tuple(int(m) for m in row_shape)
    col_shape = tuple(int(n) for n in col_shape)
    d = len(row_shape)

    if len(col_shape) != d:
        raise ValueError("row_shape and col_shape must have the same length.")

    M = int(np.prod(row_shape))
    N = int(np.prod(col_shape))
    if matrix.shape != (M, N):
        raise ValueError(
            f"Matrix shape {matrix.shape} does not match "
            f"({M}, {N}) from row/col shapes."
        )

    # Reshape matrix to tensor with separated row/col dims:
    # (m1, m2, ..., md, n1, n2, ..., nd)
    tensor = matrix.reshape(row_shape + col_shape)

    # Permute to interleave: (m1, n1, m2, n2, ..., md, nd)
    perm = []
    for k in range(d):
        perm.append(k)        # m_k
        perm.append(d + k)    # n_k
    tensor = np.transpose(tensor, perm)

    # Merge each (m_k, n_k) pair: shape (m1*n1, m2*n2, ..., md*nd)
    merged_shape = tuple(row_shape[k] * col_shape[k] for k in range(d))
    tensor = tensor.reshape(merged_shape)

    # TT-SVD on the merged tensor
    tt = tensor_to_tt(tensor, eps=eps, max_ranks=max_ranks)

    # Reshape TT cores to TTm cores
    return tt_to_ttm(tt, row_shape, col_shape)


def _truncation_rank_eps(sigma: np.ndarray, delta: float) -> int:
    """Find the smallest rank *r* such that the truncation error is <= delta.

    The truncation error for rank *r* is
    ``sqrt(sigma[r]^2 + sigma[r+1]^2 + ... + sigma[-1]^2)``.
    We want the smallest *r* such that this error is at most *delta*.

    Parameters
    ----------
    sigma : 1-D ndarray
        Singular values in descending order.
    delta : float
        Absolute error threshold.

    Returns
    -------
    int
        The truncation rank (>= 1).
    """
    # Cumulative sum of squares from the right (tail error squared)
    tail_sq = np.cumsum(sigma[::-1] ** 2)[::-1]
    # tail_sq[i] = sum_{j >= i} sigma[j]^2
    # We want the smallest r such that tail_sq[r] <= delta^2
    # i.e., discarding sigma[r:] has error <= delta
    delta_sq = delta ** 2
    for r in range(1, len(sigma)):
        if tail_sq[r] <= delta_sq:
            return r
    return len(sigma)
