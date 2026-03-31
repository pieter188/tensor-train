"""Higher-level linear algebra operations on tensors.

Includes the Higher-Order SVD (Tucker decomposition) and the
Gram-Schmidt process in TT format.
"""

from __future__ import annotations

import numpy as np
from numpy.linalg import svd

from tensortrain.utils import matricize, nmode_product


def hosvd(tensor, ranks):
    """Higher-Order SVD (Tucker decomposition).

    Computes the Tucker decomposition of a tensor using the truncated
    SVD of each mode-*n* unfolding.

    Parameters
    ----------
    tensor : ndarray
        A *d*-dimensional array.
    ranks : sequence of int
        Truncation ranks ``(R_1, ..., R_d)`` for each mode.

    Returns
    -------
    core : ndarray
        The core tensor of shape ``(R_1, R_2, ..., R_d)``.
    factors : list of ndarray
        Factor matrices ``[U_1, ..., U_d]`` where ``U_k`` has shape
        ``(I_k, R_k)``.

    Notes
    -----
    The original tensor is approximately reconstructed as:

    .. math::

        \\mathcal{X} \\approx \\mathcal{G} \\times_1 U_1 \\times_2 U_2
        \\cdots \\times_d U_d

    References
    ----------
    .. [1] L. De Lathauwer, B. De Moor, J. Vandewalle, "A Multilinear
       Singular Value Decomposition", *SIAM J. Matrix Anal. Appl.*,
       21(4), 1253-1278, 2000.

    Examples
    --------
    >>> import numpy as np
    >>> from tensortrain.linalg import hosvd
    >>> X = np.random.randn(4, 5, 3)
    >>> core, factors = hosvd(X, ranks=(2, 3, 2))
    >>> core.shape
    (2, 3, 2)
    >>> [f.shape for f in factors]
    [(4, 2), (5, 3), (3, 2)]
    """
    tensor = np.asarray(tensor, dtype=np.float64)
    d = tensor.ndim
    ranks = list(ranks)

    if len(ranks) != d:
        raise ValueError(f"Expected {d} ranks, got {len(ranks)}.")

    factors = []
    for k in range(d):
        mat = matricize(tensor, k)
        U, _, _ = svd(mat, full_matrices=False)
        r = min(ranks[k], U.shape[1])
        ranks[k] = r
        factors.append(U[:, :r])

    # Core = tensor contracted with U_k^T along each mode
    core = tensor
    for k in range(d):
        core = nmode_product(core, factors[k].T, k)

    return core, factors


def gram_schmidt_tt(
    vectors: list,
    eps: float = 1e-5,
    max_rank: int = 5,
) -> list:
    """Modified Gram-Schmidt orthogonalization on TT vectors.

    Takes a list of TT vectors and returns a list of orthonormal TT
    vectors spanning (approximately) the same subspace.

    Parameters
    ----------
    vectors : list of TensorTrain
        Input vectors, all with the same shape.
    eps : float, optional
        Rounding tolerance applied when ranks exceed *max_rank*.
    max_rank : int, optional
        Maximum TT-rank before rounding is triggered.

    Returns
    -------
    list of TensorTrain
        Orthonormal TT vectors.

    Notes
    -----
    This is an approximate algorithm: rounding is applied to keep ranks
    manageable, which introduces small deviations from exact
    orthogonality.  The resulting vectors satisfy
    ``|<q_i, q_j> - delta_ij| < O(eps)`` in practice.

    Uses optimized raw-core operations internally for speed.

    References
    ----------
    .. [1] Master thesis, P. van Klaveren — Modified Gram-Schmidt
       orthogonalization in TT format for the TNKF algorithm.
    """
    from tensortrain._core import TensorTrain
    from tensortrain._fast import (
        batched_dot,
        batched_sub_scaled,
        fast_orthogonalize_and_normalize,
        fast_round,
        max_rank as _max_rank,
    )

    n = len(vectors)

    # Work entirely with raw core lists for speed
    v = [[c.copy() for c in vec._cores] for vec in vectors]
    q_cores = [None] * n

    for i in range(n):
        # Orthogonalize to site 0 and normalize
        q_cores[i] = fast_orthogonalize_and_normalize(v[i])

        # Batched: compute all projections dot(q_i, v_j) at once
        remaining = v[i + 1:]
        if remaining:
            projs = batched_dot(q_cores[i], remaining)
            batched_sub_scaled(remaining, q_cores[i], projs)

            # Truncate any vectors whose ranks grew too large
            for j_off in range(len(remaining)):
                if _max_rank(remaining[j_off]) > max_rank:
                    remaining[j_off] = fast_round(remaining[j_off], eps)

    # Wrap results as TensorTrain objects
    return [TensorTrain(cores) for cores in q_cores]
