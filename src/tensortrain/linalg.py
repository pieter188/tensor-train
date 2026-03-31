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
