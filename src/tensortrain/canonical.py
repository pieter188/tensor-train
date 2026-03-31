"""Canonical forms for Tensor Trains.

A TT is in *site-k canonical form* (also called *k-orthogonal*) when:

- All cores to the left of site *k* are left-orthogonal
  (their mode-``(r_left, n)`` unfolding has orthonormal columns).
- All cores to the right of site *k* are right-orthogonal
  (their mode-``(n, r_right)`` unfolding has orthonormal rows).
- The Frobenius norm of the full tensor equals the Frobenius norm of core *k*.

This form is used by the rounding algorithm and for computing norms
efficiently.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from tensortrain._core import TensorTrain


def orthogonalize(tt: TensorTrain, site: int) -> TensorTrain:
    """Bring a TensorTrain into site-*k* canonical form.

    Parameters
    ----------
    tt : TensorTrain
        Input tensor train.
    site : int
        The site index (0-indexed) at which to concentrate the norm.
        Must be in ``[0, d)``.

    Returns
    -------
    TensorTrain
        A new TensorTrain in site-*k* canonical form.

    Notes
    -----
    Left orthogonalization uses QR decomposition; right orthogonalization
    uses QR decomposition as well. The ``R`` factors are absorbed into the
    adjacent core at each step.

    References
    ----------
    .. [1] I. V. Oseledets, "Tensor-Train Decomposition", *SIAM J. Sci.
       Comput.*, 33(5), 2295-2317, 2011. Section 2.2.
    """
    from tensortrain._core import TensorTrain

    d = tt.ndim
    if not 0 <= site < d:
        raise ValueError(f"site must be in [0, {d}), got {site}.")

    cores = [c.copy() for c in tt._cores]

    # Left-to-right QR sweep: cores 0, 1, ..., site-1
    for k in range(site):
        r_left, n_k, r_right = cores[k].shape
        # Reshape to (r_left * n_k, r_right) and do QR
        Q, R = np.linalg.qr(cores[k].reshape(r_left * n_k, r_right))
        r_new = Q.shape[1]
        cores[k] = Q.reshape(r_left, n_k, r_new)
        # Absorb R into the next core
        cores[k + 1] = np.einsum("ij,jkl->ikl", R, cores[k + 1])

    # Right-to-left QR sweep: cores d-1, d-2, ..., site+1
    for k in range(d - 1, site, -1):
        r_left, n_k, r_right = cores[k].shape
        # Reshape to (r_left, n_k * r_right) and do QR on the transpose
        # This is equivalent to an RQ decomposition
        Q, R = np.linalg.qr(cores[k].reshape(r_left, n_k * r_right).T)
        r_new = Q.shape[1]
        cores[k] = Q.T.reshape(r_new, n_k, r_right)
        # Absorb R^T (shape r_left x r_new) into the previous core
        # Contract last axis of cores[k-1] (r_left) with first axis of R.T
        cores[k - 1] = np.einsum("ijk,kl->ijl", cores[k - 1], R.T)

    result = TensorTrain(cores)
    result._norm_index = site
    return result
