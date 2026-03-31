"""Arithmetic operations on Tensor Trains.

All operations work directly on the TT cores without reconstructing the
full tensor.  Addition and subtraction increase the TT-ranks; use
:func:`tensortrain.rounding.tt_round` afterwards to compress.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from tensortrain._validation import check_compatible

if TYPE_CHECKING:
    from tensortrain._core import TensorTrain
    from tensortrain._matrix import TTMatrix


def add(a: TensorTrain, b: TensorTrain) -> TensorTrain:
    """Add two tensor trains: ``C = A + B``.

    The resulting TT has ranks that are the *sum* of the input ranks
    (except at the boundaries which remain 1).  Rounding is recommended
    afterwards.

    Parameters
    ----------
    a, b : TensorTrain
        Must have the same mode sizes.

    Returns
    -------
    TensorTrain

    Notes
    -----
    For core *k* (0-indexed, d cores total):

    - k = 0 (first):  ``C_0[:, i, :] = [A_0[:, i, :] , B_0[:, i, :]]``
      (horizontal concatenation along the right rank).
    - k = d-1 (last):  ``C_d[:, i, :] = [A_d[:, i, :] ; B_d[:, i, :]]``
      (vertical concatenation along the left rank).
    - otherwise:  block-diagonal ``[[A_k, 0], [0, B_k]]`` per index slice.

    References
    ----------
    .. [1] I. V. Oseledets, "Tensor-Train Decomposition", *SIAM J. Sci.
       Comput.*, 33(5), 2295-2317, 2011. Eq. (1.3).
    """
    from tensortrain._core import TensorTrain

    check_compatible(a, b)
    d = a.ndim
    cores = []

    for k in range(d):
        ca = a.core(k)  # (ra_left, n_k, ra_right)
        cb = b.core(k)  # (rb_left, n_k, rb_right)
        n_k = ca.shape[1]

        if k == 0:
            # First core: concatenate along right rank (axis 2)
            # Shape: (1, n_k, ra_right + rb_right)
            core = np.concatenate([ca, cb], axis=2)
        elif k == d - 1:
            # Last core: concatenate along left rank (axis 0)
            # Shape: (ra_left + rb_left, n_k, 1)
            core = np.concatenate([ca, cb], axis=0)
        else:
            # Middle cores: block-diagonal per index slice
            # Shape: (ra_left + rb_left, n_k, ra_right + rb_right)
            ra_l, _, ra_r = ca.shape
            rb_l, _, rb_r = cb.shape
            core = np.zeros((ra_l + rb_l, n_k, ra_r + rb_r))
            core[:ra_l, :, :ra_r] = ca
            core[ra_l:, :, ra_r:] = cb

        cores.append(core)

    return TensorTrain(cores)


def sub(a: TensorTrain, b: TensorTrain) -> TensorTrain:
    """Subtract two tensor trains: ``C = A - B``.

    Equivalent to ``add(a, -b)`` but avoids an extra copy.

    Parameters
    ----------
    a, b : TensorTrain

    Returns
    -------
    TensorTrain
    """
    return add(a, -b)


def dot(a: TensorTrain, b: TensorTrain) -> float:
    """Inner product of two tensor trains: ``<A, B>``.

    Computed by contracting corresponding cores left-to-right without
    forming the full tensors.

    Parameters
    ----------
    a, b : TensorTrain
        Must have the same mode sizes.

    Returns
    -------
    float
        The scalar inner product.

    Notes
    -----
    At each step *k*, we maintain a matrix ``M`` of shape
    ``(ra_k, rb_k)`` that accumulates the partial contraction over
    modes ``0, 1, ..., k-1``.  The contraction over mode *k* is:

    .. math::

        M'[\\alpha', \\beta'] = \\sum_{\\alpha, i_k, \\beta}
        M[\\alpha, \\beta] \\, A^{(k)}[\\alpha, i_k, \\alpha']
        \\, B^{(k)}[\\beta, i_k, \\beta']

    References
    ----------
    .. [1] I. V. Oseledets, "Tensor-Train Decomposition", *SIAM J. Sci.
       Comput.*, 33(5), 2295-2317, 2011.
    """
    check_compatible(a, b)

    # Initialize: M has shape (1, 1) = [[1]]
    M = np.ones((1, 1))

    for k in range(a.ndim):
        ca = a.core(k)  # (ra_l, n_k, ra_r)
        cb = b.core(k)  # (rb_l, n_k, rb_r)

        # Contract M with ca and cb over the left ranks and mode index
        # M[alpha, beta] * ca[alpha, i, alpha'] * cb[beta, i, beta']
        # = sum over alpha, beta, i
        # Result shape: (ra_r, rb_r)
        M = np.einsum("ab,aic,bid->cd", M, ca, cb)

    return float(M.item())


# ======================================================================
# TTMatrix arithmetic
# ======================================================================


def matvec(A: TTMatrix, x: TensorTrain) -> TensorTrain:
    """Matrix-vector product ``A @ x`` in TT format.

    Parameters
    ----------
    A : TTMatrix
        Matrix with ``A.col_shape == x.shape``.
    x : TensorTrain
        Vector.

    Returns
    -------
    TensorTrain
        Result with ``shape == A.row_shape`` and ranks that are the
        products of the input ranks.

    Notes
    -----
    For each core *k*, the contraction is:

    .. math::

        C^{(k)}[\\alpha \\beta, m, \\alpha' \\beta']
        = \\sum_n A^{(k)}[\\alpha, m, n, \\alpha'] \\, x^{(k)}[\\beta, n, \\beta']

    References
    ----------
    .. [1] I. V. Oseledets, "Tensor-Train Decomposition", *SIAM J. Sci.
       Comput.*, 33(5), 2295-2317, 2011.
    """
    from tensortrain._core import TensorTrain

    if A.col_shape != x.shape:
        raise ValueError(
            f"Column shape {A.col_shape} does not match vector shape {x.shape}."
        )

    cores = []
    for k in range(A.ndim):
        ca = A.core(k)   # (ra_l, m_k, n_k, ra_r)
        cx = x.core(k)   # (rx_l, n_k, rx_r)
        ra_l, m_k, n_k, ra_r = ca.shape
        rx_l, _, rx_r = cx.shape

        # Contract over n_k: result shape (ra_l, rx_l, m_k, ra_r, rx_r)
        result = np.einsum("amnb,cnd->acmbd", ca, cx)
        # Reshape to (ra_l * rx_l, m_k, ra_r * rx_r)
        cores.append(result.reshape(ra_l * rx_l, m_k, ra_r * rx_r))

    return TensorTrain(cores)


def matmat(A: TTMatrix, B: TTMatrix) -> TTMatrix:
    """Matrix-matrix product ``A @ B`` in TTm format.

    Parameters
    ----------
    A : TTMatrix
        Left matrix with ``A.col_shape == B.row_shape``.
    B : TTMatrix
        Right matrix.

    Returns
    -------
    TTMatrix
        Result with ``row_shape == A.row_shape``,
        ``col_shape == B.col_shape``, and ranks that are the products
        of the input ranks.
    """
    from tensortrain._matrix import TTMatrix

    if A.col_shape != B.row_shape:
        raise ValueError(
            f"A.col_shape {A.col_shape} does not match B.row_shape {B.row_shape}."
        )

    cores = []
    for k in range(A.ndim):
        ca = A.core(k)   # (ra_l, m_k, n_k, ra_r)
        cb = B.core(k)   # (rb_l, n_k, p_k, rb_r)
        ra_l, m_k, n_k, ra_r = ca.shape
        rb_l, _, p_k, rb_r = cb.shape

        # Contract over n_k: result shape (ra_l, rb_l, m_k, p_k, ra_r, rb_r)
        result = np.einsum("amnb,cnpd->acmpbd", ca, cb)
        # Reshape to (ra_l * rb_l, m_k, p_k, ra_r * rb_r)
        cores.append(result.reshape(ra_l * rb_l, m_k, p_k, ra_r * rb_r))

    return TTMatrix(cores)


def add_ttm(a: TTMatrix, b: TTMatrix) -> TTMatrix:
    """Add two TTMatrices: ``C = A + B``.

    Same block-diagonal structure as TT addition, but with 4-D cores.

    Parameters
    ----------
    a, b : TTMatrix
        Must have the same row_shape and col_shape.

    Returns
    -------
    TTMatrix
    """
    from tensortrain._matrix import TTMatrix

    if a.row_shape != b.row_shape or a.col_shape != b.col_shape:
        raise ValueError(
            f"Incompatible shapes: ({a.row_shape}, {a.col_shape}) "
            f"vs ({b.row_shape}, {b.col_shape})."
        )

    d = a.ndim
    cores = []

    for k in range(d):
        ca = a.core(k)  # (ra_l, m, n, ra_r)
        cb = b.core(k)  # (rb_l, m, n, rb_r)
        m_k = ca.shape[1]
        n_k = ca.shape[2]

        if k == 0:
            core = np.concatenate([ca, cb], axis=3)
        elif k == d - 1:
            core = np.concatenate([ca, cb], axis=0)
        else:
            ra_l, _, _, ra_r = ca.shape
            rb_l, _, _, rb_r = cb.shape
            core = np.zeros((ra_l + rb_l, m_k, n_k, ra_r + rb_r))
            core[:ra_l, :, :, :ra_r] = ca
            core[ra_l:, :, :, ra_r:] = cb

        cores.append(core)

    return TTMatrix(cores)


def sub_ttm(a: TTMatrix, b: TTMatrix) -> TTMatrix:
    """Subtract two TTMatrices: ``C = A - B``."""
    return add_ttm(a, -b)
