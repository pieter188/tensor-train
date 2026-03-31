"""Conversion functions between TT format and dense tensors/matrices.

This module handles reconstruction (TT -> full tensor) and reshaping
operations between different tensor/matrix representations.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from tensortrain._core import TensorTrain
    from tensortrain._matrix import TTMatrix


def tt_to_full(tt: TensorTrain) -> np.ndarray:
    """Reconstruct the full dense tensor from a TensorTrain.

    Contracts all cores from left to right via sequential matrix
    multiplications.

    Parameters
    ----------
    tt : TensorTrain
        The tensor train to reconstruct.

    Returns
    -------
    ndarray
        Dense tensor of shape ``tt.shape``.

    Notes
    -----
    This materialises the full tensor, which may be very large for
    high-dimensional problems.  Use only for moderate sizes or testing.

    Examples
    --------
    >>> import numpy as np
    >>> from tensortrain import TensorTrain
    >>> tt = TensorTrain.ones((3, 4, 5))
    >>> np.allclose(tt.full(), np.ones((3, 4, 5)))
    True
    """
    d = tt.ndim
    # Start with first core: shape (1, n_0, r_1) -> (n_0, r_1)
    result = tt.core(0).reshape(tt.core(0).shape[1], tt.core(0).shape[2])

    for k in range(1, d):
        core_k = tt.core(k)  # (r_k, n_k, r_{k+1})
        r_k, n_k, r_next = core_k.shape

        # result has shape (n_0 * n_1 * ... * n_{k-1}, r_k)
        # core_k reshaped to (r_k, n_k * r_{k+1})
        result = result @ core_k.reshape(r_k, n_k * r_next)
        # Now shape is (n_0 * ... * n_{k-1}, n_k * r_{k+1})
        # Reshape to (n_0 * ... * n_k, r_{k+1})
        result = result.reshape(-1, r_next)

    # Final shape: (n_0 * n_1 * ... * n_{d-1}, 1) -> tt.shape
    return result.reshape(tt.shape)


def matrix_to_tensor(matrix, shape):
    """Reshape a matrix (or vector) into a tensor of the given shape.

    Parameters
    ----------
    matrix : ndarray
        A 1-D or 2-D array whose total number of elements equals
        ``prod(shape)``.
    shape : tuple of int
        Target tensor shape.

    Returns
    -------
    ndarray
        Tensor of the requested shape.
    """
    matrix = np.asarray(matrix, dtype=np.float64)
    return matrix.reshape(shape)


def tensor_to_matrix(tensor, row_modes, col_modes):
    """Reshape a tensor into a matrix by grouping modes.

    Parameters
    ----------
    tensor : ndarray
        A *d*-dimensional array.
    row_modes : tuple of int
        Mode indices that form the row dimension.
    col_modes : tuple of int
        Mode indices that form the column dimension.

    Returns
    -------
    ndarray
        A 2-D array of shape ``(prod of row mode sizes, prod of col mode sizes)``.
    """
    tensor = np.asarray(tensor)
    perm = list(row_modes) + list(col_modes)
    tensor = np.transpose(tensor, perm)
    row_size = int(np.prod([tensor.shape[i] for i in range(len(row_modes))]))
    col_size = int(np.prod([tensor.shape[i] for i in range(len(row_modes), len(perm))]))
    return tensor.reshape(row_size, col_size)


def ttm_to_full(ttm: TTMatrix) -> np.ndarray:
    """Reconstruct the full dense matrix from a TTMatrix.

    Parameters
    ----------
    ttm : TTMatrix
        The tensor train matrix to reconstruct.

    Returns
    -------
    ndarray
        Dense matrix of shape ``ttm.shape``.

    Notes
    -----
    The algorithm converts TTm to TT (merging row/col modes), reconstructs
    the full tensor, then permutes and reshapes to recover the matrix.
    """
    d = ttm.ndim
    row_shape = ttm.row_shape
    col_shape = ttm.col_shape

    # Convert to TT by merging row/col modes
    tt = ttm_to_tt(ttm)

    # Reconstruct the full tensor with merged modes: (m1*n1, m2*n2, ...)
    tensor = tt_to_full(tt)

    # Reshape to separate row and column dims: (m1, n1, m2, n2, ...)
    interleaved_shape = []
    for m, n in zip(row_shape, col_shape):
        interleaved_shape.extend([m, n])
    tensor = tensor.reshape(interleaved_shape)

    # Permute: group all row dims first, then all col dims
    # (m1, n1, m2, n2, ...) -> (m1, m2, ..., n1, n2, ...)
    row_axes = list(range(0, 2 * d, 2))
    col_axes = list(range(1, 2 * d, 2))
    tensor = np.transpose(tensor, row_axes + col_axes)

    # Reshape to (M, N)
    M = int(np.prod(row_shape))
    N = int(np.prod(col_shape))
    return tensor.reshape(M, N)


def tt_to_ttm(
    tt: TensorTrain,
    row_shape: tuple[int, ...],
    col_shape: tuple[int, ...],
) -> TTMatrix:
    """Convert a TensorTrain to a TTMatrix by splitting each mode.

    Each TT core of shape ``(r_{k-1}, m_k * n_k, r_k)`` is reshaped to
    a TTm core of shape ``(r_{k-1}, m_k, n_k, r_k)``.

    Parameters
    ----------
    tt : TensorTrain
        Must have ``tt.shape[k] == row_shape[k] * col_shape[k]`` for all *k*.
    row_shape : tuple of int
        Row mode sizes ``(m_1, ..., m_d)``.
    col_shape : tuple of int
        Column mode sizes ``(n_1, ..., n_d)``.

    Returns
    -------
    TTMatrix
    """
    from tensortrain._matrix import TTMatrix

    d = tt.ndim
    if len(row_shape) != d or len(col_shape) != d:
        raise ValueError("row_shape and col_shape must have length d.")

    cores = []
    for k in range(d):
        r_left, mn, r_right = tt.core(k).shape
        m_k, n_k = row_shape[k], col_shape[k]
        if m_k * n_k != mn:
            raise ValueError(
                f"Core {k}: mode size {mn} != {m_k} * {n_k}."
            )
        cores.append(tt.core(k).reshape(r_left, m_k, n_k, r_right))

    return TTMatrix(cores)


def ttm_to_tt(ttm: TTMatrix) -> TensorTrain:
    """Convert a TTMatrix to a TensorTrain by merging row/col modes.

    Each TTm core of shape ``(r_{k-1}, m_k, n_k, r_k)`` is reshaped to
    a TT core of shape ``(r_{k-1}, m_k * n_k, r_k)``.

    Parameters
    ----------
    ttm : TTMatrix

    Returns
    -------
    TensorTrain
    """
    from tensortrain._core import TensorTrain

    cores = []
    for k in range(ttm.ndim):
        r_left, m_k, n_k, r_right = ttm.core(k).shape
        cores.append(ttm.core(k).reshape(r_left, m_k * n_k, r_right))

    return TensorTrain(cores)


def transpose_ttm(ttm: TTMatrix) -> TTMatrix:
    """Transpose a TTMatrix by swapping row and column modes in each core.

    Parameters
    ----------
    ttm : TTMatrix

    Returns
    -------
    TTMatrix
        The transposed matrix.
    """
    from tensortrain._matrix import TTMatrix

    cores = []
    for k in range(ttm.ndim):
        # Swap axes 1 (row) and 2 (col): (r_left, m, n, r_right) -> (r_left, n, m, r_right)
        cores.append(np.transpose(ttm.core(k), (0, 2, 1, 3)).copy())
    return TTMatrix(cores)


def vector_to_tensor(vector, shape):
    """Reshape a vector into a tensor.

    Parameters
    ----------
    vector : ndarray
        A 1-D array.
    shape : tuple of int
        Target tensor shape.

    Returns
    -------
    ndarray
    """
    vector = np.asarray(vector, dtype=np.float64)
    return vector.reshape(shape)
