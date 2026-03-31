"""Utility functions for tensor operations.

Provides fundamental tensor operations that are used throughout the library
and are also useful on their own: mode-n matricization, n-mode products,
and outer products.
"""

import numpy as np


def matricize(tensor, mode):
    """Mode-*n* matricization (unfolding) of a tensor.

    Unfolds ``tensor`` along ``mode`` into a matrix whose rows correspond to
    the given mode and whose columns cycle through the remaining modes in the
    standard (row-major for the complement) order.

    Parameters
    ----------
    tensor : ndarray
        A *d*-dimensional array.
    mode : int
        The mode along which to unfold (0-indexed).

    Returns
    -------
    ndarray
        A 2-D array of shape ``(I_mode, prod(I_k for k != mode))``.

    Examples
    --------
    >>> import numpy as np
    >>> from tensortrain.utils import matricize
    >>> X = np.arange(24).reshape(2, 3, 4)
    >>> matricize(X, 0).shape
    (2, 12)
    >>> matricize(X, 1).shape
    (3, 8)
    """
    tensor = np.asarray(tensor)
    d = tensor.ndim
    if not 0 <= mode < d:
        raise ValueError(
            f"mode must be in [0, {d}), got {mode}."
        )

    # Move the target mode to front, then flatten the rest
    # Order of remaining modes: 0, 1, ..., mode-1, mode+1, ..., d-1
    axes = [mode] + list(range(mode)) + list(range(mode + 1, d))
    return np.reshape(np.transpose(tensor, axes), (tensor.shape[mode], -1))


def nmode_product(tensor, matrix, mode):
    """N-mode product of a tensor with a matrix.

    Computes the *n*-mode product :math:`\\mathcal{X} \\times_n U`, which
    multiplies each mode-*n* fiber of the tensor by the matrix *U*.

    Parameters
    ----------
    tensor : ndarray
        A *d*-dimensional array with mode sizes ``(I_1, ..., I_d)``.
    matrix : ndarray
        A 2-D array of shape ``(J, I_mode)`` where ``I_mode`` is the size
        of ``tensor`` along ``mode``.
    mode : int
        The mode along which to multiply (0-indexed).

    Returns
    -------
    ndarray
        Array of shape ``(I_1, ..., I_{n-1}, J, I_{n+1}, ..., I_d)``.

    Examples
    --------
    >>> import numpy as np
    >>> from tensortrain.utils import nmode_product
    >>> X = np.random.randn(3, 4, 5)
    >>> U = np.random.randn(2, 4)
    >>> Y = nmode_product(X, U, mode=1)
    >>> Y.shape
    (3, 2, 5)
    """
    tensor = np.asarray(tensor)
    matrix = np.asarray(matrix)
    d = tensor.ndim

    if not 0 <= mode < d:
        raise ValueError(f"mode must be in [0, {d}), got {mode}.")
    if matrix.ndim != 2:
        raise ValueError(f"matrix must be 2-D, got ndim={matrix.ndim}.")
    if matrix.shape[1] != tensor.shape[mode]:
        raise ValueError(
            f"matrix columns ({matrix.shape[1]}) must match tensor mode "
            f"size ({tensor.shape[mode]})."
        )

    # Use np.tensordot: contract axis 1 of matrix with the given mode of tensor
    result = np.tensordot(matrix, tensor, axes=([1], [mode]))
    # tensordot places the matrix's axis-0 first; we need to move it to `mode`
    # Current order: (J, I_0, ..., I_{mode-1}, I_{mode+1}, ..., I_{d-1})
    # Target order:  (I_0, ..., I_{mode-1}, J, I_{mode+1}, ..., I_{d-1})
    if mode > 0:
        # Move axis 0 to position `mode`
        result = np.moveaxis(result, 0, mode)

    return result


def outer_product(*vectors):
    """Outer product of an arbitrary number of vectors.

    Computes the tensor :math:`\\mathcal{X}` such that
    :math:`\\mathcal{X}(i_1, i_2, \\ldots, i_d) = v_1(i_1) \\cdot v_2(i_2)
    \\cdots v_d(i_d)`.

    Parameters
    ----------
    *vectors : ndarray
        One-dimensional arrays.

    Returns
    -------
    ndarray
        A *d*-dimensional array of shape ``(len(v1), len(v2), ...)``.

    Examples
    --------
    >>> import numpy as np
    >>> from tensortrain.utils import outer_product
    >>> a = np.array([1.0, 2.0])
    >>> b = np.array([3.0, 4.0, 5.0])
    >>> outer_product(a, b).shape
    (2, 3)
    >>> np.allclose(outer_product(a, b), np.outer(a, b))
    True
    """
    if len(vectors) < 2:
        raise ValueError("Need at least 2 vectors for an outer product.")

    for k, v in enumerate(vectors):
        v = np.asarray(v)
        if v.ndim != 1:
            raise ValueError(f"Vector {k} must be 1-D, got ndim={v.ndim}.")

    # Build up via successive outer products using broadcasting
    result = np.asarray(vectors[0])
    for v in vectors[1:]:
        v = np.asarray(v)
        # Expand result with a new trailing axis, v with leading axes
        result = result[..., np.newaxis] * v[tuple(np.newaxis for _ in range(result.ndim)) + (slice(None),)]

    return result


def multi_index(shape, flat_index):
    """Convert a flat (linear) index to a multi-index.

    Parameters
    ----------
    shape : tuple of int
        The tensor shape.
    flat_index : int
        Linear index (row-major / C-order).

    Returns
    -------
    tuple of int
        The corresponding multi-index.

    Examples
    --------
    >>> from tensortrain.utils import multi_index
    >>> multi_index((3, 4, 5), 0)
    (0, 0, 0)
    >>> multi_index((3, 4, 5), 59)
    (2, 3, 4)
    """
    return tuple(int(i) for i in np.unravel_index(flat_index, shape))
