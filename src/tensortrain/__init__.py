"""tensortrain — Tensor Train decomposition library.

A pure-Python library (built on NumPy/SciPy) for constructing, manipulating,
and computing with tensors in the Tensor Train (TT) format.

Core classes
------------
TensorTrain
    Tensor Train representation of a *d*-dimensional tensor.
TTMatrix
    Tensor Train Matrix (TTm) representation of a structured matrix.

Examples
--------
>>> import numpy as np
>>> import tensortrain as tt

Decompose a tensor:

>>> X = np.random.randn(4, 5, 3)
>>> t = tt.TensorTrain.from_tensor(X, eps=1e-8)
>>> np.allclose(t.full(), X)
True

Arithmetic in TT format:

>>> a = tt.TensorTrain.random((4, 5, 3), (2, 3), rng=0)
>>> b = tt.TensorTrain.random((4, 5, 3), (2, 2), rng=1)
>>> c = a + b
>>> c_rounded = c.round(eps=1e-6)
"""

__version__ = "0.1.0"

# Core classes
from tensortrain._core import TensorTrain
from tensortrain._matrix import TTMatrix

# Decomposition
from tensortrain.decompose import matrix_to_ttm, qtt, tensor_to_tt

# Arithmetic
from tensortrain.arithmetic import add, concat_ttm, dot, matmat, matvec, sub

# Rounding
from tensortrain.rounding import tt_round, tt_round_to_ranks

# Conversion
from tensortrain.convert import (
    combine_tt_vectors,
    extract_column,
    transpose_ttm,
    tt_to_full,
    tt_to_ttm,
    ttm_to_full,
    ttm_to_tt,
)

# Canonical forms
from tensortrain.canonical import orthogonalize

# Linear algebra
from tensortrain.linalg import gram_schmidt_tt, hosvd

__all__ = [
    # Classes
    "TensorTrain",
    "TTMatrix",
    # Decomposition
    "tensor_to_tt",
    "matrix_to_ttm",
    "qtt",
    # Arithmetic
    "add",
    "sub",
    "dot",
    "matvec",
    "matmat",
    # Rounding
    "tt_round",
    "tt_round_to_ranks",
    # Conversion
    "tt_to_full",
    "ttm_to_full",
    "tt_to_ttm",
    "ttm_to_tt",
    "transpose_ttm",
    # Canonical
    "orthogonalize",
    # Linear algebra
    "hosvd",
]
