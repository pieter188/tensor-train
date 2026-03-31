"""Shared test fixtures for the tensortrain test suite."""

import numpy as np
import pytest

from tensortrain import TensorTrain


@pytest.fixture
def rng():
    """Reproducible random generator."""
    return np.random.default_rng(42)


@pytest.fixture
def tensor_3d(rng):
    """A random 4x3x2 tensor."""
    return rng.standard_normal((4, 3, 2))


@pytest.fixture
def tensor_4d(rng):
    """A random 3x4x2x5 tensor."""
    return rng.standard_normal((3, 4, 2, 5))


@pytest.fixture
def rank1_tt():
    """A rank-1 tensor train (3x4x2) built from known vectors."""
    a = np.array([1.0, 2.0, 3.0]).reshape(1, 3, 1)
    b = np.array([4.0, 5.0, 6.0, 7.0]).reshape(1, 4, 1)
    c = np.array([8.0, 9.0]).reshape(1, 2, 1)
    return TensorTrain([a, b, c])


@pytest.fixture
def random_tt(rng):
    """A random TT with shape (4, 3, 5) and ranks (1, 2, 3, 1)."""
    return TensorTrain.random((4, 3, 5), (2, 3), rng=rng)
