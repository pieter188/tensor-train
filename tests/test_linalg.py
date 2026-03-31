"""Tests for tensortrain.linalg."""

import numpy as np
import pytest

from tensortrain.linalg import hosvd
from tensortrain.utils import nmode_product


class TestHOSVD:
    """Higher-Order SVD (Tucker decomposition)."""

    def test_reconstruction(self, rng):
        """Core times all factors should approximately recover the tensor."""
        X = rng.standard_normal((4, 5, 3))
        core, factors = hosvd(X, ranks=(4, 5, 3))  # full rank = exact
        # Reconstruct
        recon = core
        for k, U in enumerate(factors):
            recon = nmode_product(recon, U, k)
        np.testing.assert_allclose(recon, X, atol=1e-12)

    def test_truncated(self, rng):
        """Truncated HOSVD gives a low-rank approximation."""
        X = rng.standard_normal((6, 8, 5))
        core, factors = hosvd(X, ranks=(3, 4, 2))
        assert core.shape == (3, 4, 2)
        assert factors[0].shape == (6, 3)
        assert factors[1].shape == (8, 4)
        assert factors[2].shape == (5, 2)

    def test_orthogonal_factors(self, rng):
        """Factor matrices should have orthonormal columns."""
        X = rng.standard_normal((4, 5, 3))
        _, factors = hosvd(X, ranks=(3, 4, 2))
        for U in factors:
            np.testing.assert_allclose(U.T @ U, np.eye(U.shape[1]), atol=1e-12)

    def test_rank_clamping(self):
        """Requested ranks larger than mode size should be clamped."""
        X = np.random.randn(3, 4, 2)
        core, factors = hosvd(X, ranks=(10, 10, 10))
        assert core.shape == (3, 4, 2)

    def test_wrong_ranks_length(self):
        with pytest.raises(ValueError, match="ranks"):
            hosvd(np.zeros((3, 4, 5)), ranks=(2, 3))
