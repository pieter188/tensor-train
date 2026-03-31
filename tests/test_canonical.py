"""Tests for tensortrain.canonical."""

import numpy as np
import pytest

from tensortrain import TensorTrain
from tensortrain.canonical import orthogonalize


class TestOrthogonalize:
    """Site-k canonical form."""

    def test_preserves_tensor(self, random_tt):
        """Orthogonalization must not change the represented tensor."""
        X = random_tt.full()
        for site in range(random_tt.ndim):
            tt_orth = orthogonalize(random_tt, site)
            np.testing.assert_allclose(tt_orth.full(), X, atol=1e-12)

    def test_norm_in_site_core(self, random_tt):
        """The Frobenius norm of the tensor equals the norm of the site core."""
        X = random_tt.full()
        tensor_norm = np.linalg.norm(X)
        for site in range(random_tt.ndim):
            tt_orth = orthogonalize(random_tt, site)
            core_norm = np.linalg.norm(tt_orth.core(site))
            np.testing.assert_allclose(core_norm, tensor_norm, rtol=1e-10)

    def test_left_orthogonal_cores(self, random_tt):
        """Cores left of the site should be left-orthogonal."""
        site = random_tt.ndim - 1  # orthogonalize to last site
        tt_orth = orthogonalize(random_tt, site)
        for k in range(site):
            core = tt_orth.core(k)
            r_left, n_k, r_right = core.shape
            # Reshape to (r_left * n_k, r_right) and check Q^T Q = I
            Q = core.reshape(r_left * n_k, r_right)
            np.testing.assert_allclose(
                Q.T @ Q, np.eye(r_right), atol=1e-12
            )

    def test_right_orthogonal_cores(self, random_tt):
        """Cores right of the site should be right-orthogonal."""
        site = 0  # orthogonalize to first site
        tt_orth = orthogonalize(random_tt, site)
        for k in range(1, random_tt.ndim):
            core = tt_orth.core(k)
            r_left, n_k, r_right = core.shape
            # Reshape to (r_left, n_k * r_right) and check Q Q^T = I
            Q = core.reshape(r_left, n_k * r_right)
            np.testing.assert_allclose(
                Q @ Q.T, np.eye(r_left), atol=1e-12
            )

    def test_middle_site(self, rng):
        """Works for a middle site."""
        tt = TensorTrain.random((4, 3, 5, 2), (2, 3, 2), rng=rng)
        X = tt.full()
        tt_orth = orthogonalize(tt, 2)
        np.testing.assert_allclose(tt_orth.full(), X, atol=1e-12)
        # Core 2 should carry the norm
        np.testing.assert_allclose(
            np.linalg.norm(tt_orth.core(2)),
            np.linalg.norm(X),
            rtol=1e-10,
        )

    def test_method_delegation(self, random_tt):
        """TensorTrain.orthogonalize() delegates correctly."""
        tt_orth = random_tt.orthogonalize(0)
        np.testing.assert_allclose(tt_orth.full(), random_tt.full(), atol=1e-12)

    def test_invalid_site(self, random_tt):
        with pytest.raises(ValueError, match="site"):
            orthogonalize(random_tt, random_tt.ndim)
        with pytest.raises(ValueError, match="site"):
            orthogonalize(random_tt, -1)

    def test_single_core(self):
        """Orthogonalizing a single-core TT to site 0 is a no-op."""
        tt = TensorTrain([np.array([[[1.0], [2.0], [3.0]]])])
        tt_orth = orthogonalize(tt, 0)
        np.testing.assert_allclose(tt_orth.full(), tt.full())
