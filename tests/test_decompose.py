"""Tests for tensortrain.decompose."""

import numpy as np
import pytest

from tensortrain import TensorTrain
from tensortrain.decompose import tensor_to_tt


class TestTensorToTT:
    """TT-SVD decomposition."""

    def test_exact_round_trip_3d(self, tensor_3d):
        """No-truncation decomposition reconstructs exactly."""
        tt = tensor_to_tt(tensor_3d)
        np.testing.assert_allclose(tt.full(), tensor_3d, atol=1e-12)

    def test_exact_round_trip_4d(self, tensor_4d):
        tt = tensor_to_tt(tensor_4d)
        np.testing.assert_allclose(tt.full(), tensor_4d, atol=1e-12)

    def test_exact_round_trip_2d(self, rng):
        """Works on matrices (2-D tensors)."""
        M = rng.standard_normal((5, 7))
        tt = tensor_to_tt(M)
        np.testing.assert_allclose(tt.full(), M, atol=1e-12)

    def test_eps_truncation_error_bound(self, tensor_3d):
        """Relative error is within the requested epsilon."""
        eps = 1e-4
        tt = tensor_to_tt(tensor_3d, eps=eps)
        recon = tt.full()
        rel_error = np.linalg.norm(recon - tensor_3d) / np.linalg.norm(tensor_3d)
        assert rel_error <= eps

    def test_eps_truncation_reduces_ranks(self, rng):
        """Epsilon truncation should produce lower ranks than exact."""
        # Build a tensor with decaying singular values (low effective rank)
        X = np.zeros((6, 8, 5))
        for _ in range(3):
            a = rng.standard_normal(6)
            b = rng.standard_normal(8)
            c = rng.standard_normal(5)
            X += a[:, None, None] * b[None, :, None] * c[None, None, :]

        tt_exact = tensor_to_tt(X)
        tt_approx = tensor_to_tt(X, eps=1e-6)

        # Approximate should have equal or smaller ranks
        for r_exact, r_approx in zip(tt_exact.ranks, tt_approx.ranks):
            assert r_approx <= r_exact

    def test_max_ranks(self, tensor_3d):
        """max_ranks limits the TT-ranks."""
        tt = tensor_to_tt(tensor_3d, max_ranks=[2, 2])
        assert tt.ranks[1] <= 2
        assert tt.ranks[2] <= 2

    def test_max_ranks_with_eps(self, rng):
        """Both eps and max_ranks applied simultaneously."""
        X = rng.standard_normal((4, 5, 3))
        tt = tensor_to_tt(X, eps=1e-2, max_ranks=[2, 2])
        assert tt.ranks[1] <= 2
        assert tt.ranks[2] <= 2

    def test_rank1_tensor(self):
        """A rank-1 tensor should decompose to rank-1 TT."""
        a = np.array([1.0, 2.0, 3.0])
        b = np.array([4.0, 5.0])
        c = np.array([6.0, 7.0, 8.0, 9.0])
        X = a[:, None, None] * b[None, :, None] * c[None, None, :]
        tt = tensor_to_tt(X, eps=1e-12)
        assert tt.ranks == (1, 1, 1, 1)
        np.testing.assert_allclose(tt.full(), X, atol=1e-10)

    def test_1d_raises(self):
        with pytest.raises(ValueError, match="at least 2"):
            tensor_to_tt(np.array([1.0, 2.0, 3.0]))

    def test_wrong_max_ranks_length(self):
        with pytest.raises(ValueError, match="max_ranks"):
            tensor_to_tt(np.zeros((3, 4, 5)), max_ranks=[2])

    def test_5d_tensor(self, rng):
        """Works for higher-order tensors."""
        X = rng.standard_normal((2, 3, 4, 3, 2))
        tt = tensor_to_tt(X, eps=1e-10)
        np.testing.assert_allclose(tt.full(), X, atol=1e-8)

    def test_shape_preserved(self, tensor_4d):
        tt = tensor_to_tt(tensor_4d, eps=1e-8)
        assert tt.shape == tensor_4d.shape


class TestFromTensor:
    """TensorTrain.from_tensor class method."""

    def test_basic(self, tensor_3d):
        tt = TensorTrain.from_tensor(tensor_3d, eps=1e-10)
        np.testing.assert_allclose(tt.full(), tensor_3d, atol=1e-8)

    def test_exact(self, tensor_3d):
        tt = TensorTrain.from_tensor(tensor_3d)
        np.testing.assert_allclose(tt.full(), tensor_3d, atol=1e-12)


class TestRoundTrips:
    """Reconstruction via .full()."""

    def test_ones_tensor(self):
        tt = TensorTrain.ones((3, 4, 5))
        np.testing.assert_allclose(tt.full(), np.ones((3, 4, 5)))

    def test_zeros_tensor(self):
        tt = TensorTrain.zeros((3, 4, 5))
        np.testing.assert_allclose(tt.full(), np.zeros((3, 4, 5)))

    def test_random_tt_reconstruction(self):
        """A random TT -> full -> TT should round-trip."""
        tt = TensorTrain.random((4, 3, 5), (2, 3), rng=42)
        X = tt.full()
        tt2 = TensorTrain.from_tensor(X, eps=1e-12)
        np.testing.assert_allclose(tt2.full(), X, atol=1e-10)
