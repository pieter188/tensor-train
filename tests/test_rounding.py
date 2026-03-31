"""Tests for tensortrain.rounding."""

import numpy as np
import pytest

from tensortrain import TensorTrain
from tensortrain.rounding import tt_round, tt_round_to_ranks


class TestRound:
    """Epsilon-based rounding."""

    def test_error_within_tolerance(self, rng):
        """Rounding error must be within the requested epsilon."""
        X = rng.standard_normal((4, 5, 3))
        tt = TensorTrain.from_tensor(X)
        eps = 1e-4
        tt_rounded = tt_round(tt, eps=eps)
        recon = tt_rounded.full()
        rel_error = np.linalg.norm(recon - X) / np.linalg.norm(X)
        assert rel_error <= eps

    def test_reduces_ranks(self, rng):
        """Rounding a bloated TT should reduce ranks."""
        # Create a low-rank tensor, decompose exactly, add noise, round
        X = np.zeros((6, 8, 5))
        for _ in range(2):
            a = rng.standard_normal(6)
            b = rng.standard_normal(8)
            c = rng.standard_normal(5)
            X += a[:, None, None] * b[None, :, None] * c[None, None, :]

        tt_exact = TensorTrain.from_tensor(X)
        # Artificially inflate ranks via addition
        tt_bloated = tt_exact + tt_exact  # doubles ranks
        tt_rounded = tt_round(tt_bloated, eps=1e-10)

        # Rounded ranks should be less than bloated ranks
        for r_bloated, r_rounded in zip(tt_bloated.ranks[1:-1], tt_rounded.ranks[1:-1]):
            assert r_rounded <= r_bloated

    def test_preserves_tensor_tight_eps(self, rng):
        """With very small eps, rounding should barely change the tensor."""
        X = rng.standard_normal((4, 3, 5))
        tt = TensorTrain.from_tensor(X)
        tt_rounded = tt_round(tt, eps=1e-14)
        np.testing.assert_allclose(tt_rounded.full(), X, atol=1e-10)

    def test_requires_eps_or_ranks(self, random_tt):
        with pytest.raises(ValueError, match="At least one"):
            tt_round(random_tt)

    def test_method_delegation(self, rng):
        """TensorTrain.round() delegates to tt_round."""
        X = rng.standard_normal((4, 3, 5))
        tt = TensorTrain.from_tensor(X)
        tt_rounded = tt.round(eps=1e-4)
        rel_error = np.linalg.norm(tt_rounded.full() - X) / np.linalg.norm(X)
        assert rel_error <= 1e-4

    def test_idempotent(self, rng):
        """Rounding an already-rounded TT should be stable."""
        X = rng.standard_normal((4, 5, 3))
        tt = TensorTrain.from_tensor(X, eps=1e-6)
        tt_r1 = tt_round(tt, eps=1e-6)
        tt_r2 = tt_round(tt_r1, eps=1e-6)
        np.testing.assert_allclose(tt_r1.full(), tt_r2.full(), atol=1e-12)
        assert tt_r2.ranks == tt_r1.ranks


class TestRoundToRanks:
    """Rank-based rounding."""

    def test_exact_ranks(self, rng):
        """Rounding to specific ranks should produce those ranks."""
        X = rng.standard_normal((4, 5, 3))
        tt = TensorTrain.from_tensor(X)
        target_ranks = [2, 2]
        tt_rounded = tt_round_to_ranks(tt, target_ranks)
        # Ranks should be at most the targets
        for i, (r, target) in enumerate(zip(tt_rounded.ranks[1:-1], target_ranks)):
            assert r <= target

    def test_preserves_shape(self, rng):
        X = rng.standard_normal((4, 5, 3))
        tt = TensorTrain.from_tensor(X)
        tt_rounded = tt_round_to_ranks(tt, [2, 2])
        assert tt_rounded.shape == tt.shape


class TestRoundMaxRanks:
    """Combined eps + max_ranks."""

    def test_both_constraints(self, rng):
        X = rng.standard_normal((4, 5, 3))
        tt = TensorTrain.from_tensor(X)
        tt_rounded = tt_round(tt, eps=1e-2, max_ranks=[2, 2])
        assert tt_rounded.ranks[1] <= 2
        assert tt_rounded.ranks[2] <= 2
