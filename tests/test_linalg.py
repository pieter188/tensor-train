"""Tests for tensortrain.linalg."""

import numpy as np
import pytest

from tensortrain import TensorTrain
from tensortrain.arithmetic import dot
from tensortrain.linalg import gram_schmidt_tt, hosvd
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


class TestGramSchmidtTT:
    """Modified Gram-Schmidt on TT vectors."""

    def test_orthonormality(self, rng):
        """Resulting vectors should be approximately orthonormal."""
        n = 5
        vecs = [TensorTrain.random((4, 3), (2,), rng=rng) for _ in range(n)]
        q = gram_schmidt_tt(vecs, eps=1e-10, max_rank=20)
        assert len(q) == n
        for i in range(n):
            for j in range(n):
                expected = 1.0 if i == j else 0.0
                np.testing.assert_allclose(
                    dot(q[i], q[j]), expected, atol=1e-6
                )

    def test_span_preserved(self, rng):
        """GS output should span approximately the same space as input."""
        # Build vectors from known dense vectors, check reconstruction
        M = rng.standard_normal((12, 4))
        Q_dense, _ = np.linalg.qr(M)  # 12x4 orthonormal

        vecs = []
        for j in range(4):
            tt = TensorTrain.from_tensor(M[:, j].reshape(3, 4))
            vecs.append(tt)

        q = gram_schmidt_tt(vecs, eps=1e-10, max_rank=20)
        # Each original vector should be expressible as linear combo of q
        for j in range(4):
            coeffs = [dot(q[i], vecs[j]) for i in range(4)]
            recon = q[0] * coeffs[0]
            for i in range(1, 4):
                recon = recon + q[i] * coeffs[i]
            recon_rounded = recon.round(eps=1e-10)
            np.testing.assert_allclose(
                recon_rounded.full(), vecs[j].full(), atol=1e-4
            )

    def test_single_vector(self):
        """GS on one vector just normalizes it."""
        v = TensorTrain.random((3, 4), (2,), rng=42)
        q = gram_schmidt_tt([v], eps=1e-10)
        np.testing.assert_allclose(dot(q[0], q[0]), 1.0, atol=1e-10)
