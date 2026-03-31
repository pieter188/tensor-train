"""Tests for tensortrain.utils."""

import numpy as np
import pytest

from tensortrain.utils import matricize, multi_index, nmode_product, outer_product


class TestMatricize:
    """Mode-n unfolding."""

    def test_mode0_shape(self):
        X = np.arange(24).reshape(2, 3, 4)
        M = matricize(X, 0)
        assert M.shape == (2, 12)

    def test_mode1_shape(self):
        X = np.arange(24).reshape(2, 3, 4)
        M = matricize(X, 1)
        assert M.shape == (3, 8)

    def test_mode2_shape(self):
        X = np.arange(24).reshape(2, 3, 4)
        M = matricize(X, 2)
        assert M.shape == (4, 6)

    def test_mode0_values(self):
        # For a 2x3 matrix, mode-0 matricization is the matrix itself
        X = np.array([[1, 2, 3], [4, 5, 6]])
        M = matricize(X, 0)
        np.testing.assert_array_equal(M, X)

    def test_mode1_values(self):
        X = np.array([[1, 2, 3], [4, 5, 6]])
        M = matricize(X, 1)
        # Mode-1: rows are column fibers
        np.testing.assert_array_equal(M, np.array([[1, 4], [2, 5], [3, 6]]))

    def test_4d_tensor(self):
        X = np.random.randn(3, 4, 2, 5)
        for mode in range(4):
            M = matricize(X, mode)
            assert M.shape[0] == X.shape[mode]
            assert M.shape[1] == X.size // X.shape[mode]

    def test_invalid_mode(self):
        X = np.zeros((2, 3, 4))
        with pytest.raises(ValueError, match="mode"):
            matricize(X, 3)
        with pytest.raises(ValueError, match="mode"):
            matricize(X, -1)

    def test_round_trip(self):
        """Unfolding and refolding should recover the original tensor."""
        rng = np.random.default_rng(42)
        X = rng.standard_normal((3, 4, 5))
        for mode in range(3):
            M = matricize(X, mode)
            # Reconstruct: the complement modes are in order 0,...,mode-1, mode+1,...,d-1
            shape = list(X.shape)
            n = shape.pop(mode)
            # Transpose back (mode, complement) -> original
            axes = list(range(mode)) + list(range(mode + 1, X.ndim))
            recon = np.reshape(M, [n] + shape)
            # Inverse of the permute [mode, 0, ..., mode-1, mode+1, ..., d-1]
            inv_axes = list(range(1, mode + 1)) + [0] + list(range(mode + 1, X.ndim))
            recon = np.transpose(recon, inv_axes)
            np.testing.assert_allclose(recon, X)


class TestNmodeProduct:
    """N-mode product."""

    def test_basic_shape(self):
        X = np.random.randn(3, 4, 5)
        U = np.random.randn(2, 4)
        Y = nmode_product(X, U, mode=1)
        assert Y.shape == (3, 2, 5)

    def test_mode0(self):
        X = np.random.randn(3, 4, 5)
        U = np.random.randn(6, 3)
        Y = nmode_product(X, U, mode=0)
        assert Y.shape == (6, 4, 5)

    def test_last_mode(self):
        X = np.random.randn(3, 4, 5)
        U = np.random.randn(7, 5)
        Y = nmode_product(X, U, mode=2)
        assert Y.shape == (3, 4, 7)

    def test_identity_matrix(self):
        rng = np.random.default_rng(0)
        X = rng.standard_normal((3, 4, 5))
        I = np.eye(4)
        Y = nmode_product(X, I, mode=1)
        np.testing.assert_allclose(Y, X)

    def test_against_matricization(self):
        """Compare n-mode product with explicit matricize-multiply-fold."""
        rng = np.random.default_rng(123)
        X = rng.standard_normal((3, 4, 5))
        U = rng.standard_normal((2, 4))
        Y = nmode_product(X, U, mode=1)

        # Explicit: unfold, multiply, refold
        M = matricize(X, 1)  # 4 x 15
        result_mat = U @ M  # 2 x 15
        # Refold: shape should be (2, 3, 5) but with mode-1 first
        result = np.reshape(result_mat, (2, 3, 5))
        result = np.transpose(result, (1, 0, 2))
        np.testing.assert_allclose(Y, result)

    def test_dimension_mismatch(self):
        X = np.zeros((3, 4, 5))
        U = np.zeros((2, 3))  # wrong: should be (2, 4) for mode=1
        with pytest.raises(ValueError, match="columns"):
            nmode_product(X, U, mode=1)


class TestOuterProduct:
    """Outer product of vectors."""

    def test_two_vectors(self):
        a = np.array([1.0, 2.0])
        b = np.array([3.0, 4.0, 5.0])
        result = outer_product(a, b)
        assert result.shape == (2, 3)
        np.testing.assert_allclose(result, np.outer(a, b))

    def test_three_vectors(self):
        a = np.array([1.0, 2.0])
        b = np.array([3.0, 4.0])
        c = np.array([5.0, 6.0, 7.0])
        result = outer_product(a, b, c)
        assert result.shape == (2, 2, 3)
        # Check specific element: result[i,j,k] = a[i]*b[j]*c[k]
        assert result[1, 0, 2] == pytest.approx(2.0 * 3.0 * 7.0)

    def test_four_vectors(self):
        vecs = [np.ones(n) for n in (2, 3, 4, 5)]
        result = outer_product(*vecs)
        assert result.shape == (2, 3, 4, 5)
        np.testing.assert_allclose(result, np.ones((2, 3, 4, 5)))

    def test_single_vector_raises(self):
        with pytest.raises(ValueError, match="at least 2"):
            outer_product(np.array([1, 2, 3]))

    def test_non_1d_raises(self):
        with pytest.raises(ValueError, match="1-D"):
            outer_product(np.ones((2, 3)), np.ones(4))


class TestMultiIndex:
    """Linear to multi-index conversion."""

    def test_first_element(self):
        assert multi_index((3, 4, 5), 0) == (0, 0, 0)

    def test_last_element(self):
        assert multi_index((3, 4, 5), 59) == (2, 3, 4)

    def test_middle_element(self):
        # For shape (3,4,5), index 1 in C-order is (0,0,1)
        assert multi_index((3, 4, 5), 1) == (0, 0, 1)
