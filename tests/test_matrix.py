"""Tests for TTMatrix class and TTm operations."""

import numpy as np
import pytest

from tensortrain import TensorTrain, TTMatrix
from tensortrain.arithmetic import add_ttm, matmat, matvec
from tensortrain.convert import transpose_ttm, tt_to_ttm, ttm_to_full, ttm_to_tt
from tensortrain.decompose import matrix_to_ttm


class TestTTMatrixConstruction:
    """TTMatrix construction and properties."""

    def test_basic(self):
        cores = [np.ones((1, 2, 3, 2)), np.ones((2, 4, 5, 1))]
        ttm = TTMatrix(cores)
        assert ttm.ndim == 2
        assert ttm.row_shape == (2, 4)
        assert ttm.col_shape == (3, 5)
        assert ttm.shape == (8, 15)
        assert ttm.ranks == (1, 2, 1)

    def test_eye(self):
        ttm = TTMatrix.eye((3, 4))
        M = ttm.full()
        np.testing.assert_allclose(M, np.eye(12))

    def test_eye_3modes(self):
        ttm = TTMatrix.eye((2, 3, 2))
        M = ttm.full()
        np.testing.assert_allclose(M, np.eye(12))

    def test_repr(self):
        ttm = TTMatrix.eye((2, 3))
        r = repr(ttm)
        assert "TTMatrix" in r
        assert "row_shape=(2, 3)" in r


class TestMatrixToTTm:
    """Decomposition of dense matrices into TTm."""

    def test_round_trip(self, rng):
        M = rng.standard_normal((6, 10))
        ttm = matrix_to_ttm(M, row_shape=(2, 3), col_shape=(2, 5))
        np.testing.assert_allclose(ttm.full(), M, atol=1e-12)

    def test_round_trip_eps(self, rng):
        M = rng.standard_normal((12, 15))
        eps = 1e-6
        ttm = matrix_to_ttm(M, row_shape=(3, 4), col_shape=(3, 5), eps=eps)
        rel_err = np.linalg.norm(ttm.full() - M) / np.linalg.norm(M)
        assert rel_err <= eps

    def test_identity(self):
        I = np.eye(12)
        ttm = matrix_to_ttm(I, row_shape=(3, 4), col_shape=(3, 4), eps=1e-12)
        np.testing.assert_allclose(ttm.full(), I, atol=1e-10)

    def test_from_matrix_classmethod(self, rng):
        M = rng.standard_normal((6, 10))
        ttm = TTMatrix.from_matrix(M, row_shape=(2, 3), col_shape=(2, 5))
        np.testing.assert_allclose(ttm.full(), M, atol=1e-12)


class TestTTmConversions:
    """TT <-> TTm conversions."""

    def test_tt_to_ttm_to_tt_roundtrip(self, rng):
        tt = TensorTrain.random((6, 20, 10), (2, 3), rng=rng)
        row_shape = (2, 4, 2)
        col_shape = (3, 5, 5)
        ttm = tt_to_ttm(tt, row_shape, col_shape)
        tt2 = ttm_to_tt(ttm)
        np.testing.assert_allclose(tt2.full(), tt.full(), atol=1e-12)

    def test_transpose(self, rng):
        M = rng.standard_normal((6, 10))
        ttm = matrix_to_ttm(M, row_shape=(2, 3), col_shape=(2, 5))
        ttm_t = transpose_ttm(ttm)
        np.testing.assert_allclose(ttm_t.full(), M.T, atol=1e-12)

    def test_transpose_property(self, rng):
        M = rng.standard_normal((6, 10))
        ttm = matrix_to_ttm(M, row_shape=(2, 3), col_shape=(2, 5))
        np.testing.assert_allclose(ttm.T.full(), M.T, atol=1e-12)


class TestMatvec:
    """TTm @ TT matrix-vector product."""

    def test_identity_matvec(self, rng):
        """I @ x == x."""
        x_data = rng.standard_normal(12)
        x_tt = TensorTrain.from_tensor(x_data.reshape(3, 4))
        I_ttm = TTMatrix.eye((3, 4))
        result = matvec(I_ttm, x_tt)
        np.testing.assert_allclose(result.full(), x_tt.full(), atol=1e-12)

    def test_matvec_against_dense(self, rng):
        """A @ x in TT format matches dense computation."""
        M = rng.standard_normal((6, 10))
        x = rng.standard_normal(10)
        expected = M @ x

        A_ttm = matrix_to_ttm(M, row_shape=(2, 3), col_shape=(2, 5))
        x_tt = TensorTrain.from_tensor(x.reshape(2, 5))
        result = matvec(A_ttm, x_tt)

        np.testing.assert_allclose(result.full().ravel(), expected, atol=1e-10)

    def test_matmul_operator(self, rng):
        """TTMatrix @ TensorTrain via __matmul__."""
        M = rng.standard_normal((6, 10))
        x = rng.standard_normal(10)

        A_ttm = matrix_to_ttm(M, row_shape=(2, 3), col_shape=(2, 5))
        x_tt = TensorTrain.from_tensor(x.reshape(2, 5))
        result = A_ttm @ x_tt

        np.testing.assert_allclose(result.full().ravel(), M @ x, atol=1e-10)

    def test_incompatible_shapes(self):
        A = TTMatrix.eye((3, 4))
        x = TensorTrain.ones((2, 3))
        with pytest.raises(ValueError, match="does not match"):
            matvec(A, x)


class TestMatmat:
    """TTm @ TTm matrix-matrix product."""

    def test_identity_matmat(self, rng):
        """I @ A == A."""
        M = rng.standard_normal((6, 10))
        A_ttm = matrix_to_ttm(M, row_shape=(2, 3), col_shape=(2, 5))
        I_ttm = TTMatrix.eye((2, 3))
        result = matmat(I_ttm, A_ttm)
        np.testing.assert_allclose(result.full(), M, atol=1e-10)

    def test_matmat_against_dense(self, rng):
        """A @ B in TTm format matches dense computation."""
        A = rng.standard_normal((6, 10))
        B = rng.standard_normal((10, 15))
        expected = A @ B

        A_ttm = matrix_to_ttm(A, row_shape=(2, 3), col_shape=(2, 5))
        B_ttm = matrix_to_ttm(B, row_shape=(2, 5), col_shape=(3, 5))
        result = matmat(A_ttm, B_ttm)

        np.testing.assert_allclose(result.full(), expected, atol=1e-9)

    def test_matmul_operator(self, rng):
        """TTMatrix @ TTMatrix via __matmul__."""
        A = rng.standard_normal((6, 10))
        B = rng.standard_normal((10, 15))

        A_ttm = matrix_to_ttm(A, row_shape=(2, 3), col_shape=(2, 5))
        B_ttm = matrix_to_ttm(B, row_shape=(2, 5), col_shape=(3, 5))
        result = A_ttm @ B_ttm

        np.testing.assert_allclose(result.full(), A @ B, atol=1e-9)


class TestTTMatrixArithmetic:
    """Addition/subtraction of TTMatrices."""

    def test_add(self, rng):
        A = rng.standard_normal((6, 10))
        B = rng.standard_normal((6, 10))
        A_ttm = matrix_to_ttm(A, row_shape=(2, 3), col_shape=(2, 5))
        B_ttm = matrix_to_ttm(B, row_shape=(2, 3), col_shape=(2, 5))
        C_ttm = A_ttm + B_ttm
        np.testing.assert_allclose(C_ttm.full(), A + B, atol=1e-10)

    def test_sub(self, rng):
        A = rng.standard_normal((6, 10))
        B = rng.standard_normal((6, 10))
        A_ttm = matrix_to_ttm(A, row_shape=(2, 3), col_shape=(2, 5))
        B_ttm = matrix_to_ttm(B, row_shape=(2, 3), col_shape=(2, 5))
        C_ttm = A_ttm - B_ttm
        np.testing.assert_allclose(C_ttm.full(), A - B, atol=1e-10)

    def test_scalar_mul(self, rng):
        A = rng.standard_normal((6, 10))
        A_ttm = matrix_to_ttm(A, row_shape=(2, 3), col_shape=(2, 5))
        np.testing.assert_allclose((3.0 * A_ttm).full(), 3.0 * A, atol=1e-10)
