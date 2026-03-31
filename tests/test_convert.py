"""Tests for tensortrain.convert."""

import numpy as np
import pytest

from tensortrain import TensorTrain, TTMatrix
from tensortrain.convert import (
    combine_tt_vectors,
    extract_column,
    matrix_to_tensor,
    tensor_to_matrix,
    vector_to_tensor,
)
from tensortrain.decompose import matrix_to_ttm


class TestMatrixToTensor:
    def test_basic(self):
        M = np.arange(24.0)
        T = matrix_to_tensor(M, (2, 3, 4))
        assert T.shape == (2, 3, 4)
        assert T[0, 0, 0] == 0.0

    def test_2d(self):
        M = np.arange(12.0).reshape(3, 4)
        T = matrix_to_tensor(M, (3, 4))
        np.testing.assert_array_equal(T, M)


class TestTensorToMatrix:
    def test_basic(self):
        T = np.arange(24.0).reshape(2, 3, 4)
        M = tensor_to_matrix(T, row_modes=(0, 1), col_modes=(2,))
        assert M.shape == (6, 4)

    def test_single_mode_rows(self):
        T = np.arange(24.0).reshape(2, 3, 4)
        M = tensor_to_matrix(T, row_modes=(1,), col_modes=(0, 2))
        assert M.shape == (3, 8)


class TestVectorToTensor:
    def test_basic(self):
        v = np.arange(24.0)
        T = vector_to_tensor(v, (2, 3, 4))
        assert T.shape == (2, 3, 4)
        assert T[0, 0, 0] == 0.0


class TestExtractColumn:
    """Extract columns from TTMatrix."""

    def test_identity_columns(self):
        """Columns of identity matrix are standard basis vectors."""
        ttm = TTMatrix.eye((3, 4))
        for i in range(12):
            col = extract_column(ttm, i)
            expected = np.zeros(12)
            expected[i] = 1.0
            np.testing.assert_allclose(col.full().ravel(), expected, atol=1e-14)

    def test_against_dense(self, rng):
        """Extract column matches dense matrix column."""
        M = rng.standard_normal((6, 10))
        ttm = matrix_to_ttm(M, row_shape=(2, 3), col_shape=(2, 5))
        for i in [0, 3, 7, 9]:
            col = extract_column(ttm, i)
            np.testing.assert_allclose(col.full().ravel(), M[:, i], atol=1e-12)

    def test_method_delegation(self, rng):
        """TTMatrix.column() works."""
        M = rng.standard_normal((6, 10))
        ttm = matrix_to_ttm(M, row_shape=(2, 3), col_shape=(2, 5))
        col = ttm.column(4)
        np.testing.assert_allclose(col.full().ravel(), M[:, 4], atol=1e-12)

    def test_shape(self, rng):
        M = rng.standard_normal((6, 10))
        ttm = matrix_to_ttm(M, row_shape=(2, 3), col_shape=(2, 5))
        col = extract_column(ttm, 0)
        assert col.shape == (2, 3)

    def test_invalid_index(self):
        ttm = TTMatrix.eye((3, 4))
        with pytest.raises(ValueError, match="col_index"):
            extract_column(ttm, 12)


class TestCombineTTVectors:
    """Combine TT vectors into TTMatrix."""

    def test_round_trip_identity(self):
        """Extract all columns then combine should recover the matrix."""
        ttm = TTMatrix.eye((3, 4))
        cols = [extract_column(ttm, i) for i in range(12)]
        recovered = combine_tt_vectors(cols, col_shape=(3, 4), eps=1e-10, max_rank=20)
        np.testing.assert_allclose(recovered.full(), np.eye(12), atol=1e-6)

    def test_round_trip_random(self, rng):
        """Extract + combine recovers a random matrix (approximately)."""
        M = rng.standard_normal((6, 10))
        ttm = matrix_to_ttm(M, row_shape=(2, 3), col_shape=(2, 5))
        cols = [extract_column(ttm, i) for i in range(10)]
        recovered = combine_tt_vectors(cols, col_shape=(2, 5), eps=1e-8, max_rank=20)
        # Lossy due to rounding during combination
        np.testing.assert_allclose(recovered.full(), M, atol=0.5)

    def test_wrong_number_of_vectors(self):
        vecs = [TensorTrain.ones((2, 3))] * 5
        with pytest.raises(ValueError, match="Expected 6"):
            combine_tt_vectors(vecs, col_shape=(2, 3))
