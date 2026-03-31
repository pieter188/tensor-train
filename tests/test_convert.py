"""Tests for tensortrain.convert."""

import numpy as np

from tensortrain.convert import matrix_to_tensor, tensor_to_matrix, vector_to_tensor


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
