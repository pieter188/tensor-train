"""Tests for tensortrain.arithmetic."""

import numpy as np
import pytest

from tensortrain import TensorTrain, TTMatrix
from tensortrain.arithmetic import add, concat_ttm, dot, sub
from tensortrain.convert import extract_column
from tensortrain.decompose import matrix_to_ttm


class TestAdd:
    """TT addition."""

    def test_add_equals_full_add(self, rng):
        """full(A + B) == full(A) + full(B)."""
        a = TensorTrain.random((4, 3, 5), (2, 3), rng=0)
        b = TensorTrain.random((4, 3, 5), (2, 2), rng=1)
        c = add(a, b)
        np.testing.assert_allclose(c.full(), a.full() + b.full(), atol=1e-12)

    def test_add_ranks(self):
        """Addition sums the internal ranks."""
        a = TensorTrain.random((4, 3, 5), (2, 3), rng=0)
        b = TensorTrain.random((4, 3, 5), (1, 2), rng=1)
        c = add(a, b)
        # Internal ranks should be sum: (2+1, 3+2) = (3, 5)
        assert c.ranks == (1, 3, 5, 1)

    def test_add_operator(self, rng):
        """a + b uses __add__."""
        a = TensorTrain.random((3, 4), (2,), rng=0)
        b = TensorTrain.random((3, 4), (3,), rng=1)
        c = a + b
        np.testing.assert_allclose(c.full(), a.full() + b.full(), atol=1e-12)

    def test_incompatible_shapes(self):
        a = TensorTrain.zeros((3, 4))
        b = TensorTrain.zeros((3, 5))
        with pytest.raises(ValueError, match="Incompatible"):
            add(a, b)

    def test_add_2d(self):
        """Works for 2-core TTs."""
        a = TensorTrain.random((5, 7), (3,), rng=0)
        b = TensorTrain.random((5, 7), (2,), rng=1)
        c = a + b
        np.testing.assert_allclose(c.full(), a.full() + b.full(), atol=1e-12)


class TestSub:
    """TT subtraction."""

    def test_sub_equals_full_sub(self):
        a = TensorTrain.random((4, 3, 5), (2, 3), rng=0)
        b = TensorTrain.random((4, 3, 5), (2, 2), rng=1)
        c = sub(a, b)
        np.testing.assert_allclose(c.full(), a.full() - b.full(), atol=1e-12)

    def test_sub_operator(self):
        a = TensorTrain.random((3, 4), (2,), rng=0)
        b = TensorTrain.random((3, 4), (3,), rng=1)
        c = a - b
        np.testing.assert_allclose(c.full(), a.full() - b.full(), atol=1e-12)

    def test_self_sub_is_zero(self):
        a = TensorTrain.random((3, 4, 5), (2, 3), rng=42)
        c = a - a
        np.testing.assert_allclose(c.full(), np.zeros((3, 4, 5)), atol=1e-12)


class TestDot:
    """TT inner product."""

    def test_dot_equals_full_dot(self):
        a = TensorTrain.random((4, 3, 5), (2, 3), rng=0)
        b = TensorTrain.random((4, 3, 5), (2, 2), rng=1)
        result = dot(a, b)
        expected = np.sum(a.full() * b.full())
        np.testing.assert_allclose(result, expected, rtol=1e-10)

    def test_norm_via_dot(self):
        """<a, a> should equal ||a||^2."""
        a = TensorTrain.random((4, 3, 5), (2, 3), rng=42)
        result = dot(a, a)
        expected = np.linalg.norm(a.full()) ** 2
        np.testing.assert_allclose(result, expected, rtol=1e-10)

    def test_norm_method(self):
        """TensorTrain.norm() should work."""
        a = TensorTrain.random((4, 3, 5), (2, 3), rng=42)
        np.testing.assert_allclose(
            a.norm(), np.linalg.norm(a.full()), rtol=1e-10
        )

    def test_dot_method(self):
        a = TensorTrain.random((3, 4), (2,), rng=0)
        b = TensorTrain.random((3, 4), (3,), rng=1)
        np.testing.assert_allclose(
            a.dot(b), np.sum(a.full() * b.full()), rtol=1e-10
        )

    def test_dot_2d(self):
        a = TensorTrain.random((5, 7), (3,), rng=0)
        b = TensorTrain.random((5, 7), (2,), rng=1)
        np.testing.assert_allclose(
            dot(a, b), np.sum(a.full() * b.full()), rtol=1e-10
        )

    def test_incompatible_shapes(self):
        a = TensorTrain.zeros((3, 4))
        b = TensorTrain.zeros((3, 5))
        with pytest.raises(ValueError, match="Incompatible"):
            dot(a, b)

    def test_dot_with_dense_array(self):
        """dot(tt, ndarray) should work."""
        a = TensorTrain.random((4, 3, 5), (2, 3), rng=0)
        b_dense = np.random.randn(4, 3, 5)
        result = dot(a, b_dense)
        expected = np.sum(a.full() * b_dense)
        np.testing.assert_allclose(result, expected, rtol=1e-10)

    def test_dot_dense_shape_mismatch(self):
        a = TensorTrain.random((3, 4), (2,), rng=0)
        with pytest.raises(ValueError, match="Shape mismatch"):
            dot(a, np.zeros((3, 5)))


class TestRoundingWithAdd:
    """Test that rounding works after addition (cross-phase test)."""

    def test_add_then_round(self):
        """Addition followed by rounding should preserve the tensor."""
        a = TensorTrain.random((4, 3, 5), (2, 3), rng=0)
        b = TensorTrain.random((4, 3, 5), (2, 2), rng=1)
        c = a + b
        c_rounded = c.round(eps=1e-10)
        np.testing.assert_allclose(
            c_rounded.full(), a.full() + b.full(), atol=1e-8
        )
        # Rounded should have smaller or equal ranks
        for r_orig, r_round in zip(c.ranks, c_rounded.ranks):
            assert r_round <= r_orig


class TestConcatTTm:
    """Block-column concatenation [A | B]."""

    def test_concat_shape(self):
        """Last col mode should be doubled."""
        I = TTMatrix.eye((3, 4))
        C = concat_ttm(I, I)
        assert C.row_shape == (3, 4)
        assert C.col_shape == (3, 8)

    def test_concat_columns_first_half_is_A(self, rng):
        """First N columns of [A|B] should be A's columns."""
        A = rng.standard_normal((6, 10))
        B = rng.standard_normal((6, 10))
        A_ttm = matrix_to_ttm(A, row_shape=(2, 3), col_shape=(2, 5))
        B_ttm = matrix_to_ttm(B, row_shape=(2, 3), col_shape=(2, 5))
        C_ttm = concat_ttm(A_ttm, B_ttm)

        # Columns of A: last mode index 0..4 (first half)
        # Columns of B: last mode index 5..9 (second half)
        # Column (j1, j2) in C: j2 < 5 → A column (j1, j2)
        #                        j2 >= 5 → B column (j1, j2-5)
        n_cols_A = 10  # prod(col_shape) of A
        col_shape_C = C_ttm.col_shape  # (2, 10)
        n_d = A_ttm.col_shape[-1]  # 5

        for flat_a in range(n_cols_A):
            # Map A's flat column index to multi-index in A's col_shape
            j1 = flat_a // 5
            j2 = flat_a % 5
            # In C, same (j1, j2) → flat = j1 * 10 + j2
            flat_c = j1 * 10 + j2
            col_c = extract_column(C_ttm, flat_c)
            col_a = extract_column(A_ttm, flat_a)
            np.testing.assert_allclose(
                col_c.full(), col_a.full(), atol=1e-10
            )

    def test_concat_columns_second_half_is_B(self, rng):
        """Second N columns of [A|B] should be B's columns."""
        A = rng.standard_normal((6, 10))
        B = rng.standard_normal((6, 10))
        A_ttm = matrix_to_ttm(A, row_shape=(2, 3), col_shape=(2, 5))
        B_ttm = matrix_to_ttm(B, row_shape=(2, 3), col_shape=(2, 5))
        C_ttm = concat_ttm(A_ttm, B_ttm)

        n_d = A_ttm.col_shape[-1]  # 5
        for flat_b in range(10):
            j1 = flat_b // 5
            j2 = flat_b % 5
            flat_c = j1 * 10 + (j2 + n_d)
            col_c = extract_column(C_ttm, flat_c)
            col_b = extract_column(B_ttm, flat_b)
            np.testing.assert_allclose(
                col_c.full(), col_b.full(), atol=1e-10
            )

    def test_mismatched_row_shape(self):
        A = TTMatrix.eye((3, 4))
        B = TTMatrix.eye((2, 6))
        with pytest.raises(ValueError, match="row_shape"):
            concat_ttm(A, B)
