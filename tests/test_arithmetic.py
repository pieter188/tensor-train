"""Tests for tensortrain.arithmetic."""

import numpy as np
import pytest

from tensortrain import TensorTrain
from tensortrain.arithmetic import add, dot, sub


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
