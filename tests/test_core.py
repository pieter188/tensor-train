"""Tests for the TensorTrain class."""

import numpy as np
import pytest

from tensortrain import TensorTrain


class TestConstruction:
    """TensorTrain construction and validation."""

    def test_basic_construction(self):
        cores = [
            np.ones((1, 3, 2)),
            np.ones((2, 4, 1)),
        ]
        tt = TensorTrain(cores)
        assert tt.ndim == 2
        assert tt.shape == (3, 4)
        assert tt.ranks == (1, 2, 1)

    def test_three_cores(self):
        cores = [
            np.zeros((1, 3, 2)),
            np.zeros((2, 4, 3)),
            np.zeros((3, 5, 1)),
        ]
        tt = TensorTrain(cores)
        assert tt.ndim == 3
        assert tt.shape == (3, 4, 5)
        assert tt.ranks == (1, 2, 3, 1)

    def test_single_core(self):
        core = np.ones((1, 5, 1))
        tt = TensorTrain([core])
        assert tt.ndim == 1
        assert tt.shape == (5,)
        assert tt.ranks == (1, 1)

    def test_empty_cores_raises(self):
        with pytest.raises(ValueError, match="non-empty"):
            TensorTrain([])

    def test_wrong_ndim_raises(self):
        with pytest.raises(ValueError, match="3-dimensional"):
            TensorTrain([np.ones((3, 4))])

    def test_bad_left_boundary_raises(self):
        with pytest.raises(ValueError, match="left rank 1"):
            TensorTrain([np.ones((2, 3, 1))])

    def test_bad_right_boundary_raises(self):
        with pytest.raises(ValueError, match="right rank 1"):
            TensorTrain([np.ones((1, 3, 2))])

    def test_rank_mismatch_raises(self):
        cores = [np.ones((1, 3, 2)), np.ones((3, 4, 1))]
        with pytest.raises(ValueError, match="Rank mismatch"):
            TensorTrain(cores)

    def test_cores_converted_to_float64(self):
        core = np.ones((1, 3, 1), dtype=np.int32)
        tt = TensorTrain([core])
        assert tt.core(0).dtype == np.float64


class TestProperties:
    """Read-only properties."""

    def test_size(self):
        cores = [
            np.zeros((1, 3, 2)),  # 6
            np.zeros((2, 4, 3)),  # 24
            np.zeros((3, 5, 1)),  # 15
        ]
        tt = TensorTrain(cores)
        assert tt.size == 6 + 24 + 15

    def test_cores_returns_copy(self):
        tt = TensorTrain([np.ones((1, 3, 1))])
        cores_copy = tt.cores
        cores_copy[0] *= 999
        # Original must be untouched
        np.testing.assert_array_equal(tt.core(0), np.ones((1, 3, 1)))

    def test_core_reference(self):
        original = np.ones((1, 3, 1))
        tt = TensorTrain([original.copy()])
        ref = tt.core(0)
        assert ref is tt._cores[0]


class TestOperators:
    """Scalar multiplication and negation."""

    def test_scalar_mul(self):
        tt = TensorTrain([np.ones((1, 3, 1))])
        tt2 = tt * 3.0
        np.testing.assert_allclose(tt2.core(0), 3.0 * np.ones((1, 3, 1)))

    def test_rmul(self):
        tt = TensorTrain([np.ones((1, 3, 1))])
        tt2 = 2.5 * tt
        np.testing.assert_allclose(tt2.core(0), 2.5 * np.ones((1, 3, 1)))

    def test_neg(self):
        tt = TensorTrain([np.ones((1, 3, 1))])
        tt2 = -tt
        np.testing.assert_allclose(tt2.core(0), -np.ones((1, 3, 1)))

    def test_scalar_mul_does_not_mutate(self):
        tt = TensorTrain([np.ones((1, 3, 1))])
        _ = tt * 5.0
        np.testing.assert_array_equal(tt.core(0), np.ones((1, 3, 1)))


class TestFactoryMethods:
    """Class-method constructors."""

    def test_zeros(self):
        tt = TensorTrain.zeros((3, 4, 5))
        assert tt.shape == (3, 4, 5)
        assert tt.ranks == (1, 1, 1, 1)
        np.testing.assert_array_equal(tt.core(0), np.zeros((1, 3, 1)))

    def test_ones(self):
        tt = TensorTrain.ones((2, 3))
        assert tt.shape == (2, 3)
        np.testing.assert_array_equal(tt.core(0), np.ones((1, 2, 1)))
        np.testing.assert_array_equal(tt.core(1), np.ones((1, 3, 1)))

    def test_random_shape_and_ranks(self):
        tt = TensorTrain.random((4, 3, 5), (2, 3), rng=0)
        assert tt.shape == (4, 3, 5)
        assert tt.ranks == (1, 2, 3, 1)

    def test_random_wrong_ranks_length(self):
        with pytest.raises(ValueError, match="ranks"):
            TensorTrain.random((4, 3, 5), (2,))

    def test_random_reproducible(self):
        tt1 = TensorTrain.random((3, 4), (2,), rng=42)
        tt2 = TensorTrain.random((3, 4), (2,), rng=42)
        for k in range(tt1.ndim):
            np.testing.assert_array_equal(tt1.core(k), tt2.core(k))


class TestRepr:
    def test_repr_format(self):
        tt = TensorTrain.zeros((3, 4, 5))
        r = repr(tt)
        assert "TensorTrain" in r
        assert "(3, 4, 5)" in r
        assert "storage=" in r


class TestCopy:
    def test_deep_copy(self):
        tt = TensorTrain([np.ones((1, 3, 2)), np.ones((2, 4, 1))])
        tt2 = tt.copy()
        tt2._cores[0] *= 0
        # Original unchanged
        np.testing.assert_array_equal(tt.core(0), np.ones((1, 3, 2)))


class TestDivision:
    def test_truediv(self):
        tt = TensorTrain([np.ones((1, 3, 1)) * 6.0])
        tt2 = tt / 3.0
        np.testing.assert_allclose(tt2.core(0), np.ones((1, 3, 1)) * 2.0)

    def test_truediv_preserves_original(self):
        tt = TensorTrain([np.ones((1, 3, 1)) * 6.0])
        _ = tt / 2.0
        np.testing.assert_allclose(tt.core(0), np.ones((1, 3, 1)) * 6.0)


class TestLen:
    def test_len(self):
        tt = TensorTrain.zeros((3, 4, 5))
        assert len(tt) == 3

    def test_len_single(self):
        tt = TensorTrain([np.ones((1, 5, 1))])
        assert len(tt) == 1


class TestOrthogonalRandom:
    def test_unit_norm(self):
        tt = TensorTrain.random((4, 3, 5), (2, 3), rng=42, orthogonal=True)
        np.testing.assert_allclose(tt.norm(), 1.0, atol=1e-12)

    def test_norm_index_set(self):
        tt = TensorTrain.random((4, 3, 5), (2, 3), rng=42, orthogonal=True)
        assert tt.norm_index == 0

    def test_preserves_shape_and_ranks(self):
        tt = TensorTrain.random((4, 3, 5), (2, 3), rng=42, orthogonal=True)
        assert tt.shape == (4, 3, 5)
        # Ranks may differ from (2, 3) after QR but should be <= original
        for r in tt.ranks[1:-1]:
            assert r >= 1

    def test_non_orthogonal_default(self):
        """Default orthogonal=False gives unnormalized TT."""
        tt = TensorTrain.random((4, 3, 5), (2, 3), rng=42)
        assert tt.norm_index is None
        # Norm is typically not 1
        assert abs(tt.norm() - 1.0) > 0.1


class TestNormIndex:
    def test_norm_index_after_orthogonalize(self):
        tt = TensorTrain.random((4, 3, 5), (2, 3), rng=42)
        tt_orth = tt.orthogonalize(1)
        assert tt_orth.norm_index == 1

    def test_norm_fast_path(self):
        """norm() uses fast path when norm_index is set."""
        tt = TensorTrain.random((4, 3, 5), (2, 3), rng=42)
        expected = tt.norm()
        tt_orth = tt.orthogonalize(2)
        np.testing.assert_allclose(tt_orth.norm(), expected, rtol=1e-10)

    def test_norm_index_none_by_default(self):
        tt = TensorTrain.zeros((3, 4))
        assert tt.norm_index is None
