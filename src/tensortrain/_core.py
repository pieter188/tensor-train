"""TensorTrain class — the central data structure of the library."""

from __future__ import annotations

import copy
from typing import TYPE_CHECKING

import numpy as np

from tensortrain._validation import validate_cores

if TYPE_CHECKING:
    from collections.abc import Sequence


class TensorTrain:
    """Tensor Train decomposition of a *d*-dimensional tensor.

    A tensor :math:`\\mathcal{X} \\in \\mathbb{R}^{I_1 \\times \\cdots \\times I_d}`
    is stored as *d* cores, where core *k* is a 3-D array of shape
    ``(r_{k-1}, I_k, r_k)`` with boundary ranks ``r_0 = r_d = 1``.

    The element at multi-index ``(i_1, ..., i_d)`` is recovered by the
    matrix product:

    .. math::

        \\mathcal{X}(i_1, \\ldots, i_d)
        = G^{(1)}[:, i_1, :] \\; G^{(2)}[:, i_2, :] \\; \\cdots \\; G^{(d)}[:, i_d, :]

    Parameters
    ----------
    cores : list of ndarray
        Each ``cores[k]`` has shape ``(r_{k-1}, I_k, r_k)``.

    Examples
    --------
    Build a simple rank-1 tensor train for a 3x4x2 tensor:

    >>> import numpy as np
    >>> from tensortrain import TensorTrain
    >>> cores = [
    ...     np.random.randn(1, 3, 1),
    ...     np.random.randn(1, 4, 1),
    ...     np.random.randn(1, 2, 1),
    ... ]
    >>> tt = TensorTrain(cores)
    >>> tt.shape
    (3, 4, 2)
    >>> tt.ranks
    (1, 1, 1, 1)
    """

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __init__(self, cores: Sequence[np.ndarray]) -> None:
        self._cores = validate_cores(list(cores))
        self._norm_index: int | None = None

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def ndim(self) -> int:
        """Number of tensor dimensions (modes)."""
        return len(self._cores)

    @property
    def shape(self) -> tuple[int, ...]:
        """Mode sizes ``(I_1, I_2, ..., I_d)``."""
        return tuple(c.shape[1] for c in self._cores)

    @property
    def ranks(self) -> tuple[int, ...]:
        """TT-ranks ``(r_0, r_1, ..., r_d)`` — length ``d + 1``."""
        return (1,) + tuple(c.shape[2] for c in self._cores)

    @property
    def size(self) -> int:
        """Total number of stored floating-point elements across all cores."""
        return sum(c.size for c in self._cores)

    @property
    def norm_index(self) -> int | None:
        """Core index where the norm is concentrated, or ``None`` if unknown.

        Set after :meth:`orthogonalize`.  When known, :meth:`norm` uses
        this for an O(1) computation instead of a full inner product.
        """
        return self._norm_index

    @property
    def cores(self) -> list[np.ndarray]:
        """Return a *copy* of the internal cores list.

        Modifying the returned list does not affect the TensorTrain.
        Use :meth:`core` for read access to individual cores without copying.
        """
        return [c.copy() for c in self._cores]

    # ------------------------------------------------------------------
    # Core access
    # ------------------------------------------------------------------

    def core(self, k: int) -> np.ndarray:
        """Return a *reference* to core *k* (0-indexed).

        Parameters
        ----------
        k : int
            Core index in ``[0, d)``.

        Returns
        -------
        ndarray
            Shape ``(r_{k-1}, I_k, r_k)``.
        """
        return self._cores[k]

    # ------------------------------------------------------------------
    # Convenience methods (delegate to module functions)
    # ------------------------------------------------------------------

    def full(self) -> np.ndarray:
        """Reconstruct the full dense tensor.

        Returns
        -------
        ndarray
            Dense tensor of shape :attr:`shape`.

        See Also
        --------
        tensortrain.convert.tt_to_full
        """
        from tensortrain.convert import tt_to_full

        return tt_to_full(self)

    def round(
        self,
        eps: float | None = None,
        max_ranks: Sequence[int] | None = None,
    ) -> TensorTrain:
        """Return a rank-reduced copy of this tensor train.

        Parameters
        ----------
        eps : float, optional
            Relative Frobenius-norm tolerance.
        max_ranks : sequence of int, optional
            Maximum allowed ranks ``(r_1, ..., r_{d-1})``.

        Returns
        -------
        TensorTrain

        See Also
        --------
        tensortrain.rounding.tt_round
        """
        from tensortrain.rounding import tt_round

        return tt_round(self, eps=eps, max_ranks=max_ranks)

    def norm(self) -> float:
        """Frobenius norm computed efficiently in TT format.

        When :attr:`norm_index` is set (after orthogonalization), this
        is an O(1) operation.  Otherwise computes via inner product.

        Returns
        -------
        float

        See Also
        --------
        tensortrain.arithmetic.dot
        """
        if self._norm_index is not None:
            return float(np.linalg.norm(self._cores[self._norm_index]))

        from tensortrain.arithmetic import dot

        return float(np.sqrt(dot(self, self)))

    def dot(self, other: TensorTrain) -> float:
        """Inner product ``<self, other>``.

        Parameters
        ----------
        other : TensorTrain

        Returns
        -------
        float

        See Also
        --------
        tensortrain.arithmetic.dot
        """
        from tensortrain.arithmetic import dot

        return dot(self, other)

    def orthogonalize(self, site: int) -> TensorTrain:
        """Return a copy in site-*k* canonical form.

        Parameters
        ----------
        site : int
            The site index at which the norm is concentrated.

        Returns
        -------
        TensorTrain

        See Also
        --------
        tensortrain.canonical.orthogonalize
        """
        from tensortrain.canonical import orthogonalize

        return orthogonalize(self, site)

    def copy(self) -> TensorTrain:
        """Return a deep copy."""
        result = TensorTrain([c.copy() for c in self._cores])
        result._norm_index = self._norm_index
        return result

    # ------------------------------------------------------------------
    # Constructors (class methods)
    # ------------------------------------------------------------------

    @classmethod
    def from_tensor(
        cls,
        tensor: np.ndarray,
        eps: float | None = None,
        max_ranks: Sequence[int] | None = None,
    ) -> TensorTrain:
        """Decompose a dense tensor into TT format via TT-SVD.

        Parameters
        ----------
        tensor : ndarray
            A *d*-dimensional array.
        eps : float, optional
            Relative Frobenius-norm error tolerance.
        max_ranks : sequence of int, optional
            Maximum TT-ranks ``(r_1, ..., r_{d-1})``.

        Returns
        -------
        TensorTrain

        See Also
        --------
        tensortrain.decompose.tensor_to_tt
        """
        from tensortrain.decompose import tensor_to_tt

        return tensor_to_tt(tensor, eps=eps, max_ranks=max_ranks)

    @classmethod
    def zeros(cls, shape: Sequence[int]) -> TensorTrain:
        """TT representation of the zero tensor.

        Parameters
        ----------
        shape : sequence of int
            Mode sizes ``(I_1, ..., I_d)``.

        Returns
        -------
        TensorTrain
            All-zero tensor with rank-1 cores.
        """
        cores = [np.zeros((1, int(n), 1)) for n in shape]
        return cls(cores)

    @classmethod
    def ones(cls, shape: Sequence[int]) -> TensorTrain:
        """TT representation of the all-ones tensor.

        Parameters
        ----------
        shape : sequence of int
            Mode sizes ``(I_1, ..., I_d)``.

        Returns
        -------
        TensorTrain
            All-ones tensor with rank-1 cores.
        """
        cores = [np.ones((1, int(n), 1)) for n in shape]
        return cls(cores)

    @classmethod
    def from_vector(
        cls,
        vector: np.ndarray,
        eps: float | None = None,
        max_ranks: Sequence[int] | None = None,
    ) -> TensorTrain:
        """Decompose a vector into QTT (Quantized Tensor Train) format.

        Automatically factorizes the vector length into prime factors and
        reshapes the vector into a tensor before applying TT-SVD.

        Parameters
        ----------
        vector : ndarray
            A 1-D array.
        eps : float, optional
            Relative Frobenius-norm tolerance.
        max_ranks : sequence of int, optional
            Maximum TT-ranks.

        Returns
        -------
        TensorTrain

        See Also
        --------
        tensortrain.decompose.qtt
        """
        from tensortrain.decompose import qtt

        return qtt(vector, eps=eps, max_ranks=max_ranks)

    @classmethod
    def random(
        cls,
        shape: Sequence[int],
        ranks: Sequence[int],
        rng: np.random.Generator | int | None = None,
        orthogonal: bool = False,
    ) -> TensorTrain:
        """Random TT with specified mode sizes and ranks.

        Parameters
        ----------
        shape : sequence of int
            Mode sizes ``(I_1, ..., I_d)``.
        ranks : sequence of int
            Internal ranks ``(r_1, ..., r_{d-1})``.  Length must be ``d - 1``.
        rng : Generator, int, or None
            NumPy random generator or seed.
        orthogonal : bool, optional
            If True, orthogonalize cores and normalize to unit norm.

        Returns
        -------
        TensorTrain
        """
        shape = tuple(int(n) for n in shape)
        ranks = tuple(int(r) for r in ranks)
        d = len(shape)
        if len(ranks) != d - 1:
            raise ValueError(
                f"Expected {d - 1} ranks for {d} modes, got {len(ranks)}."
            )
        if isinstance(rng, (int, np.integer)):
            rng = np.random.default_rng(rng)
        elif rng is None:
            rng = np.random.default_rng()

        full_ranks = (1,) + ranks + (1,)
        cores = [
            rng.standard_normal((full_ranks[k], shape[k], full_ranks[k + 1]))
            for k in range(d)
        ]
        tt = cls(cores)
        if orthogonal:
            tt = tt.orthogonalize(0)
            # Normalize to unit Frobenius norm
            n = float(np.linalg.norm(tt._cores[0]))
            if n > 0:
                tt._cores[0] = tt._cores[0] / n
                tt._norm_index = 0
        return tt

    # ------------------------------------------------------------------
    # Operator overloads
    # ------------------------------------------------------------------

    def __add__(self, other: TensorTrain) -> TensorTrain:
        if not isinstance(other, TensorTrain):
            return NotImplemented
        from tensortrain.arithmetic import add

        return add(self, other)

    def __sub__(self, other: TensorTrain) -> TensorTrain:
        if not isinstance(other, TensorTrain):
            return NotImplemented
        from tensortrain.arithmetic import sub

        return sub(self, other)

    def __neg__(self) -> TensorTrain:
        cores = [c.copy() for c in self._cores]
        cores[0] = -cores[0]
        return TensorTrain(cores)

    def __mul__(self, scalar: float) -> TensorTrain:
        scalar = float(scalar)
        cores = [c.copy() for c in self._cores]
        cores[0] = cores[0] * scalar
        return TensorTrain(cores)

    def __rmul__(self, scalar: float) -> TensorTrain:
        return self.__mul__(scalar)

    def __truediv__(self, scalar: float) -> TensorTrain:
        return self * (1.0 / float(scalar))

    def __len__(self) -> int:
        return self.ndim

    # ------------------------------------------------------------------
    # Representation
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        return (
            f"TensorTrain(shape={self.shape}, ranks={self.ranks}, "
            f"storage={self.size})"
        )
