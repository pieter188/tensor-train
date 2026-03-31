"""TTMatrix class — Tensor Train representation of a structured matrix."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from tensortrain._validation import validate_ttm_cores

if TYPE_CHECKING:
    from collections.abc import Sequence

    from tensortrain._core import TensorTrain


class TTMatrix:
    """Tensor Train Matrix (TTm) representation of a structured matrix.

    A matrix :math:`A \\in \\mathbb{R}^{M \\times N}` is represented in
    TT-matrix format by factoring its rows and columns into *d* mode
    pairs ``(m_k, n_k)`` such that ``M = m_1 * m_2 * ... * m_d`` and
    ``N = n_1 * n_2 * ... * n_d``.

    Each core *k* is a 4-D array of shape ``(r_{k-1}, m_k, n_k, r_k)``
    with boundary ranks ``r_0 = r_d = 1``.

    Parameters
    ----------
    cores : list of ndarray
        Each ``cores[k]`` has shape ``(r_{k-1}, m_k, n_k, r_k)``.

    Examples
    --------
    >>> import numpy as np
    >>> from tensortrain import TTMatrix
    >>> cores = [np.random.randn(1, 2, 3, 2), np.random.randn(2, 4, 5, 1)]
    >>> ttm = TTMatrix(cores)
    >>> ttm.row_shape
    (2, 4)
    >>> ttm.col_shape
    (3, 5)
    >>> ttm.shape
    (8, 15)
    """

    def __init__(self, cores: Sequence[np.ndarray]) -> None:
        self._cores = validate_ttm_cores(list(cores))

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def ndim(self) -> int:
        """Number of mode pairs."""
        return len(self._cores)

    @property
    def row_shape(self) -> tuple[int, ...]:
        """Row mode sizes ``(m_1, ..., m_d)``."""
        return tuple(c.shape[1] for c in self._cores)

    @property
    def col_shape(self) -> tuple[int, ...]:
        """Column mode sizes ``(n_1, ..., n_d)``."""
        return tuple(c.shape[2] for c in self._cores)

    @property
    def shape(self) -> tuple[int, int]:
        """Full matrix shape ``(M, N)``."""
        M = int(np.prod(self.row_shape))
        N = int(np.prod(self.col_shape))
        return (M, N)

    @property
    def ranks(self) -> tuple[int, ...]:
        """TT-ranks ``(r_0, r_1, ..., r_d)``."""
        return (1,) + tuple(c.shape[3] for c in self._cores)

    @property
    def size(self) -> int:
        """Total stored elements."""
        return sum(c.size for c in self._cores)

    @property
    def cores(self) -> list[np.ndarray]:
        """Return a *copy* of the cores list."""
        return [c.copy() for c in self._cores]

    def core(self, k: int) -> np.ndarray:
        """Reference to core *k* (0-indexed)."""
        return self._cores[k]

    # ------------------------------------------------------------------
    # Convenience methods
    # ------------------------------------------------------------------

    def full(self) -> np.ndarray:
        """Reconstruct the full dense matrix.

        Returns
        -------
        ndarray
            Dense matrix of shape :attr:`shape`.

        See Also
        --------
        tensortrain.convert.ttm_to_full
        """
        from tensortrain.convert import ttm_to_full

        return ttm_to_full(self)

    def transpose(self) -> TTMatrix:
        """Return the transpose (swap row and column modes in each core).

        Returns
        -------
        TTMatrix

        See Also
        --------
        tensortrain.convert.transpose_ttm
        """
        from tensortrain.convert import transpose_ttm

        return transpose_ttm(self)

    @property
    def T(self) -> TTMatrix:
        """Alias for :meth:`transpose`."""
        return self.transpose()

    def to_tt(self) -> TensorTrain:
        """Convert to a TensorTrain by merging row/col modes.

        Returns
        -------
        TensorTrain

        See Also
        --------
        tensortrain.convert.ttm_to_tt
        """
        from tensortrain.convert import ttm_to_tt

        return ttm_to_tt(self)

    def matvec(self, x: TensorTrain) -> TensorTrain:
        """Matrix-vector product ``A @ x`` in TT format.

        Parameters
        ----------
        x : TensorTrain
            Vector in TT format with ``x.shape == self.col_shape``.

        Returns
        -------
        TensorTrain

        See Also
        --------
        tensortrain.arithmetic.matvec
        """
        from tensortrain.arithmetic import matvec

        return matvec(self, x)

    def matmat(self, other: TTMatrix) -> TTMatrix:
        """Matrix-matrix product ``A @ B`` in TTm format.

        Parameters
        ----------
        other : TTMatrix

        Returns
        -------
        TTMatrix

        See Also
        --------
        tensortrain.arithmetic.matmat
        """
        from tensortrain.arithmetic import matmat

        return matmat(self, other)

    def copy(self) -> TTMatrix:
        """Deep copy."""
        return TTMatrix([c.copy() for c in self._cores])

    # ------------------------------------------------------------------
    # Operator overloads
    # ------------------------------------------------------------------

    def __matmul__(self, other):
        from tensortrain._core import TensorTrain

        if isinstance(other, TensorTrain):
            return self.matvec(other)
        if isinstance(other, TTMatrix):
            return self.matmat(other)
        return NotImplemented

    def __add__(self, other: TTMatrix) -> TTMatrix:
        if not isinstance(other, TTMatrix):
            return NotImplemented
        from tensortrain.arithmetic import add_ttm

        return add_ttm(self, other)

    def __sub__(self, other: TTMatrix) -> TTMatrix:
        if not isinstance(other, TTMatrix):
            return NotImplemented
        from tensortrain.arithmetic import sub_ttm

        return sub_ttm(self, other)

    def __neg__(self) -> TTMatrix:
        cores = [c.copy() for c in self._cores]
        cores[0] = -cores[0]
        return TTMatrix(cores)

    def __mul__(self, scalar: float) -> TTMatrix:
        scalar = float(scalar)
        cores = [c.copy() for c in self._cores]
        cores[0] = cores[0] * scalar
        return TTMatrix(cores)

    def __rmul__(self, scalar: float) -> TTMatrix:
        return self.__mul__(scalar)

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_matrix(
        cls,
        matrix: np.ndarray,
        row_shape: Sequence[int],
        col_shape: Sequence[int],
        eps: float | None = None,
        max_ranks: Sequence[int] | None = None,
    ) -> TTMatrix:
        """Decompose a dense matrix into TTm format.

        Parameters
        ----------
        matrix : ndarray
            A 2-D array of shape ``(prod(row_shape), prod(col_shape))``.
        row_shape : sequence of int
            Row mode sizes ``(m_1, ..., m_d)``.
        col_shape : sequence of int
            Column mode sizes ``(n_1, ..., n_d)``.
        eps : float, optional
            Relative Frobenius-norm tolerance.
        max_ranks : sequence of int, optional
            Maximum TT-ranks.

        Returns
        -------
        TTMatrix

        See Also
        --------
        tensortrain.decompose.matrix_to_ttm
        """
        from tensortrain.decompose import matrix_to_ttm

        return matrix_to_ttm(matrix, row_shape, col_shape, eps=eps, max_ranks=max_ranks)

    @classmethod
    def eye(cls, mode_shape: Sequence[int]) -> TTMatrix:
        """TTm identity matrix with ``row_shape == col_shape == mode_shape``.

        Parameters
        ----------
        mode_shape : sequence of int
            Mode sizes ``(n_1, ..., n_d)``.

        Returns
        -------
        TTMatrix
            Represents the identity matrix of size
            ``prod(mode_shape) x prod(mode_shape)``.
        """
        cores = []
        for n in mode_shape:
            core = np.zeros((1, n, n, 1))
            for i in range(n):
                core[0, i, i, 0] = 1.0
            cores.append(core)
        return cls(cores)

    # ------------------------------------------------------------------
    # Repr
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        return (
            f"TTMatrix(row_shape={self.row_shape}, col_shape={self.col_shape}, "
            f"ranks={self.ranks}, storage={self.size})"
        )
