"""Internal validation helpers for tensor train cores and inputs."""

import numpy as np


def validate_cores(cores):
    """Validate a list of TT cores for consistent shapes.

    Each core must be a 3D ndarray of shape ``(r_{k-1}, n_k, r_k)`` with
    ``r_0 = r_d = 1`` and matching ranks between adjacent cores.

    Parameters
    ----------
    cores : list of ndarray
        Candidate TT cores.

    Returns
    -------
    list of ndarray
        The validated cores (converted to float64 contiguous arrays).

    Raises
    ------
    ValueError
        If the cores violate any structural constraint.
    """
    if not cores:
        raise ValueError("cores must be a non-empty list of ndarrays.")

    out = []
    for k, core in enumerate(cores):
        core = np.asarray(core, dtype=np.float64)
        if core.ndim != 3:
            raise ValueError(
                f"Core {k} must be 3-dimensional, got ndim={core.ndim}."
            )
        out.append(np.ascontiguousarray(core))

    # Boundary ranks
    if out[0].shape[0] != 1:
        raise ValueError(
            f"First core must have left rank 1, got {out[0].shape[0]}."
        )
    if out[-1].shape[2] != 1:
        raise ValueError(
            f"Last core must have right rank 1, got {out[-1].shape[2]}."
        )

    # Adjacent rank compatibility
    for k in range(len(out) - 1):
        r_right = out[k].shape[2]
        r_left_next = out[k + 1].shape[0]
        if r_right != r_left_next:
            raise ValueError(
                f"Rank mismatch between core {k} (right rank {r_right}) "
                f"and core {k + 1} (left rank {r_left_next})."
            )

    return out


def validate_ttm_cores(cores):
    """Validate a list of TT-matrix cores for consistent shapes.

    Each core must be a 4D ndarray of shape ``(r_{k-1}, m_k, n_k, r_k)``.

    Parameters
    ----------
    cores : list of ndarray
        Candidate TTm cores.

    Returns
    -------
    list of ndarray
        The validated cores (converted to float64 contiguous arrays).

    Raises
    ------
    ValueError
        If the cores violate any structural constraint.
    """
    if not cores:
        raise ValueError("cores must be a non-empty list of ndarrays.")

    out = []
    for k, core in enumerate(cores):
        core = np.asarray(core, dtype=np.float64)
        if core.ndim != 4:
            raise ValueError(
                f"Core {k} must be 4-dimensional, got ndim={core.ndim}."
            )
        out.append(np.ascontiguousarray(core))

    if out[0].shape[0] != 1:
        raise ValueError(
            f"First core must have left rank 1, got {out[0].shape[0]}."
        )
    if out[-1].shape[3] != 1:
        raise ValueError(
            f"Last core must have right rank 1, got {out[-1].shape[3]}."
        )

    for k in range(len(out) - 1):
        r_right = out[k].shape[3]
        r_left_next = out[k + 1].shape[0]
        if r_right != r_left_next:
            raise ValueError(
                f"Rank mismatch between core {k} (right rank {r_right}) "
                f"and core {k + 1} (left rank {r_left_next})."
            )

    return out


def check_compatible(a, b):
    """Check that two TensorTrain objects have the same mode sizes.

    Parameters
    ----------
    a, b : TensorTrain
        The two tensor trains to compare.

    Raises
    ------
    ValueError
        If shapes do not match.
    """
    if a.shape != b.shape:
        raise ValueError(
            f"Incompatible shapes: {a.shape} vs {b.shape}."
        )
