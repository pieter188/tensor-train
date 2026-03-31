"""Fast internal operations on raw core lists.

These functions bypass TensorTrain construction and validation for
performance-critical inner loops (e.g., Gram-Schmidt).  They operate
on plain ``list[np.ndarray]`` instead of TensorTrain objects.

Not part of the public API.
"""

from __future__ import annotations

import numpy as np
from numpy.linalg import qr, svd


def fast_dot(cores_a: list[np.ndarray], cores_b: list[np.ndarray]) -> float:
    """Inner product of two TTs given as raw core lists."""
    M = np.ones((1, 1))
    for ca, cb in zip(cores_a, cores_b):
        M = np.einsum("ab,aic,bid->cd", M, ca, cb)
    return float(M.ravel()[0])


def fast_sub_scaled(
    v_cores: list[np.ndarray],
    q_cores: list[np.ndarray],
    scalar: float,
) -> list[np.ndarray]:
    """Compute ``v - scalar * q`` as block-diagonal TT cores.

    Returns new core list without creating TensorTrain objects.
    Equivalent to ``sub(v, q * scalar)`` but with zero object overhead.
    """
    d = len(v_cores)
    out = []
    for k in range(d):
        ca = v_cores[k]
        cb = q_cores[k]
        ra_l, n_k, ra_r = ca.shape
        rb_l, _, rb_r = cb.shape

        if k == 0:
            # First core: concat along right rank; negate and scale q
            # ca shape (1, n, ra_r), cb shape (1, n, rb_r)
            new = np.empty((1, n_k, ra_r + rb_r))
            new[:, :, :ra_r] = ca
            new[:, :, ra_r:] = -scalar * cb
            out.append(new)
        elif k == d - 1:
            # Last core: concat along left rank
            new = np.empty((ra_l + rb_l, n_k, 1))
            new[:ra_l, :, :] = ca
            new[ra_l:, :, :] = cb
            out.append(new)
        else:
            # Middle: block diagonal
            new = np.zeros((ra_l + rb_l, n_k, ra_r + rb_r))
            new[:ra_l, :, :ra_r] = ca
            new[ra_l:, :, ra_r:] = cb
            out.append(new)
    return out


def fast_normalize(cores: list[np.ndarray]) -> tuple[list[np.ndarray], float]:
    """Normalize a TT by dividing the first core by the full norm.

    Returns (normalized_cores, norm). Modifies cores in-place.
    """
    norm = fast_dot(cores, cores)
    norm = np.sqrt(max(norm, 0.0))
    if norm > 0:
        cores[0] = cores[0] / norm
    return cores, norm


def fast_orthogonalize_and_normalize(
    cores: list[np.ndarray],
) -> list[np.ndarray]:
    """Left-QR sweep to site 0, then normalize. Returns new cores."""
    d = len(cores)
    out = [c.copy() for c in cores]

    # Right-to-left QR sweep (makes cores 1..d-1 right-orthogonal)
    for k in range(d - 1, 0, -1):
        r_left, n_k, r_right = out[k].shape
        mat = out[k].reshape(r_left, n_k * r_right)
        Q, R = qr(mat.T)
        r_new = Q.shape[1]
        out[k] = Q.T.reshape(r_new, n_k, r_right)
        # Absorb R^T into previous core
        out[k - 1] = np.einsum("ijk,kl->ijl", out[k - 1], R.T)

    # Normalize first core
    norm = np.linalg.norm(out[0])
    if norm > 0:
        out[0] = out[0] / norm
    return out


def fast_round(
    cores: list[np.ndarray],
    eps: float,
) -> list[np.ndarray]:
    """TT-rounding on raw core lists.

    Left-QR sweep then right-SVD sweep with epsilon truncation.
    """
    d = len(cores)
    out = [c.copy() for c in cores]

    # Left-to-right QR sweep
    for k in range(d - 1):
        r_left, n_k, r_right = out[k].shape
        Q, R = qr(out[k].reshape(r_left * n_k, r_right))
        r_new = Q.shape[1]
        out[k] = Q.reshape(r_left, n_k, r_new)
        out[k + 1] = np.einsum("ij,jkl->ikl", R, out[k + 1])

    # Compute delta from norm (concentrated in last core)
    tt_norm = np.linalg.norm(out[-1].ravel())
    delta = (eps / np.sqrt(max(d - 1, 1))) * tt_norm
    delta_sq = delta * delta

    # Right-to-left SVD sweep with truncation
    for k in range(d - 1, 0, -1):
        r_left, n_k, r_right = out[k].shape
        mat = out[k].reshape(r_left, n_k * r_right)
        U, sigma, Vt = svd(mat, full_matrices=False)

        # Truncation rank
        r = len(sigma)
        if delta_sq > 0:
            tail_sq = np.cumsum(sigma[::-1] ** 2)[::-1]
            for r_try in range(1, len(sigma)):
                if tail_sq[r_try] <= delta_sq:
                    r = r_try
                    break

        r = max(r, 1)
        Vt_trunc = Vt[:r, :]
        out[k] = Vt_trunc.reshape(r, n_k, r_right)
        US = U[:, :r] * sigma[:r]
        out[k - 1] = np.einsum("ijk,kl->ijl", out[k - 1], US)

    return out


def max_rank(cores: list[np.ndarray]) -> int:
    """Maximum TT-rank of a core list."""
    r = 0
    for c in cores:
        r = max(r, c.shape[0], c.shape[2])
    return r


def batched_dot(
    q_cores: list[np.ndarray],
    v_list: list[list[np.ndarray]],
) -> np.ndarray:
    """Compute dot(q, v[j]) for all j simultaneously.

    Uses vectorized einsum when all vectors have the same core shapes.
    Falls back to individual dots for heterogeneous ranks.

    Parameters
    ----------
    q_cores : list of ndarray
        Single TT core list (the fixed vector).
    v_list : list of list of ndarray
        List of TT core lists to dot with q.

    Returns
    -------
    ndarray
        1-D array of length ``len(v_list)`` with the dot products.
    """
    n_batch = len(v_list)
    if n_batch == 0:
        return np.array([])

    d = len(q_cores)

    # Check if all vectors have the same core shapes (can batch)
    ref_shapes = tuple(v_list[0][k].shape for k in range(d))
    homogeneous = all(
        tuple(v_list[j][k].shape for k in range(d)) == ref_shapes
        for j in range(1, n_batch)
    )

    if homogeneous:
        # Fast path: stack and do batched einsum
        stacked = [np.array([v_list[j][k] for j in range(n_batch)]) for k in range(d)]
        M = np.ones((n_batch, 1, 1))
        for k in range(d):
            M = np.einsum("zab,aic,zbid->zcd", M, q_cores[k], stacked[k])
        return M.reshape(n_batch)
    else:
        # Fallback: individual dots
        result = np.empty(n_batch)
        for j in range(n_batch):
            result[j] = fast_dot(q_cores, v_list[j])
        return result


def batched_sub_scaled(
    v_list: list[list[np.ndarray]],
    q_cores: list[np.ndarray],
    scalars: np.ndarray,
) -> list[list[np.ndarray]]:
    """Compute ``v[j] - scalars[j] * q`` for all j simultaneously.

    Modifies v_list in-place and returns it.
    """
    d = len(q_cores)
    n_batch = len(v_list)

    for j in range(n_batch):
        s = scalars[j]
        if abs(s) < 1e-30:
            continue
        v_cores = v_list[j]
        for k in range(d):
            ca = v_cores[k]
            cb = q_cores[k]
            ra_l, n_k, ra_r = ca.shape
            rb_l, _, rb_r = cb.shape

            if k == 0:
                new = np.empty((1, n_k, ra_r + rb_r))
                new[:, :, :ra_r] = ca
                new[:, :, ra_r:] = -s * cb
            elif k == d - 1:
                new = np.empty((ra_l + rb_l, n_k, 1))
                new[:ra_l, :, :] = ca
                new[ra_l:, :, :] = cb
            else:
                new = np.zeros((ra_l + rb_l, n_k, ra_r + rb_r))
                new[:ra_l, :, :ra_r] = ca
                new[ra_l:, :, ra_r:] = cb
            v_cores[k] = new
    return v_list
