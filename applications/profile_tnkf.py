"""Profile the TNKF to identify bottlenecks.

Runs a single sparse frame and times each major step.
"""

import time
from pathlib import Path

import numpy as np

import tensortrain as tt
from tensortrain.arithmetic import add, concat_ttm, dot, matvec, sub
from tensortrain.canonical import orthogonalize
from tensortrain.convert import combine_tt_vectors, extract_column, ttm_to_tt, tt_to_ttm
from tensortrain.decompose import tensor_to_tt
from tensortrain.linalg import gram_schmidt_tt
from tensortrain.rounding import tt_round

from tnkf import (
    load_video_grayscale,
    mask_video,
    init_matrix_from_dense,
    init_covariance,
)


def timed(label):
    """Context manager that prints elapsed time."""
    class Timer:
        def __enter__(self):
            self.t0 = time.perf_counter()
            return self
        def __exit__(self, *args):
            elapsed = time.perf_counter() - self.t0
            self.elapsed = elapsed
            print(f"  {label}: {elapsed:.3f}s")
    return Timer()


def profile_one_frame():
    video = load_video_grayscale(
        Path(__file__).parent.parent / "matlab-og" / "original_video_TC.mp4",
        (9, 16),
    )
    masked, mask = mask_video(video, alpha=0.95, beta=0.01)

    Row = (3, 3, 4, 2, 2)
    Col = (3, 3, 4, 2, 2)
    h, w = 9, 16
    n = h * w  # 144
    n_cols = int(np.prod(Col))
    eps = 1e-4

    y = masked.reshape(n, -1)

    print("=" * 60)
    print("TNKF Profiling — single sparse frame")
    print("=" * 60)

    # --- Initialization ---
    with timed("Init C_ttm (identity)"):
        C_ttm = init_matrix_from_dense(np.eye(n), Row, Col, h, eps=1e-6)

    with timed("Init L_ttm (identity)"):
        L_ttm = init_matrix_from_dense(np.eye(n), Row, Col, h, eps=1e-6)

    with timed("Init Q_ttm (covariance)"):
        Q_ttm = init_covariance(masked, Row, Col, eps=1e-6)

    with timed("Init x_tt (first frame)"):
        x_tt = tensor_to_tt(y[:, 0].reshape(Row), eps=1e-6)

    # Use frame 3 (first sparse frame)
    k = 3
    observed = np.nonzero(y[:, k])[0]
    y_values = y[observed, k]
    print(f"\nFrame {k}: {len(observed)} observed pixels out of {n}")

    # ============================================================
    # TIME PROPAGATION
    # ============================================================
    print("\n--- TIME PROPAGATION ---")

    with timed("concat_ttm([L | Q])"):
        L_aug = concat_ttm(L_ttm, Q_ttm)

    with timed("Transpose L_aug"):
        L_aug_t = L_aug.T

    with timed("Canonical form (orthogonalize)"):
        L_aug_tt = ttm_to_tt(L_aug_t)
        L_aug_tt = orthogonalize(L_aug_tt, 0)
        L_aug_t = tt_to_ttm(L_aug_tt, L_aug_t.row_shape, L_aug_t.col_shape)

    with timed(f"Extract {n_cols} columns from L_aug_t"):
        tp_vectors = [extract_column(L_aug_t, i) for i in range(n_cols)]

    with timed(f"Gram-Schmidt on {n_cols} vectors (TIME PROP)"):
        q_vectors = gram_schmidt_tt(tp_vectors, eps=eps, max_rank=5)

    with timed(f"combine_tt_vectors ({n_cols} vectors -> TTm)"):
        Q_mat = combine_tt_vectors(q_vectors, Col, eps=eps, max_rank=5)

    with timed("R = Q^T @ L_aug_t (matmat)"):
        R_ttm = Q_mat.T @ L_aug_t

    with timed("Round R + transpose -> L_ttm"):
        R_tt = ttm_to_tt(R_ttm)
        R_tt = tt_round(R_tt, eps=eps)
        R_ttm = tt_to_ttm(R_tt, R_ttm.row_shape, R_ttm.col_shape)
        L_ttm = R_ttm.T

    # ============================================================
    # MEASUREMENT UPDATE (single pixel only)
    # ============================================================
    print(f"\n--- MEASUREMENT UPDATE (1 pixel out of {len(observed)}) ---")
    pixel_idx = observed[0]

    with timed("Extract column c from C_ttm"):
        c_tt = extract_column(C_ttm, pixel_idx)

    with timed("Residual v = y - <c, x>"):
        v = y_values[0] - dot(c_tt, x_tt)

    with timed("L^T @ c (matvec)"):
        Lt_c = matvec(L_ttm.T, c_tt)

    with timed("L @ (L^T @ c) (matvec)"):
        L_Lt_c = matvec(L_ttm, Lt_c)

    with timed("s = <c, L L^T c> (dot)"):
        s = dot(c_tt, L_Lt_c)

    with timed("K = L L^T c * (v/s) (scalar mul)"):
        K_tt = L_Lt_c * (v / s)

    with timed("x = x + K (add + orthogonalize)"):
        x_tt = add(x_tt, K_tt)
        x_tt = orthogonalize(x_tt, 0)

    # Covariance measurement update
    print("\n--- COVARIANCE MEASUREMENT UPDATE (1 pixel) ---")

    with timed("Transpose L_ttm"):
        L_t = L_ttm.T

    with timed(f"Extract {n_cols} columns from L^T"):
        mu_vectors = [Lt_c]
        for i in range(n_cols):
            mu_vectors.append(extract_column(L_t, i))

    with timed(f"Gram-Schmidt on {n_cols} vectors (MEAS UPDATE)"):
        q_mu = gram_schmidt_tt(mu_vectors[:n_cols], eps=eps, max_rank=5)

    with timed(f"combine_tt_vectors ({n_cols} vectors -> TTm)"):
        Q_mu_mat = combine_tt_vectors(q_mu, Col, eps=eps, max_rank=5)

    with timed("R2 = Q^T @ L_t (matmat, no augmentation)"):
        R2 = Q_mu_mat.T @ L_t
        R2_tt = ttm_to_tt(R2)
        R2_tt = tt_round(R2_tt, eps=eps)
        R2 = tt_to_ttm(R2_tt, R2.row_shape, R2.col_shape)

    with timed("Transpose R2 + shift + combine L_new"):
        zero_tt = tt.TensorTrain.zeros(Row)
        R2_t = R2.T
        r2t_vectors = []
        for i in range(1, n_cols):
            r2t_vectors.append(extract_column(R2_t, i))
        r2t_vectors.append(zero_tt)
        L_ttm = combine_tt_vectors(r2t_vectors, Col, eps=eps, max_rank=5)

    print("\n" + "=" * 60)
    print("SUMMARY: Estimated time for full frame with")
    print(f"  {len(observed)} observed pixels:")
    print(f"  Time prop happens once, meas update x{len(observed)}")
    print("=" * 60)


if __name__ == "__main__":
    profile_one_frame()
