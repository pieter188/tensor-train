"""Tensor-Networked Square-Root Kalman Filter (TNKF) for video inpainting.

Recovers missing pixels in video frames using a Kalman filter where all
state-space matrices are stored in Tensor Train (TT) format.  This
dramatically reduces memory and computation for high-dimensional state
estimation.

Based on: MSc Thesis by Pieter van Klaveren — "Tensor Networks and
Tensor-Networked Kalman Filters"

Usage
-----
    python applications/tnkf.py

The script loads ``original_video_TC.mp4``, downscales to 9x16 grayscale,
masks 95% of pixels, and recovers them via the TNKF.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

import tensortrain as tt
from tensortrain.arithmetic import add, concat_ttm, dot, matvec, sub
from tensortrain.canonical import orthogonalize
from tensortrain.convert import combine_tt_vectors, extract_column, ttm_to_tt, tt_to_ttm
from tensortrain.decompose import matrix_to_ttm, tensor_to_tt
from tensortrain.linalg import gram_schmidt_tt
from tensortrain.rounding import tt_round


# ======================================================================
# Video preprocessing
# ======================================================================


def load_video_grayscale(path: str | Path, resolution: tuple[int, int] = (9, 16)) -> np.ndarray:
    """Load a video, resize, and convert to grayscale.

    Parameters
    ----------
    path : str or Path
        Path to video file.
    resolution : (height, width)
        Target resolution after downscaling.

    Returns
    -------
    ndarray
        Grayscale video of shape ``(height, width, num_frames)``, dtype float64,
        values in [0, 255].
    """
    import imageio.v3 as iio

    frames = []
    for frame in iio.imiter(str(path)):
        # frame shape: (H, W, 3) or (H, W) uint8
        if frame.ndim == 3 and frame.shape[2] >= 3:
            # RGB to grayscale: standard ITU-R 601 weights
            gray = (
                0.2989 * frame[:, :, 0].astype(np.float64)
                + 0.5870 * frame[:, :, 1].astype(np.float64)
                + 0.1140 * frame[:, :, 2].astype(np.float64)
            )
        else:
            gray = frame.astype(np.float64)

        # Resize using simple nearest-neighbor (avoid heavy dependencies)
        h, w = resolution
        H, W = gray.shape
        row_idx = (np.arange(h) * H / h).astype(int)
        col_idx = (np.arange(w) * W / w).astype(int)
        resized = gray[np.ix_(row_idx, col_idx)]
        frames.append(resized)

    video = np.stack(frames, axis=-1)  # (height, width, num_frames)
    return video


def mask_video(
    video: np.ndarray,
    alpha: float = 0.95,
    beta: float = 0.01,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Create sparse observations by randomly masking pixels.

    Parameters
    ----------
    video : ndarray
        Shape ``(height, width, num_frames)``.
    alpha : float
        Fraction of missing pixels per frame (0.95 = 95% missing).
    beta : float
        Fraction of initial frames that are fully observed.
    rng : Generator, optional
        Random generator for reproducibility.

    Returns
    -------
    masked : ndarray
        Same shape, with masked pixels set to 0.
    mask : ndarray (bool)
        True where pixel is observed.
    """
    if rng is None:
        rng = np.random.default_rng(42)

    h, w, T = video.shape
    n_pixels = h * w
    n_unmasked_frames = max(1, int(beta * T))

    mask = np.zeros((h, w, T), dtype=bool)

    # First frames: fully observed
    mask[:, :, :n_unmasked_frames] = True

    # Remaining frames: observe (1 - alpha) fraction
    n_observed = max(1, int((1 - alpha) * n_pixels))
    for t in range(n_unmasked_frames, T):
        observed_idx = rng.choice(n_pixels, size=n_observed, replace=False)
        flat_mask = np.zeros(n_pixels, dtype=bool)
        flat_mask[observed_idx] = True
        mask[:, :, t] = flat_mask.reshape(h, w)

    masked = video * mask
    return masked, mask


# ======================================================================
# State space initialization
# ======================================================================


def find_split_index(row_partition: tuple[int, ...], n_rows: int) -> int:
    """Find split point where prod(row_partition[:i]) == n_rows."""
    prod = 1
    for i, r in enumerate(row_partition):
        prod *= r
        if prod == n_rows:
            return i + 1
    raise ValueError(
        f"Cannot split partition {row_partition} to produce {n_rows}."
    )


def kronecker_ttm(a: tt.TTMatrix, b: tt.TTMatrix) -> tt.TTMatrix:
    """Kronecker product of two TTMatrices by concatenating cores.

    Parameters
    ----------
    a, b : TTMatrix
        The result represents ``kron(a.full(), b.full())``.

    Returns
    -------
    TTMatrix
    """
    cores = [c.copy() for c in a._cores] + [c.copy() for c in b._cores]
    return tt.TTMatrix(cores)


def init_matrix_from_dense(
    mat: np.ndarray,
    row_partition: tuple[int, ...],
    col_partition: tuple[int, ...],
    n_rows_first: int,
    eps: float = 1e-6,
) -> tt.TTMatrix:
    """Initialize a TTMatrix via Kronecker product A2 ⊗ A1.

    Splits the matrix into two blocks based on the image row/column
    structure and builds a Kronecker TTm.

    Parameters
    ----------
    mat : ndarray
        Square matrix of size ``n x n``.
    row_partition, col_partition : tuple of int
        TT mode partitions (must satisfy ``prod == n``).
    n_rows_first : int
        Number of image rows (height); used to split the partition.
    eps : float
        TT-SVD tolerance.

    Returns
    -------
    TTMatrix
    """
    n = mat.shape[0]
    split = find_split_index(row_partition, n_rows_first)

    # First block: rows of image (size n_rows_first)
    row2 = row_partition[:split]
    col2 = col_partition[:split]
    mat2 = mat[:n_rows_first, :n_rows_first]
    ttm2 = matrix_to_ttm(mat2, row_shape=row2, col_shape=col2, eps=eps)

    # Second block: columns of image (size n_cols)
    n_cols = n // n_rows_first
    row1 = row_partition[split:]
    col1 = col_partition[split:]
    mat1 = mat[:n_cols, :n_cols]
    ttm1 = matrix_to_ttm(mat1, row_shape=row1, col_shape=col1, eps=eps)

    return kronecker_ttm(ttm2, ttm1)


def init_covariance(
    data: np.ndarray,
    row_partition: tuple[int, ...],
    col_partition: tuple[int, ...],
    eps: float = 1e-6,
) -> tt.TTMatrix:
    """Build process noise covariance Q from spatial correlation.

    Computes pixel-wise correlation from the video data and builds a
    band-diagonal covariance matrix in TTm format via Kronecker structure.

    Parameters
    ----------
    data : ndarray
        Video data, shape ``(height, width, num_frames)``.
    row_partition, col_partition : tuple of int
    eps : float

    Returns
    -------
    TTMatrix
    """
    h, w, T = data.shape
    n = h * w

    # Compute correlation structure
    data_flat = data.reshape(n, T)
    # Use correlation of pixel time series
    if n <= 40 * 40:
        W = np.corrcoef(data_flat)
    else:
        # Random subset for large images
        rng = np.random.default_rng(0)
        idx = rng.choice(n, size=min(n, 1600), replace=False)
        W = np.corrcoef(data_flat[idx])

    # Compute mean correlation vs Manhattan distance
    max_dist = h + w - 2
    corr_by_dist = {}
    for i in range(min(W.shape[0], h)):
        for j in range(min(W.shape[0], w)):
            for ii in range(i, min(W.shape[0], h)):
                for jj in range(j, min(W.shape[0], w)):
                    d = abs(i - ii) + abs(j - jj)
                    if d not in corr_by_dist:
                        corr_by_dist[d] = []
                    idx1 = jj * h + ii if (jj * h + ii) < W.shape[0] else None
                    idx2 = j * h + i if (j * h + i) < W.shape[0] else None
                    if idx1 is not None and idx2 is not None:
                        corr_by_dist[d].append(W[idx1, idx2])

    corr_mean = {}
    for d, vals in corr_by_dist.items():
        corr_mean[d] = np.mean(vals)

    # Find bandwidth (max distance where correlation > 0.2)
    bandwidth = max(
        (d for d, c in corr_mean.items() if c > 0.2),
        default=1,
    )
    bandwidth = max(bandwidth, 1)

    # Build band-diagonal covariance matrices
    Q1 = np.eye(w)  # width dimension
    Q2 = np.eye(h)  # height dimension
    for i in range(1, bandwidth):
        decay = 1.0 - (i + 1) / 16.0
        if decay <= 0:
            break
        if i < w:
            Q1 += np.diag(np.full(w - i, decay), -i) + np.diag(np.full(w - i, decay), i)
        if i < h:
            Q2 += np.diag(np.full(h - i, decay), -i) + np.diag(np.full(h - i, decay), i)

    # Convert to TTm format
    split = find_split_index(row_partition, h)
    row_q2, col_q2 = row_partition[:split], col_partition[:split]
    row_q1, col_q1 = row_partition[split:], col_partition[split:]

    ttm_q2 = matrix_to_ttm(Q2, row_shape=row_q2, col_shape=col_q2, eps=eps)
    ttm_q1 = matrix_to_ttm(Q1, row_shape=row_q1, col_shape=col_q1, eps=eps)

    return kronecker_ttm(ttm_q2, ttm_q1)


# ======================================================================
# TNKF Algorithm
# ======================================================================


def tnkf(
    original_video: np.ndarray,
    masked_video: np.ndarray,
    observation_mask: np.ndarray,
    row_partition: tuple[int, ...],
    col_partition: tuple[int, ...],
    eps: float = 1e-5,
    max_rank_tp: int = 5,
    max_rank_mu: int = 5,
    max_frames: int | None = None,
    max_pixels_per_frame: int | None = None,
) -> np.ndarray:
    """Run the Tensor-Networked Square-Root Kalman Filter.

    Parameters
    ----------
    original_video : ndarray
        Shape ``(h, w, T)`` — used only for RMSE reporting.
    masked_video : ndarray
        Shape ``(h, w, T)`` — sparse observations.
    observation_mask : ndarray (bool)
        Shape ``(h, w, T)`` — True where observed.
    row_partition, col_partition : tuple of int
        TT mode partitions.
    eps : float
        TT-SVD and rounding tolerance.
    max_rank_tp : int
        Max rank threshold during time propagation QR.
    max_rank_mu : int
        Max rank threshold during measurement update QR.
    max_frames : int, optional
        Process only first N frames (for testing).

    Returns
    -------
    ndarray
        Estimated video, shape ``(h, w, T)``.
    """
    h, w, T = masked_video.shape
    n = h * w

    if max_frames is not None:
        T = min(T, max_frames)

    # Vectorize observations
    y = masked_video[:, :, :T].reshape(n, T)

    print(f"State dimension: {n}")
    print(f"Partition: Row={row_partition}, Col={col_partition}")
    print(f"Number of frames: {T}")
    print()

    # Initialize state space matrices
    print("Initializing state space matrices...")
    C_ttm = init_matrix_from_dense(
        np.eye(n), row_partition, col_partition, h, eps=1e-6
    )
    L_ttm = init_matrix_from_dense(
        np.eye(n), row_partition, col_partition, h, eps=1e-6
    )
    Q_ttm = init_covariance(
        masked_video[:, :, :T], row_partition, col_partition, eps=1e-6
    )

    # Initial state estimate = first frame
    x_est = np.zeros((n, T))
    x_est[:, 0] = y[:, 0]
    x_tt = tensor_to_tt(x_est[:, 0].reshape(row_partition), eps=1e-6)

    # Storage for state estimates
    estimates = [x_est[:, 0].copy()]

    print("Starting TNKF iterations...")
    for k in range(T - 1):
        # Find observed pixels in frame k+1
        observed = np.nonzero(y[:, k + 1])[0]
        y_values = y[observed, k + 1]

        if len(observed) == n:
            # All pixels observed — set state directly
            x_tt = tensor_to_tt(y[:, k + 1].reshape(row_partition), eps=eps)
            estimates.append(y[:, k + 1].copy())
            print(f"  Frame {k + 1}: fully observed")
            continue

        # ---- TIME PROPAGATION ----
        # x_{k+1|k} = x_k (identity state transition)
        # L_{k+1|k} from QR of [L | Q]^T

        import time as _time
        _t0 = _time.perf_counter()

        L_aug = concat_ttm(L_ttm, Q_ttm)
        L_aug_t = L_aug.T

        # Bring to canonical form
        L_aug_tt = ttm_to_tt(L_aug_t)
        L_aug_tt = orthogonalize(L_aug_tt, 0)
        L_aug_t = tt_to_ttm(L_aug_tt, L_aug_t.row_shape, L_aug_t.col_shape)

        # Extract columns for Gram-Schmidt QR
        n_cols = int(np.prod(col_partition))
        tp_vectors = [extract_column(L_aug_t, i) for i in range(n_cols)]

        print(f"    [TP] extract: {_time.perf_counter()-_t0:.1f}s, "
              f"L ranks={L_ttm.ranks}", flush=True)
        _t0 = _time.perf_counter()

        # Modified Gram-Schmidt
        q_vectors = gram_schmidt_tt(tp_vectors, eps=eps, max_rank=max_rank_tp)

        print(f"    [TP] GS: {_time.perf_counter()-_t0:.1f}s", flush=True)
        _t0 = _time.perf_counter()

        # Combine Q vectors into matrix
        Q_mat = combine_tt_vectors(q_vectors, col_partition, eps=eps, max_rank=max_rank_tp)

        # R = Q^T @ L_aug^T
        R_ttm = Q_mat.T @ L_aug_t
        R_tt = ttm_to_tt(R_ttm)
        R_tt = tt_round(R_tt, eps=eps)
        R_ttm = tt_to_ttm(R_tt, R_ttm.row_shape, R_ttm.col_shape)
        L_ttm = R_ttm.T

        print(f"    [TP] combine+R: {_time.perf_counter()-_t0:.1f}s, "
              f"L_new ranks={L_ttm.ranks}", flush=True)

        # ---- MEASUREMENT UPDATE (per observed pixel) ----
        max_pixels = max_pixels_per_frame
        pixels_to_process = observed[:max_pixels] if max_pixels else observed
        values_to_process = y_values[:max_pixels] if max_pixels else y_values

        for l_idx, pixel_idx in enumerate(pixels_to_process):
            _t0 = _time.perf_counter()

            # Extract observation vector c (column of C = identity)
            c_tt = extract_column(C_ttm, pixel_idx)

            # Residual: v = y_obs - <c, x>
            v = values_to_process[l_idx] - dot(c_tt, x_tt)

            # Innovation covariance: s = c^T P c = c^T L L^T c
            Lt_c = matvec(L_ttm.T, c_tt)
            L_Lt_c = matvec(L_ttm, Lt_c)
            s = dot(c_tt, L_Lt_c)

            if abs(s) < 1e-15:
                continue

            # Kalman gain: K = L L^T c * (v / s)
            K_tt = L_Lt_c * (v / s)

            # State update: x = x + K
            x_tt = add(x_tt, K_tt)
            x_tt = orthogonalize(x_tt, 0)

            # ---- Covariance measurement update ----
            # QR of [c^T L; L]^T
            L_t = L_ttm.T

            # Extract columns of L^T + prepend L^T c
            mu_vectors = [Lt_c]  # First vector: L^T c
            for i in range(n_cols):
                mu_vectors.append(extract_column(L_t, i))

            # Gram-Schmidt on the first n_cols vectors
            q_mu = gram_schmidt_tt(mu_vectors[:n_cols], eps=eps, max_rank=max_rank_mu)

            # Combine into Q matrix
            Q_mu_mat = combine_tt_vectors(q_mu, col_partition, eps=eps, max_rank=max_rank_mu)

            # R2 = Q^T @ L^T directly — no need to build augmented L
            R2 = Q_mu_mat.T @ L_t
            R2_tt = ttm_to_tt(R2)
            R2_tt = tt_round(R2_tt, eps=eps)
            R2 = tt_to_ttm(R2_tt, R2.row_shape, R2.col_shape)
            R2_t = R2.T

            # Remove first column, append zero column
            zero_tt = tt.TensorTrain.zeros(row_partition)
            r2t_vectors = []
            for i in range(1, n_cols):
                r2t_vectors.append(extract_column(R2_t, i))
            r2t_vectors.append(zero_tt)

            L_ttm = combine_tt_vectors(r2t_vectors, col_partition, eps=eps, max_rank=max_rank_mu)

            print(f"    [MU] pixel {l_idx+1}/{len(pixels_to_process)}: "
                  f"{_time.perf_counter()-_t0:.1f}s, "
                  f"L ranks={L_ttm.ranks}", flush=True)

        # Reconstruct state estimate for this frame
        x_full = x_tt.full().ravel()
        estimates.append(x_full.copy())

        # Compute RMSE vs original
        if original_video is not None:
            orig_frame = original_video[:, :, k + 1].ravel()
            rmse = np.sqrt(np.mean((x_full - orig_frame) ** 2))
            n_obs = len(observed)
            print(f"  Frame {k + 1}/{T - 1}: {n_obs} observed pixels, RMSE = {rmse:.2f}")
        else:
            print(f"  Frame {k + 1}/{T - 1}: {len(observed)} observed pixels")

    # Stack into video
    result = np.array(estimates).T.reshape(h, w, T)
    return result


# ======================================================================
# Visualization
# ======================================================================


def visualize_results(
    original: np.ndarray,
    masked: np.ndarray,
    estimated: np.ndarray,
    frame_indices: list[int] | None = None,
    save_path: str | None = None,
):
    """Plot original, masked, and estimated frames side by side.

    Parameters
    ----------
    original, masked, estimated : ndarray
        Shape ``(h, w, T)``.
    frame_indices : list of int, optional
        Which frames to show. Defaults to 5 evenly spaced.
    save_path : str, optional
        If given, save figure to this path instead of showing.
    """
    import matplotlib.pyplot as plt

    T = original.shape[2]
    if frame_indices is None:
        frame_indices = np.linspace(1, T - 1, min(5, T - 1), dtype=int).tolist()

    n_frames = len(frame_indices)
    fig, axes = plt.subplots(3, n_frames, figsize=(3 * n_frames, 9))

    if n_frames == 1:
        axes = axes[:, np.newaxis]

    vmin, vmax = original.min(), original.max()

    for col, t in enumerate(frame_indices):
        axes[0, col].imshow(original[:, :, t], cmap="gray", vmin=vmin, vmax=vmax)
        axes[0, col].set_title(f"Original t={t}")
        axes[0, col].axis("off")

        axes[1, col].imshow(masked[:, :, t], cmap="gray", vmin=vmin, vmax=vmax)
        axes[1, col].set_title(f"Masked t={t}")
        axes[1, col].axis("off")

        if t < estimated.shape[2]:
            axes[2, col].imshow(estimated[:, :, t], cmap="gray", vmin=vmin, vmax=vmax)
            axes[2, col].set_title(f"Estimated t={t}")
        axes[2, col].axis("off")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Figure saved to {save_path}")
    else:
        plt.show()


# ======================================================================
# Main
# ======================================================================


def main():
    """Run the full TNKF video inpainting pipeline."""
    # --- Configuration (matches thesis) ---
    video_path = Path(__file__).parent.parent / "matlab-og" / "original_video_TC.mp4"
    resolution = (9, 16)
    row_partition = (3, 3, 4, 2, 2)
    col_partition = (3, 3, 4, 2, 2)
    alpha = 0.95
    beta = 0.01
    eps = 1e-5
    max_frames = None  # None = all frames
    max_pixels_per_frame = 3  # Limit pixel updates per frame (None = all)

    print("=" * 60)
    print("Tensor-Networked Square-Root Kalman Filter (TNKF)")
    print("=" * 60)
    print(f"Video: {video_path}")
    print(f"Resolution: {resolution}")
    print(f"Missing pixels: {alpha * 100:.0f}%")
    print(f"Max frames: {max_frames}")
    print(f"Max pixels/frame: {max_pixels_per_frame}")
    print()

    # Load and preprocess
    print("Loading video...")
    video = load_video_grayscale(video_path, resolution)
    print(f"Video shape: {video.shape}")

    print("Masking video...")
    masked, mask = mask_video(video, alpha=alpha, beta=beta)
    n_obs_avg = mask[:, :, 1:].mean()
    print(f"Average observation rate: {n_obs_avg * 100:.1f}%")
    print()

    # Run TNKF
    estimated = tnkf(
        original_video=video,
        masked_video=masked,
        observation_mask=mask,
        row_partition=row_partition,
        col_partition=col_partition,
        eps=eps,
        max_rank_tp=5,
        max_rank_mu=5,
        max_frames=max_frames,
        max_pixels_per_frame=max_pixels_per_frame,
    )

    # Overall RMSE
    T = estimated.shape[2]
    rmse_total = np.sqrt(np.mean((estimated[:, :, :T] - video[:, :, :T]) ** 2))
    print(f"\nOverall RMSE: {rmse_total:.2f}")

    # Visualize
    print("\nGenerating visualization...")
    visualize_results(
        video[:, :, :T],
        masked[:, :, :T],
        estimated,
        save_path="tnkf_results.png",
    )

    print("\nDone!")


if __name__ == "__main__":
    main()
