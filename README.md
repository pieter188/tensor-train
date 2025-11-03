# Tensor Train Decomposition

A Python implementation of Tensor Train (TT) decomposition methods, evolved from a Master's thesis on Tensor-Networked Kalman Filters.

## Overview

This repository contains implementations of Tensor Train decomposition and related algorithms. The initial implementation was done in MATLAB as part of a Master's thesis, but is now being rebuilt in Python with focus on performance, efficiency, and modern best practices.

## About Tensor Networks

As stated in the introduction of the original thesis, this work covers the fundamental methods of tensor networks, which are important to be able to understand and improve the existing technique: Tensor-Networked Kalman Filter (TNKF). Therefore, the focus is on Tensor Networks (TN).

TN are efficient in data storage and are intuitive in their application for matrix-matrix and matrix-vector computations. The implementation covers:

- **Basic tensor definitions and operations** - Fundamental tensor operations including inner products, matricization, and n-mode products
- **Tensor decompositions** - Primarily the Tensor-Train (TT) decomposition using TT-SVD algorithm
- **TT and Tensor-Train matrix (TTm) operations** - Efficient operations on tensor train representations
- **Core algorithms** - TT-rounding, contraction, and various tensor network functions

### Key Concepts

A **tensor** can be described as a type of array, where:
- A vector has 1 index
- A matrix has 2 indices
- A tensor can have d ≥ 3 indices (also known as a d-way or dth order tensor)

Tensors are visualized using **tensor diagrams** (as presented by Penrose), where each branch represents a single index of the tensor.

### Tensor-Train Decomposition

A d-dimensional tensor **X** ∈ ℝ^(I₁×I₂×...×Iₐ) written in TT-format can be decomposed into d cores, where each core is a 3-way tensor:

```
X(i₁, i₂, ..., iₐ) = X⁽¹⁾(i₁)X⁽²⁾(i₂)···X⁽ᵈ⁾(iₐ)
```

Each core has dimension X⁽ᵏ⁾ ∈ ℝ^(Rₖ₋₁×Iₖ×Rₖ), where Rₖ are the TT-ranks and R₀ = Rₐ = 1.

### Storage Complexity Comparison

| Decomposition | Storage Complexity |
|---------------|-------------------|
| Full Tensor   | O(Iᵈ) |
| TT            | O(dIR²) |
| TTm           | O(dI²R²) |
| CPD           | O(dIR) |
| Tucker        | O(dIR + Rᴺ) |

where I = max{I₁, I₂, ..., Iₐ}, R = max{R₁, R₂, ..., Rₐ₋₁}

## Repository Structure

```
tensor-train.git/
├── matlab-og/           # Original MATLAB implementation from thesis
├── papers/              # Research papers and references
├── python/              # Modern Python implementation (work in progress)
├── starter_exercises/   # Educational exercises and examples
└── README.md
```

## Features (Python Implementation Goals)

- **TT-SVD Algorithm** - Singular Value Decomposition-based tensor train construction
- **TT Operations**:
  - Addition and concatenation
  - Scalar multiplication and transpose
  - Matrix-vector and matrix-matrix products
  - Inner products and contractions
- **TT-Rounding** - Rank reduction and approximation control
- **Performance Optimization** - Leveraging modern Python libraries (NumPy, SciPy, etc.)
- **Efficient Storage** - Significant memory reduction compared to full tensor storage

## Mathematical Operations

### Core Operations Implemented

1. **Matricization** - Unfolding tensors into matrix representations
2. **n-mode Product** - Tensor-matrix multiplication over specific indices
3. **TT-SVD** - Decomposition algorithm with controlled approximation error
4. **TT-Rounding** - QR-factorization based rank truncation
5. **Contraction** - Optimal combination of connected tensor trains

### Time Complexities

- **Dot Product**: O((d-2)(I(RS)²))
- **Matrix-Vector Product**: O(dIJ(RS)²)
- **Matrix-Matrix Product**: O(dIJL(RS)²)
- **TT-Rounding**: O(dIR³)

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/tensor-train.git
cd tensor-train

# Install dependencies (once Python implementation is ready)
pip install -r python/requirements.txt
```

## Usage Example

```python
# Example usage (to be implemented)
import tensor_train as tt

# Create a tensor
X = ...  # Your high-dimensional tensor

# Decompose into TT format
tt_cores = tt.tt_svd(X, epsilon=1e-6)

# Perform operations in TT format
Y = tt.matvec(A_tt, x_tt)  # Matrix-vector product
Z = tt.matmat(A_tt, B_tt)  # Matrix-matrix product

# Round to reduce ranks
Y_rounded = tt.round(Y, epsilon=1e-6)
```

## Background

This work originated from a Master's thesis focusing on Tensor-Networked Kalman Filters (TNKF). The thesis explored efficient methods for handling high-dimensional data through tensor network decompositions, with particular emphasis on the Tensor-Train format for its simplicity and efficiency in:

- Low-rank matrix approximations
- Simple and intuitive tensor operations
- Efficient data storage
- Matrix and vector computations

## References

Key papers and resources can be found in the `papers/` directory. The implementation is based on foundational work in tensor networks, including:

- TT-SVD algorithms by Oseledets
- Tensor decomposition methods
- Efficient tensor operations and contractions

## Development Status

🚧 **Work in Progress** 🚧

The Python implementation is currently under active development. The MATLAB implementation serves as a reference for the improved Python version.

**Goals:**
- ✅ Document mathematical foundations
- 🔄 Implement core TT-SVD algorithm in Python
- 🔄 Implement TT operations (addition, multiplication, contraction)
- 🔄 Optimize performance using NumPy/SciPy
- 📋 Add comprehensive test suite
- 📋 Create example notebooks
- 📋 Benchmark against MATLAB implementation

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## License

[Add your license here]

## Citation

If you use this code in your research, please cite:

```bibtex
@mastersthesis{vanklaveren_tensor_networks,
  author = {P. van Klaveren},
  title = {Tensor Networks and Tensor-Networked Kalman Filters},
  school = {[Your University]},
  year = {[Year]},
}
```

## Contact

[Your contact information]

---

*This implementation aims to make tensor train decomposition accessible, efficient, and easy to use for researchers and practitioners working with high-dimensional data.*
