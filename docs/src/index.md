# KrylovComplexity.jl

**KrylovComplexity.jl** is a Julia package for computing Krylov complexity (spread complexity) in quantum many-body systems, built around the random SpinXY4 model.

## What is Krylov Complexity?

Given a Hamiltonian ``H`` and an initial state ``|\psi_0\rangle``, the Lanczos algorithm generates an orthonormal basis for the Krylov subspace:

```math
\mathcal{K}_K = \text{span}\{|\psi_0\rangle,\, H|\psi_0\rangle,\, H^2|\psi_0\rangle,\, \ldots\}
```

The time-evolved state ``|\psi(t)\rangle = e^{-iHt}|\psi_0\rangle`` spreads through this basis. **Krylov complexity** measures how far it has spread:

```math
C_K(t) = \sum_{n=0}^{K-1} n\, |c_n(t)|^2, \qquad c_n(t) = \langle K_n | \psi(t)\rangle
```

This quantity grows linearly at late times in maximally chaotic systems, and has been proposed as a holographic dual to bulk complexity.

## Package Features

- Efficient sparse construction of the random **SpinXY4 Hamiltonian** (related to SYK)
- **Parity block diagonalization** to reduce matrix sizes by 2×
- **BLAS-accelerated Lanczos** with full reorthogonalization for stable coefficient extraction
- **Thermofield double (TFD)** state construction at arbitrary temperature
- Fast **Krylov complexity** evaluation using eigendecomposition
- **Level spacing statistics** with Wigner–Dyson comparison for quantum chaos diagnostics

## Quick Installation

```julia
using Pkg
Pkg.add(url="https://github.com/nickeltingupta/KrylovComplexity.jl")
```

## Minimal Example

```julia
using KrylovComplexity, Combinatorics

N = 8
validindices = [(i,j,k,l) for (i,j,k,l) in combinations(1:2N, 4) if i < j < k < l]

H      = build_SpinXY4_bitwise_efficient(42, N, validindices)
blocks = build_parity_blocks(H, N)
ψ0     = build_tfd_state(blocks.E_plus, blocks.V_plus, 0.0)
result = extract_lanczos_coefficients_fro_blas(blocks.H_plus, ψ0, 128)

times = collect(0.0:0.05:30.0)
C_K   = compute_krylov_complexity(result.α, result.β, times)
```

## Contents

```@contents
Pages = ["api.md", "physics.md", "examples.md"]
Depth = 2
```
