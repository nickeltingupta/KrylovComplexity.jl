# KrylovComplexity.jl

[![Julia](https://img.shields.io/badge/Julia-1.9%2B-blue?logo=julia)](https://julialang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

A Julia package for computing **Krylov complexity** (also called spread complexity) in disordered quantum spin systems. Starting from a Hamiltonian and an initial state, the package runs the Lanczos algorithm, builds the Krylov basis, and evaluates the time-dependent complexity

$$C_K(t) = \sum_n n\,|c_n(t)|^2.$$

Included out of the box:
- A disordered **SpinXY4 Hamiltonian** (4-body σˣ/σʸ interactions with Gaussian random couplings, related to the SYK family of models)
- **Parity-sector block diagonalization**
- **Nearest-neighbor level-spacing statistics** with Wigner–Dyson GUE/GOE overlay
- **Thermofield Double (TFD) state** construction at arbitrary inverse temperature β
- **Lanczos algorithm** with full re-orthogonalization (BLAS-accelerated)
- **Krylov / Spread complexity** time evolution via tridiagonal eigendecomposition

---

## Table of Contents

- [Background](#background)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [API Reference](#api-reference)
- [Examples](#examples)
- [Project Structure](#project-structure)
- [Contributing](#contributing)
- [License](#license)
- [References](#references)

---

## Background

Krylov complexity quantifies quantum operator or state growth under time evolution. Given an initial state |ψ₀⟩ and a Hamiltonian H, the Lanczos algorithm constructs an orthonormal **Krylov basis** {|Kₙ⟩} such that H is tridiagonal with diagonal entries αₙ and off-diagonal entries βₙ. The time-evolved state decomposes as

$$|\psi(t)\rangle = \sum_n c_n(t)\,|K_n\rangle,$$

and the Krylov complexity measures how far the state has spread along the chain. In chaotic systems governed by random matrix statistics the Lanczos coefficients βₙ grow linearly, driving linear growth in C_K(t), while integrable systems show slower, bounded growth.

---

## Installation

This package is not yet registered in the Julia General Registry. Install directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/nickeltingupta/KrylovComplexity.jl")
```

Or in the REPL package manager (press `]`):

```
pkg> add https://github.com/nickeltingupta/KrylovComplexity.jl
```

**Julia version:** 1.9 or later is recommended.

---

## Quick Start

```julia
using KrylovComplexity, SparseArrays
using Combinatorics   # for combinations()

# 1. System parameters
N = 8   # number of spins  →  Hilbert-space dim = 2^N = 256

# 2. Build the list of valid 4-body interaction sites
validindices = [(i,j,k,l)
    for (i,j,k,l) in combinations(1:2N, 4)
    if i < j < k < l]

# 3. Build the Hamiltonian (disorder seed = 42)
H = build_SpinXY4_bitwise_efficient(42, N, validindices)

# 4. Block-diagonalize into parity sectors
blocks = build_parity_blocks(H, N)
# blocks.E_plus, blocks.V_plus   → even-parity energies & eigenvectors
# blocks.E_minus, blocks.V_minus → odd-parity energies & eigenvectors

# 5. Check level statistics (chaotic ↔ Wigner–Dyson GUE)
p = level_spacing_plot(vcat(blocks.E_plus, blocks.E_minus))
display(p)

# 6. Build the initial state (TFD at β=0 → infinite temperature)
ψ0 = build_tfd_state(blocks.E_plus, blocks.V_plus, 0.0)

# 7. Run the Lanczos algorithm to get Krylov-chain coefficients
KDim  = size(blocks.H_plus, 1)
coeff = extract_lanczos_coefficients_fro_blas(
            sparse(blocks.H_plus), ψ0, KDim)

# 8. Compute Krylov complexity C_K(t) over time
times      = collect(range(0.0, 20.0, length=500))
complexity = compute_krylov_complexity(coeff.α, coeff.β, times)

# 9. Plot
using Plots
plot(times, complexity,
     xlabel="t", ylabel="C_K(t)",
     title="Krylov Complexity  (N=$N, β=0)", lw=2, legend=false)
```

---

## API Reference

### Hamiltonian Construction

#### `build_SpinXY4_bitwise_efficient(disorder_seed, N, validindices)`

Builds the disordered SpinXY4 Hamiltonian using pre-allocated vectors and a single `sparse(I, J, V)` call — the fastest available method.

```julia
H = build_SpinXY4_bitwise_efficient(
        disorder_seed :: Int,
        N             :: Int,
        validindices  :: Vector{NTuple{4,Int64}}
    ) → SparseMatrixCSC{ComplexF64}
```

| Argument | Type | Description |
|---|---|---|
| `disorder_seed` | `Int` | RNG seed for the Gaussian couplings J ~ N(0,1) |
| `N` | `Int` | Number of spins; Hilbert-space dimension is 2ᴺ |
| `validindices` | `Vector{NTuple{4,Int64}}` | All ordered 4-tuples (i<j<k<l) from {1,…,2N} |

The Hamiltonian is

$$H = \sqrt{\frac{3}{4N^3}} \sum_{i<j<k<l} J_{ijkl}\, i^{H_F}\, \mathcal{O}_i \mathcal{O}_j \mathcal{O}_k \mathcal{O}_l$$

where odd-labelled sites carry σˣ and even-labelled sites carry σʸ. The phase factor i^{H_F} (H_F counts adjacent XY pairs) ensures Hermiticity.

#### `build_SpinXY4_bitwise_old(disorder_seed, N, validindices)`

Functionally identical to the efficient version but builds the sparse matrix by accumulating entries in a loop. Useful for debugging.

---

### Block Diagonalization

#### `build_parity_blocks(H, N)`

Decomposes H into even (Γ=+1) and odd (Γ=−1) parity sectors, where the parity operator Γ = (−1)^n̂ counts the number of up-spins. Returns a named tuple:

```julia
blocks = build_parity_blocks(
    H :: AbstractMatrix{ComplexF64},
    N :: Int
)
# → (H_plus, H_minus, E_plus, V_plus, E_minus, V_minus)
```

| Field | Size | Description |
|---|---|---|
| `H_plus` / `H_minus` | 2^(N-1) × 2^(N-1) | Dense sub-block matrices |
| `E_plus` / `E_minus` | 2^(N-1) | Sorted eigenvalues |
| `V_plus` / `V_minus` | 2^(N-1) × 2^(N-1) | Eigenvector matrices |

Working inside a single parity sector halves the Hilbert-space dimension, making both the Lanczos run and exact diagonalization significantly cheaper.

---

### Level Spacing Statistics

#### `level_spacing_plot(eigenvalues; num_bins=10, verbose=false)`

Unfolds the spectrum, computes nearest-neighbor level spacings, evaluates the level-spacing ratio r, and overlays the Wigner–Dyson GUE and GOE distributions.

```julia
p = level_spacing_plot(
        eigenvalues :: AbstractVector;
        num_bins    :: Int  = 10,
        verbose     :: Bool = false
    ) → Plots.Plot
```

The computed r value appears in the histogram legend. Reference values: GUE → r ≈ 0.6027, GOE → r ≈ 0.5359, Poisson (integrable) → r ≈ 0.3863.

---

### Thermofield Double State

#### `build_tfd_state(eigenvals, eigenvecs, β)`

Constructs the Thermofield Double state

$$|\text{TFD}(\beta)\rangle = \frac{1}{\sqrt{Z(\beta)}} \sum_n e^{-\beta E_n/2}\,|n\rangle, \qquad Z(\beta)=\sum_n e^{-\beta E_n}$$

At β=0 this reduces to the equal-weight infinite-temperature state.

```julia
ψ = build_tfd_state(
        eigenvals :: Vector{Float64},
        eigenvecs :: Matrix{ComplexF64},
        β         :: Float64 = 0.0
    ) → Vector{ComplexF64}
```

---

### Lanczos Algorithm

#### `extract_lanczos_coefficients_fro_blas(H, ψ0, KDim; verbose=false)`

Implements the **Lanczos algorithm with Full Re-Orthogonalization (FRO)** using two-pass Gram–Schmidt and BLAS matrix-vector products. Returns the tridiagonal Krylov-chain coefficients.

```julia
coeff = extract_lanczos_coefficients_fro_blas(
            H     :: SparseMatrixCSC{ComplexF64},
            ψ0    :: Vector{ComplexF64},
            KDim  :: Int;
            verbose :: Bool = false
        )
# → (α = Vector{Float64}, β = Vector{Float64})
```

| Return field | Length | Description |
|---|---|---|
| `coeff.α` | ≤ KDim | Diagonal Lanczos coefficients |
| `coeff.β` | ≤ KDim−1 | Off-diagonal Lanczos coefficients |

Early termination occurs when βⱼ ≤ √ε_machine (Krylov chain has been exhausted).

---

### Krylov Complexity

#### `compute_krylov_complexity(αn, βn, times; max_dim, show_progress)`

Diagonalizes the K×K tridiagonal Krylov matrix (O(K³)), precomputes the amplitude matrix (O(K²)), then evaluates C_K(t) for every requested time point via a matrix–vector product (O(K) per time point).

```julia
C = compute_krylov_complexity(
        αn           :: Vector{Float64},
        βn           :: Vector{Float64},
        times        :: Vector{Float64};
        max_dim      :: Int  = length(αn),
        show_progress :: Bool = false
    ) → Vector{Float64}
```

| Argument | Description |
|---|---|
| `αn`, `βn` | Lanczos coefficients from `extract_lanczos_coefficients_fro_blas` |
| `times` | Time points at which to evaluate C_K(t) |
| `max_dim` | Truncate Krylov dimension (useful to test convergence) |
| `show_progress` | Show a ProgressBar during the loop |

---

## Examples

### Disorder-averaged Krylov complexity

```julia
using KrylovComplexity, LinearAlgebra, SparseArrays, Combinatorics, Statistics, Plots

N            = 10
n_disorders  = 50
validindices = [(i,j,k,l) for (i,j,k,l) in combinations(1:2N, 4) if i < j < k < l]
times        = collect(range(0.0, 30.0, length=800))
β            = 1.0

all_C = Matrix{Float64}(undef, length(times), n_disorders)
for seed in 1:n_disorders
    H      = build_SpinXY4_bitwise_efficient(seed, N, validindices)
    blocks = build_parity_blocks(H, N)
    ψ0     = build_tfd_state(blocks.E_plus, blocks.V_plus, β)
    coeff  = extract_lanczos_coefficients_fro_blas(
                 sparse(blocks.H_plus), ψ0, size(blocks.H_plus, 1))
    all_C[:, seed] = compute_krylov_complexity(coeff.α, coeff.β, times)
end

plot(times, mean(all_C; dims=2)[:],
     ribbon = std(all_C; dims=2)[:],
     xlabel="t", ylabel="⟨C_K(t)⟩",
     title="Disorder-averaged Krylov Complexity  (N=$N, β=$β, n=$n_disorders)",
     lw=2, fillalpha=0.2, legend=false)
```

### Scanning inverse temperature β

```julia
using KrylovComplexity, SparseArrays, Combinatorics, Plots

N            = 8
validindices = [(i,j,k,l) for (i,j,k,l) in combinations(1:2N, 4) if i < j < k < l]
H            = build_SpinXY4_bitwise_efficient(1, N, validindices)
blocks       = build_parity_blocks(H, N)
times        = collect(range(0.0, 25.0, length=600))

pl = plot(xlabel="t", ylabel="C_K(t)",
          title="Krylov Complexity vs β  (N=$N)")
for β in [0.0, 0.5, 1.0, 2.0, 5.0]
    ψ0    = build_tfd_state(blocks.E_plus, blocks.V_plus, β)
    coeff = extract_lanczos_coefficients_fro_blas(
                sparse(blocks.H_plus), ψ0, size(blocks.H_plus, 1))
    C     = compute_krylov_complexity(coeff.α, coeff.β, times)
    plot!(pl, times, C, label="β = $β", lw=2)
end
display(pl)
```

### Lanczos coefficient profile (diagnostic)

```julia
using KrylovComplexity, SparseArrays, Combinatorics, Plots

N            = 10
validindices = [(i,j,k,l) for (i,j,k,l) in combinations(1:2N, 4) if i < j < k < l]
H            = build_SpinXY4_bitwise_efficient(1, N, validindices)
blocks       = build_parity_blocks(H, N)
ψ0           = build_tfd_state(blocks.E_plus, blocks.V_plus, 0.0)
coeff        = extract_lanczos_coefficients_fro_blas(
                   sparse(blocks.H_plus), ψ0, size(blocks.H_plus, 1))

p1 = plot(coeff.α, xlabel="n", ylabel="αₙ", title="Diagonal Lanczos coefficients", lw=2)
p2 = plot(coeff.β, xlabel="n", ylabel="βₙ", title="Off-diagonal Lanczos coefficients", lw=2)
plot(p1, p2, layout=(1,2), size=(900,350))
```

---

## Project Structure

```
KrylovComplexity.jl/
├── src/
│   ├── KrylovComplexity.jl     # Module entry point & exports
│   ├── hamiltonian.jl          # SpinXY4 Hamiltonian builders
│   ├── block_diagonalize.jl    # Parity-sector decomposition
│   ├── level_spacing.jl        # Level spacing statistics & plots
│   ├── lanczos.jl              # Lanczos-FRO algorithm (BLAS)
│   ├── TFD.jl                  # Thermofield Double state
│   └── spreadcomplexity.jl     # Krylov / spread complexity
├── test/
│   └── runtests.jl
├── examples/
│   ├── basic_usage.jl
│   └── disorder_average.jl
├── Project.toml
├── Manifest.toml
├── CHANGELOG.md
├── CONTRIBUTING.md
├── LICENSE
└── README.md
```

---

## Contributing

Contributions, bug reports, and feature requests are welcome. Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on how to get involved.

---

## License

This project is licensed under the **MIT License** — see [LICENSE](LICENSE) for details.

---

## References

1. Parker, D. E., Cao, X., Avery, A., Scaffidi, T., & Altman, E. (2019). *A Universal Operator Growth Hypothesis*. [arXiv:1812.02666](https://arxiv.org/abs/1812.02666)
2. Balasubramanian, V., Caputa, P., Magan, J. M., & Wu, Q. (2022). *Quantum chaos and the complexity of spread complexity*. [arXiv:2202.06957](https://arxiv.org/abs/2202.06957)
3. Erdmenger, J., Grosvenor, K. T., & Volber, R. (2023). *Spread complexity in free fermion models*. [arXiv:2302.11384](https://arxiv.org/abs/2302.11384)
4. Sachdev, S., & Ye, J. (1993). *Gapless spin-fluid ground state in a random quantum Heisenberg magnet*. Phys. Rev. Lett. **70**, 3339.
5. Kitaev, A. (2015). *A simple model of quantum holography*. KITP talks. http://online.kitp.ucsb.edu/online/entangled15/kitaev/
