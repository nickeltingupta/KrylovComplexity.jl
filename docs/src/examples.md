# Examples

A runnable version of the examples below is in [`examples/krylov_complexity_example.jl`](https://github.com/nickeltingupta/KrylovComplexity.jl/blob/main/examples/krylov_complexity_example.jl).

---

## Example 1: Krylov Complexity for a Single Disorder Realization

```julia
using KrylovComplexity, Combinatorics, Plots, LaTeXStrings

N = 8
validindices = [(i,j,k,l) for (i,j,k,l) in combinations(1:2N, 4) if i < j < k < l]

# Build Hamiltonian and block-diagonalize
H      = build_SpinXY4_bitwise_efficient(42, N, validindices)
blocks = build_parity_blocks(H, N)

# Build infinite-temperature TFD state and run Lanczos
ψ0     = build_tfd_state(blocks.E_plus, blocks.V_plus, 0.0)
result = extract_lanczos_coefficients_fro_blas(blocks.H_plus, ψ0, 128)

# Compute and plot
times = collect(range(0.0, 30.0, length=600))
C_K   = compute_krylov_complexity(result.α, result.β, times)

plot(times, C_K,
     xlabel=L"t", ylabel=L"C_K(t)",
     title="Krylov Complexity — N=$N",
     lw=2, legend=false)
```

---

## Example 2: Temperature Dependence

```julia
using KrylovComplexity, Combinatorics, Plots, LaTeXStrings

N = 8
validindices = [(i,j,k,l) for (i,j,k,l) in combinations(1:2N, 4) if i < j < k < l]
H      = build_SpinXY4_bitwise_efficient(1, N, validindices)
blocks = build_parity_blocks(H, N)

βs    = [0.0, 0.5, 1.0, 2.0, 5.0]
times = collect(range(0.0, 40.0, length=800))

p = plot(xlabel=L"t", ylabel=L"C_K(t)",
         title="Temperature Dependence — N=$N", framestyle=:box)

for β in βs
    ψ0     = build_tfd_state(blocks.E_plus, blocks.V_plus, β)
    result = extract_lanczos_coefficients_fro_blas(blocks.H_plus, ψ0, 128)
    C_K    = compute_krylov_complexity(result.α, result.β, times)
    plot!(p, times, C_K, lw=2, label=L"\beta = %$β")
end

display(p)
```

---

## Example 3: Disorder Averaging

```julia
using KrylovComplexity, Combinatorics, Statistics, Plots, LaTeXStrings

N = 8
n_disorders = 20
validindices = [(i,j,k,l) for (i,j,k,l) in combinations(1:2N, 4) if i < j < k < l]
times = collect(range(0.0, 30.0, length=500))

C_K_all = zeros(length(times), n_disorders)

for seed in 1:n_disorders
    H      = build_SpinXY4_bitwise_efficient(seed, N, validindices)
    blocks = build_parity_blocks(H, N)
    ψ0     = build_tfd_state(blocks.E_plus, blocks.V_plus, 0.0)
    result = extract_lanczos_coefficients_fro_blas(blocks.H_plus, ψ0, 128)
    C_K_all[:, seed] = compute_krylov_complexity(result.α, result.β, times)
end

μ = mean(C_K_all, dims=2) |> vec
σ = std(C_K_all,  dims=2) |> vec

plot(times, μ, ribbon=σ, fillalpha=0.3, lw=2, label="Mean ± std",
     xlabel=L"t", ylabel=L"\langle C_K(t) \rangle",
     title="Disorder-averaged Krylov Complexity — N=$N, $(n_disorders) realizations")
```

---

## Example 4: Lanczos Coefficient Spectrum

The Lanczos coefficients ``\{\alpha_n, \beta_n\}`` themselves carry information about quantum chaos. For chaotic systems, the ``\beta_n`` coefficients grow linearly before saturating — a behavior known as "Lanczos coefficient universality."

```julia
using KrylovComplexity, Combinatorics, Plots, LaTeXStrings

N = 10
validindices = [(i,j,k,l) for (i,j,k,l) in combinations(1:2N, 4) if i < j < k < l]
H      = build_SpinXY4_bitwise_efficient(1, N, validindices)
blocks = build_parity_blocks(H, N)
ψ0     = build_tfd_state(blocks.E_plus, blocks.V_plus, 0.0)
result = extract_lanczos_coefficients_fro_blas(blocks.H_plus, ψ0, size(blocks.H_plus,1))

ns = 1:length(result.α)

p1 = scatter(ns, result.α, ms=2, label=L"\alpha_n",
             xlabel="n", ylabel="Lanczos coefficient", title="α coefficients")
p2 = scatter(1:length(result.β), result.β, ms=2, color=:red, label=L"\beta_n",
             xlabel="n", title="β coefficients")

plot(p1, p2, layout=(1,2), size=(900,350))
```

---

## Example 5: Level Spacing Statistics

```julia
using KrylovComplexity, Combinatorics, Plots

N = 10
validindices = [(i,j,k,l) for (i,j,k,l) in combinations(1:2N, 4) if i < j < k < l]
H      = build_SpinXY4_bitwise_efficient(1, N, validindices)
blocks = build_parity_blocks(H, N)

# Analyze each parity sector separately
p_plus  = level_spacing_plot(blocks.E_plus;  num_bins=30, verbose=true)
p_minus = level_spacing_plot(blocks.E_minus; num_bins=30)

title!(p_plus,  "Level Spacing — +1 sector, N=$N")
title!(p_minus, "Level Spacing — −1 sector, N=$N")

plot(p_plus, p_minus, layout=(1,2), size=(1000,400))
```
