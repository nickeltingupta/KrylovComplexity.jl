#!/usr/bin/env julia
# examples/krylov_complexity_example.jl
#
# Full worked example: compute Krylov complexity for the SpinXY4 model,
# averaged over multiple disorder realizations.
#
# Run with: julia --project=. examples/krylov_complexity_example.jl

using KrylovComplexity
using Combinatorics, LinearAlgebra, Statistics, Plots, LaTeXStrings

# ============================================================
# Parameters
# ============================================================
N             = 8          # Number of spins (Hilbert space: 2^N = 256)
n_disorders   = 5          # Number of disorder realizations to average
β_temp        = 0.0        # Inverse temperature (0 = infinite temperature TFD)
KDim          = 128        # Krylov dimension (≤ 2^(N-1) = 128 for each parity block)
t_max         = 30.0
n_times       = 400

println("="^60)
println("KrylovComplexity.jl — Example: N=$N, $(n_disorders) disorder realizations")
println("="^60)

# ============================================================
# Precompute valid indices (same for all disorder realizations)
# ============================================================
validindices = [(i,j,k,l) for (i,j,k,l) in combinations(1:2N, 4) if i < j < k < l]
println("Number of interaction terms: $(length(validindices))")
println("Hilbert space dimension: 2^$N = $(1 << N)")
println()

# ============================================================
# Main loop over disorder realizations
# ============================================================
times = collect(range(0.0, t_max, length=n_times))
C_K_all = zeros(n_times, n_disorders)  # store each realization

for seed in 1:n_disorders
    println("─── Disorder realization $seed / $n_disorders ───")

    # 1. Build Hamiltonian
    print("  Building SpinXY4 Hamiltonian... ")
    H = build_SpinXY4_bitwise_efficient(seed, N, validindices)
    println("done  (nnz = $(nnz(H)))")

    # 2. Block-diagonalize by parity
    print("  Block-diagonalizing... ")
    blocks = build_parity_blocks(H, N)
    println("done  (+sector: $(length(blocks.E_plus)) states)")

    # 3. Build TFD initial state in the +1 parity sector
    ψ0 = build_tfd_state(blocks.E_plus, blocks.V_plus, β_temp)
    println("  TFD state norm: $(round(norm(ψ0), sigdigits=6))")

    # 4. Extract Lanczos coefficients
    print("  Running Lanczos (KDim=$KDim)... ")
    result = extract_lanczos_coefficients_fro_blas(blocks.H_plus, ψ0, KDim)
    println("done  (K = $(length(result.α)))")

    # 5. Compute Krylov complexity
    print("  Computing Krylov complexity... ")
    C_K_all[:, seed] = compute_krylov_complexity(result.α, result.β, times)
    println("done")

    println()
end

# ============================================================
# Average over disorders and plot
# ============================================================
C_K_mean = mean(C_K_all; dims=2) |> vec
C_K_std  = std(C_K_all;  dims=2) |> vec

p1 = plot(
    times, C_K_mean,
    ribbon=C_K_std,
    fillalpha=0.3,
    lw=2,
    label="Mean ± std ($(n_disorders) realizations)",
    xlabel=L"t",
    ylabel=L"C_K(t)",
    title="Krylov Complexity — SpinXY4, N=$N, β=$β_temp",
    framestyle=:box,
    legend=:topleft,
)

# Also show individual realizations in gray
for i in 1:n_disorders
    plot!(p1, times, C_K_all[:, i], lw=0.8, alpha=0.4, color=:gray, label=false)
end

savefig(p1, "krylov_complexity_N$(N).pdf")
println("Saved: krylov_complexity_N$(N).pdf")

# ============================================================
# Level spacing statistics (for one realization)
# ============================================================
H_ref = build_SpinXY4_bitwise_efficient(1, N, validindices)
blocks_ref = build_parity_blocks(H_ref, N)
# Use only one parity sector for clean statistics
p2 = level_spacing_plot(blocks_ref.E_plus; num_bins=20, verbose=true)
title!(p2, "Level Spacing — SpinXY4 +1 sector, N=$N")

savefig(p2, "level_spacing_N$(N).pdf")
println("Saved: level_spacing_N$(N).pdf")
