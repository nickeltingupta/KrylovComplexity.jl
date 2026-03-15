using KrylovComplexity
using Test, LinearAlgebra, SparseArrays, Combinatorics

@testset "KrylovComplexity.jl" begin

    # ── Shared setup ──────────────────────────────────────────────────────────
    N            = 8
    validindices = [(i,j,k,l)
                    for (i,j,k,l) in combinations(1:2N, 4)
                    if i < j < k < l]
    dim          = 1 << N   # 64

    # ── Hamiltonian ───────────────────────────────────────────────────────────
    @testset "Hamiltonian" begin
        H = build_SpinXY4_bitwise_efficient(1, N, validindices)

        @test size(H) == (dim, dim)
        @test issparse(H)

        # Hermiticity  (H should equal H†)
        @test norm(H - H') < 1e-10

        # Reproducibility — same seed → same matrix
        H2 = build_SpinXY4_bitwise_efficient(1, N, validindices)
        @test H == H2

        # Different seeds → different matrices
        H3 = build_SpinXY4_bitwise_efficient(2, N, validindices)
        @test H != H3

        # Consistency between implementations
        H_old = build_SpinXY4_bitwise_old(1, N, validindices)
        @test norm(H - H_old) < 1e-10
    end

    # ── Block diagonalization ─────────────────────────────────────────────────
    @testset "Parity blocks" begin
        H      = build_SpinXY4_bitwise_efficient(1, N, validindices)
        blocks = build_parity_blocks(H, N)

        half = dim ÷ 2
        @test size(blocks.H_plus)  == (half, half)
        @test size(blocks.H_minus) == (half, half)
        @test length(blocks.E_plus)  == half
        @test length(blocks.E_minus) == half

        # Eigenvalues should be real-valued sorted arrays
        @test issorted(blocks.E_plus)
        @test issorted(blocks.E_minus)

        # Eigenvectors should be unitary
        @test norm(blocks.V_plus' * blocks.V_plus - I(half))  < 1e-10
        @test norm(blocks.V_minus' * blocks.V_minus - I(half)) < 1e-10
    end

    # ── TFD state ─────────────────────────────────────────────────────────────
    @testset "TFD state" begin
        H      = build_SpinXY4_bitwise_efficient(1, N, validindices)
        blocks = build_parity_blocks(H, N)
        half   = dim ÷ 2

        # β = 0 → equal-weight superposition → all |coeffs| should be equal
        ψ_inf = build_tfd_state(blocks.E_plus, blocks.V_plus, 0.0)
        @test length(ψ_inf) == half
        @test abs(norm(ψ_inf) - 1.0) < 1e-12

        coeffs_in_eigenbasis = abs.(blocks.V_plus' * ψ_inf)
        @test maximum(coeffs_in_eigenbasis) - minimum(coeffs_in_eigenbasis) < 1e-12

        # Finite β → normalized
        ψ_thermal = build_tfd_state(blocks.E_plus, blocks.V_plus, 2.0)
        @test abs(norm(ψ_thermal) - 1.0) < 1e-12

        # Higher β → state dominated by ground state
        ψ_cold = build_tfd_state(blocks.E_plus, blocks.V_plus, 100.0)
        @test abs(norm(ψ_cold) - 1.0) < 1e-12
        @test abs2(ψ_cold' * blocks.V_plus[:, 1]) > 0.999
    end

    # ── Lanczos algorithm ─────────────────────────────────────────────────────
    @testset "Lanczos-FRO" begin
        H      = build_SpinXY4_bitwise_efficient(1, N, validindices)
        blocks = build_parity_blocks(H, N)
        half   = dim ÷ 2

        ψ0    = build_tfd_state(blocks.E_plus, blocks.V_plus, 0.0)
        coeff = extract_lanczos_coefficients_fro_blas(
                    sparse(blocks.H_plus), ψ0, half)

        @test length(coeff.α) == half
        @test length(coeff.β) == half - 1

        # α coefficients must be real (they are Float64 already, just sanity-check)
        @test all(isfinite, coeff.α)
        @test all(isfinite, coeff.β)

        # β coefficients must be non-negative
        @test all(>=(0), coeff.β)
    end

    # ── Krylov complexity ─────────────────────────────────────────────────────
    @testset "Krylov complexity" begin
        H      = build_SpinXY4_bitwise_efficient(1, N, validindices)
        blocks = build_parity_blocks(H, N)
        half   = dim ÷ 2

        ψ0    = build_tfd_state(blocks.E_plus, blocks.V_plus, 0.0)
        coeff = extract_lanczos_coefficients_fro_blas(
                    sparse(blocks.H_plus), ψ0, half)

        times = [0.0, 1.0, 5.0, 10.0]
        C     = compute_krylov_complexity(coeff.α, coeff.β, times)

        @test length(C) == length(times)
        @test all(isfinite, C)

        # C_K(0) = 0  (state starts at |K_0⟩)
        @test abs(C[1]) < 1e-10

        # Complexity is non-negative at all times
        @test all(>=(0), C)

        # Complexity is bounded by K-1 (maximum basis index)
        @test all(<=(half - 1 + 1e-10), C)

        # max_dim truncation is respected
        C_trunc = compute_krylov_complexity(coeff.α, coeff.β, times; max_dim=10)
        @test all(<=(9 + 1e-10), C_trunc)
    end

end
