function extract_lanczos_coefficients_fro_blas(H::SparseMatrixCSC{ComplexF64, Int64}, ψ0::Vector{ComplexF64}, KDim::Int; verbose=false)

    dim = length(ψ0)
    bn_pr = sqrt(eps(Float64))
    
    # Storage
    Q = Matrix{ComplexF64}(undef, dim, KDim + 1)
    w = Vector{ComplexF64}(undef, dim)
    
    # Pre-allocations
    αn = Vector{Float64}(undef, KDim)
    βn = Vector{Float64}(undef, KDim-1)
    overlaps = Vector{ComplexF64}(undef, KDim)

    # Initial vector
    Q[:, 1] .= ψ0 ./ norm(ψ0)
    
    for j in 1:KDim
        q_j = view(Q, :, j)

        # Apply Hamiltonian
        mul!(w, H, q_j)
        
        # Compute α
        αn[j] = real(dot(q_j, w))
        
        # w = w - α_j * q_j
        BLAS.axpy!(-αn[j], q_j, w)
        
        # Full reorthogonalization
        if j > 1
            Q_prev = view(Q, :, 1:j-1)
            overlaps_view = view(overlaps, 1:j-1)
            
            # First pass: overlaps = Q_prev' * w
            # Use proper BLAS call for complex matrices
            BLAS.gemv!('C', ComplexF64(1.0), Q_prev, w, ComplexF64(0.0), overlaps_view)
            
            # w = w - Q_prev * overlaps
            BLAS.gemv!('N', ComplexF64(-1.0), Q_prev, overlaps_view, ComplexF64(1.0), w)
            
            # Second pass for stability
            BLAS.gemv!('C', ComplexF64(1.0), Q_prev, w, ComplexF64(0.0), overlaps_view)
            BLAS.gemv!('N', ComplexF64(-1.0), Q_prev, overlaps_view, ComplexF64(1.0), w)
        end
        
        if j < KDim
            βn[j] = BLAS.nrm2(length(w), w, 1)
            
            # Convergence check
            if βn[j] ≤ bn_pr
                resize!(αn, j)
                resize!(βn, j-1)
                verbose && println("✓ Converged at step $j")
                break
            end
            
            # Normalize next vector
            q_next = view(Q, :, j+1)
            q_next .= w ./ βn[j]
            
            if verbose && j % 50 == 0
                println("Step $j/$KDim: β=$(round(βn[j], sigdigits=4))")
            end
        end
    end
    
    return (α = αn, β = βn)
end




# # LANCZOS-FRO: FULL REORTHOGONALIZATION ALGORITHM
# function extract_lanczos_coefficients_fro(H::SparseMatrixCSC{ComplexF64, Int64}, ψ0::Vector{ComplexF64}, KDim::Int; verbose=false)

#     dim = length(ψ0)
#     bn_pr = sqrt(eps(Float64))  # Tolerance ≈ 1e-8
    
#     # Storage for all Lanczos vectors (for full reorthogonalization)
#     Q = Matrix{ComplexF64}(undef, dim, KDim + 1)
    
#     # Working vectors
#     w = Vector{ComplexF64}(undef, dim)
    
#     # Initial vector
#     Q[:, 1] = ψ0 / norm(ψ0)
    
#     # Pre-allocate coefficient arrays (β has one less element than α)
#     αn = Vector{Float64}(undef, KDim)
#     βn = Vector{Float64}(undef, KDim-1)
    
#     # Progress tracking
#     #verbose = KDim > 50
#     #verbose = false
    
#     # Main Lanczos-FRO loop
#     for j in 1:KDim
#         # Apply Hamiltonian: w = H * q_j
#         mul!(w, H, view(Q, :, j))
        
#         # Compute diagonal element: α_j = q_j^T * w
#         αn[j] = real(dot(view(Q, :, j), w))
        
#         # w = w - α_j * q_j
#         axpy!(-αn[j], view(Q, :, j), w)
        
#         # Full reorthogonalization: w = w - Q * (Q^T * w)
#         # This ensures w ⊥ span{q_1, ..., q_j}
#         if j > 1
#             # Compute overlaps: overlaps = Q[:, 1:j]^T * w
#             overlaps = Vector{ComplexF64}(undef, j)
#             for i in 1:j-1
#                 overlaps[i] = dot(view(Q, :, i), w)
#             end
            
#             # Remove overlaps: w = w - Q[:, 1:j] * overlaps
#             for i in 1:j-1
#                 axpy!(-overlaps[i], view(Q, :, i), w)
#             end
            
#             # Second Gram-Schmidt Step
#             for i in 1:j-1
#                 overlap = dot(view(Q, :, i), w)
#                 axpy!(-overlap, view(Q, :, i), w)
#             end
#         end
        
#         # Compute off-diagonal element (only if not the last step)
#         if j < KDim
#             βn[j] = real(norm(w))
            
#             # Check for convergence or breakdown
#             if βn[j] <= bn_pr
#                 resize!(αn, j)
#                 resize!(βn, j-1)
#                 verbose && println("\n✓ Lanczos-FRO converged at step $j (β = $(βn[j]))")
#                 break
#             end
            
#             # Normalize and store next vector: q_{j+1} = w / β_j
#             Q[:, j+1] = w / βn[j]
            
#             # Progress update
#             if verbose && j % 50 == 0
#                 println("Lanczos-FRO step $j/$KDim, β = $(round(βn[j], sigdigits=4)), α = $(round(αn[j], sigdigits=4))")
#             end
#         else
#             # Last step: just compute the final α, no β needed
#             verbose && println("\n✓ Lanczos-FRO completed all $j steps")
#         end
#     end
    
#     # Verify orthogonality (optional diagnostic)
#     if verbose && length(αn) > 1
#         n_computed = length(αn)
#         orthogonality_check = norm(Q[:, 1:n_computed]' * Q[:, 1:n_computed] - I(n_computed))
#         println("Orthogonality error: $(round(orthogonality_check, sigdigits=6))")
#     end

#     return (α = αn, β = βn)
# end