function build_parity_blocks(H::AbstractMatrix{ComplexF64}, N::Int)

################################################################################################

    function separate_parity_states(N::Int)
        dim = 1 << N
        dim_half = dim ÷ 2  # Each parity sector has exactly dim/2 states
        
        indices_plus = Vector{Int}(undef, dim_half)
        indices_minus = Vector{Int}(undef, dim_half)
        
        idx_plus = 1
        idx_minus = 1
        
        for s in 0:dim-1
            if count_ones(s) % 2 == 0
                indices_plus[idx_plus] = s
                idx_plus += 1
            else
                indices_minus[idx_minus] = s
                idx_minus += 1
            end
        end
        
        return indices_plus, indices_minus
    end

################################################################################################

    indices_plus, indices_minus = separate_parity_states(N)

    H_plus = H[indices_plus .+ 1, indices_plus .+ 1]
    H_minus = H[indices_minus .+ 1, indices_minus .+ 1]
    
################################################################################################
    
    E_plus, V_plus   = eigen(Hermitian(Matrix(H_plus)))
    E_minus, V_minus = eigen(Hermitian(Matrix(H_minus)))

    return (
        H_plus = H_plus, 
        H_minus = H_minus, 
        #indices_plus = indices_plus, 
        #indices_minus = indices_minus,
        E_plus  = E_plus,  V_plus  = V_plus,
        E_minus = E_minus, V_minus = V_minus
    )
end


###############################################################################################
################################################################################################
################################################################################################


# # #Constructing the parity operator (simple way, this is not the most efficient way because memory needs to be assigned and so forth)
# # Γ = spzeros(Float64, 1 << N, 1 << N)
# # for s in 0:(1 << N) - 1
# #     Γ[s+1, s+1] = count_ones(s) & 1 == 1 ? -1.0 : 1.0
# # end

# # Leveraging the fact that (-1)^k = 1 if k is even, -1 if k is odd
# Γ = spdiagm([(-1.0)^count_ones(s) for s in 0:(1<<N)-1]);

# @test Γ * Γ == sparse(I, 1<<N, 1<<N) # Parity operator squares to identity
# @test all(Γ * results[i].H == results[i].H * Γ for i in 1:n_disorders) # Parity operator commutes with Hamiltonian

# # Build parity eigenstates corresponding to Γ for a single disorder realization using Hamiltonian and Parity Operator - gets very expensive for large disorders and N
# function build_parity_eigenstates_usingH(H, Γ)
#     eigenvalues, eigenvectors = eigen(Matrix(H))

#     # Parity eigenvectors should satisfy Γ * |ψ⟩ = ±|ψ⟩, can also test using dot product
#     parity_eigenvector_check = collect(Γ * eigenvectors[:,i] ≈ 1.0 * eigenvectors[:,i] || Γ * eigenvectors[:,i] ≈ -1.0 * eigenvectors[:,i] for i in 1:1<<N)
#     @test count(x -> x == 0, parity_eigenvector_check) == 0 # Counts and tests the number of eigenvectors that are not parity eigenvectors, should be 0

#     # Calculate parity eigenvalues for each eigenvector
#     parity_eigenvalues = real.(([dot(eigenvectors[:,i], Γ * eigenvectors[:,i]) for i in 1:1<<N]))

#     # Separate based on parity eigenvalue
#     plus_indices = findall(x -> x ≈ 1.0, parity_eigenvalues)  # +1 eigenvalue (even parity)
#     minus_indices = findall(x -> x ≈ -1.0, parity_eigenvalues)   # -1 eigenvalue (odd parity)

#     @test length(plus_indices)+length(minus_indices) == 1<<N # Check that all eigenvectors are accounted for

#     parity_plus = eigenvectors[:, plus_indices]
#     parity_minus = eigenvectors[:, minus_indices]

#     return (
#         parity_plus = parity_plus, 
#         parity_minus = parity_minus,
#         plus_indices = plus_indices,
#         minus_indices = minus_indices
#         )
# end
