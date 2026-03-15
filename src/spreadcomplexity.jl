function compute_krylov_complexity(αn::Vector{Float64}, βn::Vector{Float64}, times::Vector{Float64}; max_dim::Int=length(αn), show_progress=false)
    """
    Efficient Krylov complexity computation using eigendecomposition
    C_K(t) = ∑_{n=0}^{K-1} n * |c_n(t)|^2

    Optimization: Precompute the matrix U_{n,k} = eigenvecs[n, k] * conj(eigenvecs[1, k])
    This reduces the inner amplitude calculation to a simple vector-vector product.
    """
    K = min(max_dim, length(αn))
    βn_trunc = length(βn) >= K-1 ? βn[1:K-1] : βn

    # 1. Build and diagonalize the tridiagonal matrix (O(K^2) to O(K^3))
    T = Tridiagonal(βn_trunc, αn[1:K], βn_trunc) # More efficient than building a full matrix
    eigen_decomp = eigen(T)
    eigenvals = eigen_decomp.values
    eigenvecs = eigen_decomp.vectors # K x K matrix

    # 2. PRECOMPUTATION: This is the key step (O(K^2))
    # Precompute U[n, k] = eigenvecs[n, k] * conj(eigenvecs[1, k])
    # This matrix captures the initial condition and its projection onto all eigenvectors.
    U = eigenvecs .* conj.(eigenvecs[1, :]') # (K x K) .* (1 x K) via broadcasting -> K x K

    # 3. Precompute the time evolution factors for all k and all t (O(K * T))
    # This creates a matrix `exp_vals` of size (length(times), K)
    # where exp_vals[t_idx, k] = exp(-1im * eigenvals[k] * times[t_idx])
    exp_vals = Matrix{ComplexF64}(undef, length(times), K)
    for (t_idx, t) in enumerate(times)
        for k in 1:K
            exp_vals[t_idx, k] = exp(-1im * eigenvals[k] * t)
        end
    end

    # 4. Compute the amplitudes and complexity for all times (O(K^2 * T))
    # This is still O(K^2 * T) but the inner loop is now a fast BLAS operation.
    complexity = zeros(Float64, length(times))

    # Precompute the 'n' weights for the sum
    n_weights = collect(0:(K-1))

    progress_iter = show_progress ? ProgressBar(1:length(times)) : 1:length(times)
    for t_idx in progress_iter
        # For this time point, get the vector of exp(-iλt) for all k
        exp_vec_t = @view exp_vals[t_idx, :]

        # Calculate the amplitude for all n simultaneously.
        # This is a matrix-vector multiplication: amplitude_n = U * exp_vec_t
        # But we need to do U_{n,k} * exp_vec_t[k], which is exactly:
        amplitude_n = U * exp_vec_t # This is a K-length vector

        # Now compute |c_n(t)|^2 for all n
        cn_squared = abs2.(amplitude_n)

        # Compute complexity: dot product with the weights [0, 1, 2, ..., K-1]
        complexity[t_idx] = dot(n_weights, cn_squared)

        if show_progress && t_idx % 100 == 0
            set_description(progress_iter, "Time Evolution $t_idx/$(length(times)), t=$(round(times[t_idx], sigdigits=3))")
        end
    end

    return complexity
end


# function compute_krylov_complexity_old(αn::Vector{Float64}, βn::Vector{Float64}, times::Vector{Float64}; max_dim::Int=length(αn), show_progress=false)
#     """
#     Efficient Krylov complexity computation using eigendecomposition
#     C_K(t) = ∑_{n=0}^{K-1} n * |c_n(t)|^2
#     """
#     K = min(max_dim, length(αn))
#     βn_trunc = length(βn) >= K-1 ? βn[1:K-1] : βn
    
#     #println("Building tridiagonal matrix (K=$K)...")
#     # Build tridiagonal matrix
#     T = zeros(Float64, K, K)
#     for i in 1:K
#         T[i,i] = αn[i]
#     end
#     for i in 1:length(βn_trunc)
#         T[i, i+1] = βn_trunc[i]
#         T[i+1, i] = βn_trunc[i]
#     end
    
#     #println("Computing eigendecomposition...")
#     # Eigendecomposition: T = Q * Λ * Q†
#     eigenvals, eigenvecs = eigen(T)
    
#     complexity = zeros(Float64, length(times))
    
#     #println("Computing time evolution...")
#     # Add progress tracking for time evolution
#     #show_progress = length(times) > 100
#     #show_progress = false
#     progress_iter = show_progress ? ProgressBar(enumerate(times)) : enumerate(times)
    
#     for (t_idx, t) in progress_iter
#         # |c_n(t)|^2 = |∑_k Q[n,k] * exp(-i*λ_k*t) * Q*[k,0]|^2
#         cn_squared = zeros(Float64, K)
        
#         for n in 1:K
#             amplitude = sum(eigenvecs[n, k] * exp(-1im * eigenvals[k] * t) * conj(eigenvecs[1, k]) for k in 1:K)
#             cn_squared[n] = abs2(amplitude)
#         end
        
#         # Compute complexity
#         complexity[t_idx] = sum((n-1) * cn_squared[n] for n in 1:K)
        
#         # Update progress every 100 time points
#         if show_progress && t_idx % 100 == 0
#             set_description(progress_iter, "Time Evolution $t_idx/$(length(times)), t=$(round(t, sigdigits=3))")
#         end
#     end
    
#     return complexity #, T, eigenvals, eigenvecs
# end