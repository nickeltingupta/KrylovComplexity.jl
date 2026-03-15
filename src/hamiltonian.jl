# Build the SpinXY4 Hamiltonian using bitwise operations - OPTIMIZED VERSION
function build_SpinXY4_bitwise_efficient(disorder_seed::Int, N::Int, validindices::Vector{NTuple{4, Int64}})

    Random.seed!(disorder_seed)

    ################################################################################################

    function pauli_x(s::Int, i::Int)
        return s ⊻ (1 << i)
    end

    function pauli_y(s::Int, i::Int)
        s_prime = s ⊻ (1 << i)
        phase = ((s >> i) & 1) == 1 ? -1.0im : 1.0im
        return (s_prime, phase)
    end

    function HFCalc(i, j, k, l)
        HF = 0
        (i & 1 == 1) && (j == i + 1) && (HF += 1)
        (j & 1 == 1) && (k == j + 1) && (HF += 1)
        (k & 1 == 1) && (l == k + 1) && (HF += 1)
        return HF
    end

    ################################################################################################ 
    
    terms_in_ham = length(validindices)
    HF_vals = [HFCalc(x...) for x in validindices]
    
    dim = 1 << N
    # Pre-allocate arrays for the sparse matrix construction.
    # The number of non-zero elements is exactly dim * terms_in_ham.
    # This is an upper bound; some might be zero but it's predictable.
    nnz = dim * terms_in_ham
    I_vec = Vector{Int}(undef, nnz) # Row indices
    J_vec = Vector{Int}(undef, nnz) # Column indices
    V_vec = Vector{ComplexF64}(undef, nnz) # Values

    J_vals = rand(Normal(0, 1), terms_in_ham)

    ################################################################################################    
    index_counter = 1
    for s in 0:dim-1
        for (term_idx, sites) in enumerate(validindices) 

            amp = 1.0 + 0.0im
            final_state = s

            for site in sites 
                if site & 1 == 1  # Odd index → σₓ
                    bit_pos = (site - 1) ÷ 2
                    final_state = pauli_x(final_state, bit_pos)
                else  # Even index → σᵧ
                    bit_pos = (site - 2) ÷ 2
                    final_state, phase = pauli_y(final_state, bit_pos)
                    amp *= phase
                end
            end
            
            # Instead of inserting into the matrix, store the data.
            I_vec[index_counter] = s + 1 # 1-based row index
            J_vec[index_counter] = final_state + 1 # 1-based column index
            # Calculate the value for this matrix element
            V_vec[index_counter] = J_vals[term_idx] * (1.0im^HF_vals[term_idx]) * amp
            
            index_counter += 1
        end
    end
    
    # Construct the Hamiltonian in one step.
    Ham = sparse(I_vec, J_vec, V_vec, dim, dim)
    
    return (Ham * sqrt(0.75 / N^3))
end


################################################################################################
################################################################################################
################################################################################################


# Build the SpinXY4 Hamiltonian using bitwise operations
function build_SpinXY4_bitwise_old(disorder_seed::Int, N::Int, validindices::Vector{NTuple{4, Int64}})

    # Use the worker-specific RNG
    #r = deepcopy(rng)  # Make a local copy to avoid contention
    #Random.seed!(r, disorder_seed)  # Further seed for this realization
    
    Random.seed!(disorder_seed)

    ################################################################################################

    function pauli_x(s::Int, i::Int)
        return s ⊻ (1 << i)
    end

    function pauli_y(s::Int, i::Int)
        s_prime = s ⊻ (1 << i)
        phase = ((s >> i) & 1) == 1 ? -1.0im : 1.0im
        return (s_prime, phase)
    end

    function HFCalc_optimized(i, j, k, l)
        HF = 0
        (i & 1 == 1) && (j == i + 1) && (HF += 1)
        (j & 1 == 1) && (k == j + 1) && (HF += 1)
        (k & 1 == 1) && (l == k + 1) && (HF += 1)
        return HF
    end

    ################################################################################################ 
    
    terms_in_ham = length(validindices)
    HF = [(x, HFCalc_optimized(x...)) for x in validindices];
    
    
    dim = 1 << N
    Ham = spzeros(ComplexF64, dim, dim)

    #J = rand(r, Normal(0, 1), terms_in_ham)
    J = rand(Normal(0, 1), terms_in_ham)

    ################################################################################################    

    for s in 0:dim-1
        for i in 1:terms_in_ham # keeping 1-based indexing because this has nothing to do with bitwise operations

            amp = 1.0 + 0.0im
            final_state = s

            for j in 1:4 # keeping 1-based indexing
                site = validindices[i][j]

                if site & 1 == 1  # Odd index → σₓ
                    bit_pos = (site - 1) ÷ 2  # '-1' is to make sure that the bitwise operations act normally starting from 0 rather than 1
                    final_state = pauli_x(final_state, bit_pos)
                else  # Even index → σᵧ
                    bit_pos = (site - 2) ÷ 2  # '-1' is to make sure that the bitwise operations act normally starting from 0 rather than 1
                    final_state, phase = pauli_y(final_state, bit_pos)
                    amp *= phase
                end
            end
            
            Ham[s+1, final_state+1] += J[i] * (1.0im^HF[i][2]) * amp # +1 to convert from 0-based to 1-based indexing
        end
    end
    
    return (Ham * sqrt(0.75 / N^3))
    #return zchop.(Ham,eps(Float64))
end


################################################################################################
################################################################################################
################################################################################################


# Build the SpinXY4 Hamiltonian without bitwise operations
function build_SpinXY4_regular(disorder_seed::Int, N::Int, validindices::Vector{NTuple{4, Int64}})


    # Use the worker-specific RNG
    r = deepcopy(rng)  # Make a local copy to avoid contention
    Random.seed!(r, disorder_seed)  # Further seed for this realization

    ################################################################################################

    function pauli_op()
        
        σx = [0.0 1.0; 1.0 0.0]
        σy = [0.0 -im; im 0.0]
        #σz = [1.0 0.0; 0.0 -1.0]

        return sparse(σx), sparse(σy)#, sparse(σz)
    end

    function lattice_term(N::Int, site::Int, interaction) #interactionsize = Int(log2(size(interaction)[1])) #number of paulis in the interaction, =2 for hopping and =1 for fields

        @assert site <= N
        left = sparse(I, 2^(site-1), 2^(site-1)) 
        right = sparse(I, 2^(N-site), 2^(N-site))

        return kron(left, interaction, right)
    end

    function HFCalc(i, j, k, l)
        HF = 0
        (i & 1 == 1) && (j == i + 1) && (HF += 1)
        (j & 1 == 1) && (k == j + 1) && (HF += 1)
        (k & 1 == 1) && (l == k + 1) && (HF += 1)
        return HF
    end

    ################################################################################################

    pauli_x, pauli_y = pauli_op()

    validindices = [(i, j, k, l) for (i, j, k, l) in combinations(1:2N, 4) if i < j < k < l]
    terms_in_ham = length(validindices)
    HF = [(x, HFCalc(x...)) for x in validindices];


    dim = 1 << N
    Ham = spzeros(ComplexF64, dim, dim)

    J = rand(r, Normal(0, 1), terms_in_ham)

    ################################################################################################

    for i in 1:terms_in_ham 

        dum = sparse(I, dim, dim)
        
        for j in 1:4 
            site = validindices[i][j]

            if site & 1 == 1  # Odd index → σₓ
                dum *= lattice_term(N, Int((site+1)/2), pauli_x)
            else  # Even index → σᵧ
                dum *= lattice_term(N, Int(site/2), pauli_y)
            end
        end
        
        Ham += J[i] * dum * (1.0im^HF[i][2])

    end


    return (Ham * sqrt(0.75 / N^3))
    #return zchop.(Ham,eps(Float64))
end