function build_tfd_state(eigenvals::Vector{Float64}, eigenvecs::Matrix{ComplexF64}, β::Float64=0.0)
    """
    Build TFD state: |TFD⟩ = Σ_n e^(-βE_n/2) |n⟩ / √Z
    For β=0: equal superposition of all energy eigenstates
    """
    
    if β == 0.0
        # Equal superposition (infinite temperature)
        coeffs = ones(length(eigenvals)) / sqrt(length(eigenvals))
    else
        # Thermal weights
        weights = exp.(-β * eigenvals / 2.0)
        Z = sum(weights.^2)  # Partition function
        coeffs = weights / sqrt(Z)
    end
    
    #println("The Partition function Z = $(round(Z, sigdigits=4))\n\n\n")

    # Build TFD state in original basis
    #tfd_state = eigenvecs * coeffs

    return eigenvecs * coeffs
end