function level_spacing_plot(eigenvalues; num_bins=10, verbose=false)

    RandMat_eigenvals = sort(real(eigenvalues))
    
        
    # Computing level spacing ration to see wakanda stats we have
    rstat(evals) = begin
        s = diff(sort(real(evals)))
        mean(min.(s[1:end-1], s[2:end]) ./ max.(s[1:end-1], s[2:end]))
    end
        
    r = rstat(RandMat_eigenvals)

    (verbose) && (println("Eigenvalue range: [$(RandMat_eigenvals[1]), $(RandMat_eigenvals[end])]"))

    # Spectrum unfolding - simple
    RandMat_unfolded = (RandMat_eigenvals .- RandMat_eigenvals[1]) ./ (RandMat_eigenvals[end] - RandMat_eigenvals[1]) .* length(RandMat_eigenvals)
    
    # Spectrum unfolding - polynomial fit
    # fit = Polynomials.fit(RandMat_eigenvals, collect(1:length(RandMat_eigenvals)), 3)
    # RandMat_unfolded = fit.(RandMat_eigenvals)
    
    #RandMat_spacings = diff(RandMat_eigenvals)

    RandMat_spacings = diff(RandMat_unfolded)
    (verbose) && (println("Number of spacings: $(length(RandMat_spacings))"))
    
    keep_spacings = round(Int, 0.9*length(RandMat_spacings))
    (verbose) && (println("Trying to keep: $keep_spacings spacings"))

    truncated_RandMat_spacings = (sort(RandMat_spacings))[1:keep_spacings]
    truncated_RandMat_spacings = truncated_RandMat_spacings / mean(truncated_RandMat_spacings) # Normalize spacings

    (verbose) && (println("Spacings range: [$(truncated_RandMat_spacings[1]), $(truncated_RandMat_spacings[end])]"))

    p = histogram(
        truncated_RandMat_spacings,
        bins=num_bins,
        normalize=:pdf,
        label="Level spacings, r = $(round(r,  digits=4))", 
        title="Nearest Neighbor Level Spacings",
        framestyle=:box)

    s_theory = 0:0.01:3.0
    wigner_dyson_goe = (π/2) .* s_theory .* exp.(-π/4 .* s_theory.^2)
    wigner_dyson_gue = (32/π^2) .* s_theory.^2 .* exp.(-4 .* s_theory.^2 ./ π)

    plot!(p, s_theory, wigner_dyson_gue, linewidth=2, label=L"Wigner-Dyson-GUE, $\textbf{r}_{\textbf{th.}} \mathbf{≈}\ \mathbf{0.6027}$", color=:red)
    plot!(p, s_theory, wigner_dyson_goe, linewidth=2, label=L"Wigner-Dyson-GOE, $\textbf{r}_{\textbf{th.}} \mathbf{≈}\ \mathbf{0.5359}$", color=:green)

    return p
end