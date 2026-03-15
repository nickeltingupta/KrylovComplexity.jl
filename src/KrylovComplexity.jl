module KrylovComplexity

using Combinatorics, DataFrames, Distributions, DelimitedFiles, LaTeXStrings, LinearAlgebra, LinearAlgebra.BLAS, Plots, Polynomials, PrettyTables, ProgressBars, ProgressMeter, Random, SparseArrays, SpecialFunctions, Statistics, StatsBase, Test, ZChop

# Export public functions
export build_SpinXY4_bitwise_efficient, build_SpinXY4_bitwise_old, build_SpinXY4_regular, build_parity_blocks, level_spacing_plot, extract_lanczos_coefficients_fro_blas, build_tfd_state, compute_krylov_complexity

# Include all the files
include("hamiltonian.jl")
include("block_diagonalize.jl")
include("level_spacing.jl")
include("lanczos.jl")
include("TFD.jl")
include("spreadcomplexity.jl")

end