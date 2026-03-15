using Documenter
using KrylovComplexity

makedocs(
    sitename  = "KrylovComplexity.jl",
    authors   = "Nitin",
    modules   = [KrylovComplexity],
    format    = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = "https://nickeltingupta/KrylovComplexity.jl",
    ),
    pages = [
        "Home"    => "index.md",
        "Physics" => "physics.md",
        "Examples"=> "examples.md",
    ],
    checkdocs = :none,
    warnonly  = true,
)

deploydocs(
    repo   = "github.com/nickeltingupta/KrylovComplexity.jl.git",
    branch = "gh-pages",
    push_preview = true,
)