import Pkg
Pkg.add("Documenter")
using Documenter
using ProteinSecondaryStructures 
push!(LOAD_PATH,"../src/")
makedocs(
    modules=[ProteinSecondaryStructures],
    sitename="ProteinSecondaryStructures.jl",
    warnonly = [:missing_docs],
    pages = [
        "Home" => "index.md",
        "Overview" => "overview.md",
        "User guide" => "user_guide.md",
    ]
)
deploydocs(
    repo = "github.com/BioJulia/ProteinSecondaryStructures.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "main", 
    versions = ["stable" => "v^", "v#.#" ],
)
