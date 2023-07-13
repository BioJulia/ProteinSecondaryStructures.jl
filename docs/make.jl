import Pkg
Pkg.add("Documenter")
using Documenter
using ProteinSecondaryStructures 
push!(LOAD_PATH,"../src/")
makedocs(
    modules=[ProteinSecondaryStructures],
    sitename="ProteinSecondaryStructures.jl",
    pages = [
        "Home" => "index.md",
        "Explanation" => "explanation.md",
        "Single PDB files" => "single_pdb.md",
        "MD Trajectories" => "trajectories.md",
        "Reference" => "reference.md",
        "Citations" => "citations.md",
    ]
)
deploydocs(
    repo = "github.com/m3g/ProteinSecondaryStructures.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "main", 
    versions = ["stable" => "v^", "v#.#" ],
)
