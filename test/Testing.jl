module Testing

    export examples

    struct ExampleData
        filename::String
        ss_composition::Dict{Symbol,Dict{String, Int}}
    end

    data_dir = joinpath(@__DIR__, "data")

    examples = [
        ExampleData(
            joinpath(data_dir, "pdb", "pdb1fmc.pdb"),
            Dict(
                :stride => Dict("bend" => 0, "kappa helix" => 0, "beta strand" => 77, "strand" => 81, "loop" => 0, "310 helix" => 21, "turn" => 70, "helix" => 263, "beta bridge" => 4, "alpha helix" => 242, "pi helix" => 0, "coil" => 96),
                :dssp => Dict("bend" => 26, "kappa helix" => 0, "beta strand" => 72, "strand" => 76, "loop" => 128, "310 helix" => 18, "turn" => 52, "helix" => 228, "beta bridge" => 4, "alpha helix" => 200, "pi helix" => 10, "coil" => 0),
            )
        ),
        ExampleData(
            joinpath(data_dir, "pdb", "pdb1mjh.pdb"),
            Dict(
                :stride => Dict("bend" => 0, "kappa helix" => 0, "beta strand" => 64, "strand" => 64, "loop" => 0, "310 helix" => 3, "turn" => 30, "helix" => 135, "beta bridge" => 0, "alpha helix" => 132, "pi helix" => 0, "coil" => 58), 
                :dssp =>  Dict("bend" => 20, "kappa helix" => 7, "beta strand" => 58, "strand" => 58, "loop" => 67, "310 helix" => 3, "turn" => 14, "helix" => 128, "beta bridge" => 0, "alpha helix" => 118, "pi helix" => 0, "coil" => 0),
            )
        ),
        ExampleData(
            joinpath(data_dir, "pdb", "pdb1res.pdb"),
            Dict(
                :stride => Dict("bend" => 0, "kappa helix" => 0, "beta strand" => 0, "strand" => 0, "loop" => 0, "310 helix" => 0, "turn" => 6, "helix" => 24, "beta bridge" => 0, "alpha helix" => 24, "pi helix" => 0, "coil" => 13),
                :dssp => Dict("bend" => 5, "kappa helix" => 0, "beta strand" => 0, "strand" => 0, "loop" => 11, "310 helix" => 0, "turn" => 6, "helix" => 21, "beta bridge" => 0, "alpha helix" => 21, "pi helix" => 0, "coil" => 0), 
            )
        ),
        ExampleData(
            joinpath(data_dir, "pdb", "pdb2btm.pdb"),
            Dict(
                :stride => Dict("bend" => 0, "kappa helix" => 0, "beta strand" => 75, "strand" => 79, "loop" => 2, "310 helix" => 18, "turn" => 77, "helix" => 251, "beta bridge" => 4, "alpha helix" => 233, "pi helix" => 0, "coil" => 93), 
                :dssp => Dict("bend" => 39, "kappa helix" => 11, "beta strand" => 77, "strand" => 81, "loop" => 108, "310 helix" => 18, "turn" => 46, "helix" => 228, "beta bridge" => 4, "alpha helix" => 199, "pi helix" => 0, "coil" => 0),
            )
        ),
    ]

end