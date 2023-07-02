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
                :stride => Dict("bend" => 0, "kapa helix" => 0, "beta strand" => 77, "strand" => 81, "loop" => 0, "310 helix" => 21, "turn" => 70, "helix" => 263, "beta bridge" => 4, "alpha helix" => 242, "pi helix" => 0, "coil" => 96),
                :dssp => Dict("bend" => 26, "kapa helix" => 0, "beta strand" => 74, "strand" => 78, "loop" => 98, "310 helix" => 20, "turn" => 58, "helix" => 250, "beta bridge" => 4, "alpha helix" => 220, "pi helix" => 10, "coil" => 0),
            )
        ),
        ExampleData(
            joinpath(data_dir, "pdb", "pdb1mjh.pdb"),
            Dict(
                :stride => Dict("bend" => 0, "kapa helix" => 0, "beta strand" => 64, "strand" => 64, "loop" => 0, "310 helix" => 3, "turn" => 30, "helix" => 135, "beta bridge" => 0, "alpha helix" => 132, "pi helix" => 0, "coil" => 58), 
                :dssp =>  Dict("bend" => 22, "kapa helix" => 7, "beta strand" => 64, "strand" => 64, "loop" => 51, "310 helix" => 3, "turn" => 14, "helix" => 136, "beta bridge" => 0, "alpha helix" => 126, "pi helix" => 0, "coil" => 0),
            )
        ),
        ExampleData(
            joinpath(data_dir, "pdb", "pdb1res.pdb"),
            Dict(
                :stride => Dict("bend" => 0, "kapa helix" => 0, "beta strand" => 0, "strand" => 0, "loop" => 0, "310 helix" => 0, "turn" => 6, "helix" => 24, "beta bridge" => 0, "alpha helix" => 24, "pi helix" => 0, "coil" => 13),
                :dssp => Dict("bend" => 5, "kapa helix" => 0, "beta strand" => 0, "strand" => 0, "loop" => 9, "310 helix" => 0, "turn" => 6, "helix" => 23, "beta bridge" => 0, "alpha helix" => 23, "pi helix" => 0, "coil" => 0), 
            )
        ),
        ExampleData(
            joinpath(data_dir, "pdb", "pdb2btm.pdb"),
            Dict(
                :stride => Dict("bend" => 0, "kapa helix" => 0, "beta strand" => 75, "strand" => 79, "loop" => 2, "310 helix" => 18, "turn" => 77, "helix" => 251, "beta bridge" => 4, "alpha helix" => 233, "pi helix" => 0, "coil" => 93), 
                :dssp => Dict("bend" => 43, "kapa helix" => 11, "beta strand" => 77, "strand" => 81, "loop" => 80, "310 helix" => 18, "turn" => 54, "helix" => 242, "beta bridge" => 4, "alpha helix" => 213, "pi helix" => 0, "coil" => 0),
            )
        ),
    ]

end