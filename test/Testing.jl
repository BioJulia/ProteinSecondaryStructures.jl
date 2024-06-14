module Testing

export examples

struct ExampleData
    filename::String
    ss_composition::Dict{Symbol,Dict{String,Int}}
end

data_dir = joinpath(@__DIR__, "data")

examples = [
    ExampleData(# 1
        joinpath(data_dir, "pdb", "pdb1fmc.pdb"),
        Dict(
            :stride => Dict("310 helix" => 21, "bend" => 0, "turn" => 70, "beta bridge" => 4, "kappa helix" => 0, "pi helix" => 0, "beta strand" => 77, "alpha helix" => 242, "loop" => 0, "coil" => 96),
            :dssp => Dict("310 helix" => 20, "bend" => 26, "turn" => 58, "beta bridge" => 4, "kappa helix" => 0, "pi helix" => 10, "beta strand" => 74, "alpha helix" => 220, "loop" => 98, "coil" => 0)
        )
    ),
    ExampleData(# 2
        joinpath(data_dir, "pdb", "pdb1mjh.pdb"),
        Dict(
            :stride => Dict("310 helix" => 3, "bend" => 0, "turn" => 30, "beta bridge" => 0, "kappa helix" => 0, "pi helix" => 0, "beta strand" => 64, "alpha helix" => 132, "loop" => 0, "coil" => 58),
            :dssp => Dict("310 helix" => 3, "bend" => 20, "turn" => 14, "beta bridge" => 0, "kappa helix" => 7, "pi helix" => 0, "beta strand" => 64, "alpha helix" => 126, "loop" => 53, "coil" => 0),
        )
    ),
    ExampleData(# 3
        joinpath(data_dir, "pdb", "pdb1res.pdb"),
        Dict(
            :stride => Dict("310 helix" => 0, "bend" => 0, "turn" => 6, "beta bridge" => 0, "kappa helix" => 0, "pi helix" => 0, "beta strand" => 0, "alpha helix" => 24, "loop" => 0, "coil" => 13),
            :dssp => Dict("310 helix" => 0, "bend" => 5, "turn" => 6, "beta bridge" => 0, "kappa helix" => 0, "pi helix" => 0, "beta strand" => 0, "alpha helix" => 23, "loop" => 9, "coil" => 0),
        )

    ),
    ExampleData(# 4
        joinpath(data_dir, "pdb", "pdb2btm.pdb"),
        Dict(
            :stride => Dict("310 helix" => 18, "bend" => 0, "turn" => 77, "beta bridge" => 4, "kappa helix" => 0, "pi helix" => 0, "beta strand" => 75, "alpha helix" => 233, "loop" => 0, "coil" => 93),
            :dssp => Dict("310 helix" => 18, "bend" => 43, "turn" => 54, "beta bridge" => 4, "kappa helix" => 11, "pi helix" => 0, "beta strand" => 77, "alpha helix" => 213, "loop" => 80, "coil" => 0),
        )
    ),
    ExampleData(# 5
        joinpath(data_dir, "pdb", "bad.pdb"),
        Dict(
            :stride => Dict("310 helix" => 0, "bend" => 0, "turn" => 0, "beta bridge" => 0, "kappa helix" => 0, "pi helix" => 0, "beta strand" => 0, "alpha helix" => 0, "loop" => 0, "coil" => 0),
            :dssp => Dict("310 helix" => 0, "bend" => 0, "turn" => 0, "beta bridge" => 0, "kappa helix" => 0, "pi helix" => 0, "beta strand" => 0, "alpha helix" => 0, "loop" => 0, "coil" => 0),
        )
    ),
    ExampleData(# 6
        joinpath(data_dir, "pdb", "stride_b.pdb"),
        Dict(
            :stride => Dict("310 helix" => 0, "bend" => 0, "turn" => 2, "beta bridge" => 5, "kappa helix" => 0, "pi helix" => 0, "beta strand" => 0, "alpha helix" => 0, "loop" => 0, "coil" => 8),
            :dssp => Dict("310 helix" => 0, "bend" => 1, "turn" => 2, "beta bridge" => 2, "kappa helix" => 3, "pi helix" => 0, "beta strand" => 0, "alpha helix" => 0, "loop" => 6, "coil" => 0),
        )
    ),
    ExampleData(# 7
        joinpath(data_dir, "cif", "6e55.cif"),
        Dict(
            :dssp => Dict("310 helix" => 9, "bend" => 21, "turn" => 69, "beta bridge" => 2, "kappa helix" => 1, "pi helix" => 0, "beta strand" => 162, "alpha helix" => 97, "loop" => 119, "coil" => 0)
        )
    ),
]

end # module Testing