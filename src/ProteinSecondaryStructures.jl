module ProteinSecondaryStructures

using TestItems: @testitem

# Basic data structure for each residue
export SSData

# Secondary structure calculation methods
export stride_run
export dssp_run

# Secondary structure summary
export ss_composition

# Representation conversion functions
export ss_name
export ss_code
export ss_number

"""
    SSData

A struct to hold secondary structure data for a single residue. The fields are:

* `resname::String` - residue name
* `chain::String` - chain identifier
* `resnum::Int` - residue number
* `sscode::String` - secondary structure code
* `phi::Float64` - phi dihedral angle
* `psi::Float64` - psi dihedral angle

"""
struct SSData
    resname::String
    chain::String
    resnum::Int
    sscode::String # as String for good reasons*
    phi::Float64
    psi::Float64
end
#=
*We don´t want the user to have an error if initializing this code
with a singe-character string, or using any of the other functions,
like `ss_name` with a one-character string. It is easier just
to support strings here and convert characters to strings when 
needed. Also, the Char type is printed in a too-technical way. 
And, finally, we may want to add more secondary structure types
represented by more than one character in the future.
=#

# Add proper header and empty chain identifier to the pdb file,
# for each secondary structure calculation method
function adjust_pdb_file(pdbfile; header::AbstractString, empty_chain_identifier::AbstractString)
    tmp_file = tempname() * ".pdb"
    open(tmp_file, "w") do io
        println(io, header)
        for line in readlines(pdbfile)
            if startswith(line, "ATOM")
                if line[22] == ' '
                    line = line[1:21]*empty_chain_identifier*line[23:end]
                end
                println(io, line)
            end
        end
        print(io, "END")
    end
    return tmp_file
end

# Stride interface
include("./stride.jl")

# DSSP interface
include("./dssp.jl")

# Common interface 
#=
    ss_classes

Secondary structure classes dictionary. 

=#
const ss_classes = Dict{String,String}(
    "G" => "310 helix",
    "H" => "alpha helix",
    "I" => "pi helix",
    "P" => "kappa helix",
    "T" => "turn",
    "E" => "beta strand",
    "B" => "beta bridge",
    "S" => "bend",
    "C" => "coil",
    " " => "loop",
)
const name_to_code = Dict{String,String}(ss_classes[key] => key for key in keys(ss_classes))

const code_to_number = Dict{String,Int}(
    "G" => 1, # 310 helix
    "H" => 2, # alpha helix
    "I" => 3, # pi helix
    "P" => 4, # kappa helix
    "T" => 5, # turn
    "E" => 6, # beta strand
    "B" => 7, # beta bridge
    "S" => 8, # bend
    "C" => 9, # coil
    " " => 10, # loop
)
const number_to_code = Dict{Int,String}(code_to_number[key] => key for key in keys(code_to_number))

"""
    ss_number(code::Union{SSData,AbstractString,AbstractChar})

Returns the secondary structure number code. The input may be a
secondary structure `String` code, a secondary structure name
(`"310 helix"`, `"alpha helix"`, ..., `"coil"`), or a `SSData` object.

The classification follows the DSSP standard classes, described 
in the ProteinSecondaryStructures.jl documentation.

# Example

```jldoctest
julia> using ProteinSecondaryStructures

julia> ss_number("H")
2

julia> ss_number('B')
7

julia> ss_number("beta bridge")
7

julia> ss = SSData("ARG", "A", 1, "H", 0.0, 0.0)
SSData("ARG", "A", 1, "H", 0.0, 0.0)

julia> ss_number(ss)
2

```

"""
ss_number(code::SSData) = code_to_number[code.sscode]
ss_number(code::AbstractChar) = code_to_number[string(code)]
function ss_number(code::AbstractString) 
    if length(code) == 1
        return code_to_number[string(code)]
    else
        return code_to_number[name_to_code[code]]
    end
end

@testitem "ss_number" begin
    @test ss_number("H") == 2
    @test ss_number('B') == 7
    @test ss_number("beta bridge") == 7
    @test ss_number(SSData("ARG", "A", 1, "H", 0.0, 0.0)) == 2
end

"""
    ss_code(code::Union{SSData,String,Integer})

Returns the one-letter secondary structure code. The input may be a
secondary structure `Integer` code, a secondary structure name
(`"310 helix"`, `"alpha helix"`, ..., `"coil"`), or a `SSData` object.

The classification follows the DSSP standard classes, described 
in the ProteinSecondaryStructures.jl documentation.

# Example

```jldoctest
julia> using ProteinSecondaryStructures

julia> ss_code(2)
"H"

julia> ss_code("beta bridge")
"B"

julia> ss = SSData("ARG", "A", 1, "H", 0.0, 0.0)
SSData("ARG", "A", 1, "H", 0.0, 0.0)

julia> ss_code(ss)
"H"

```

"""
ss_code(code::Integer) = number_to_code[code]
ss_code(code::AbstractString) = name_to_code[code]
ss_code(code::SSData) = code.sscode 

@testitem "ss_code" begin
    @test ss_code(2) == "H"
    @test ss_code("beta bridge") == "B"
    @test ss_code(SSData("ARG", "A", 1, "H", 0.0, 0.0)) == "H"
end

"""
    ss_name(ss::Union{SSData, Integer, String, Char})

Return the secondary structure name. The input may be a `SSData` object, 
a secondary structure `Integer` code (1-10) or a secondary 
structure code (`G, H, ..., C`).

The classification follows the DSSP standard classes, described 
in the ProteinSecondaryStructures.jl documentation.

## Example

```jldoctest
julia> using ProteinSecondaryStructures

julia> ss_name("H")
"alpha helix"

julia> ss_name(1)
"310 helix"

julia> ss = SSData("ARG", "A", 1, "H", 0.0, 0.0)
SSData("ARG", "A", 1, "H", 0.0, 0.0)

julia> ss_name(ss)
"alpha helix"
```

"""
ss_name(ss::SSData) = ss_classes[ss.sscode]
ss_name(code_number::Int) = ss_classes[number_to_code[code_number]]
ss_name(sscode::Union{AbstractString,AbstractChar}) = ss_classes[string(sscode)]

"""
    ss_composition(data::AbstractVector{<:SSData})

Calculate the secondary structure composition of the data. Returns a dictionary of
the secondary structure types and their counts.

"""
function ss_composition(data::AbstractVector{<:SSData})
    sscomposition = Dict{String,Int}()
    for sscode in keys(ss_classes)
        sscomposition[ss_name(sscode)] = count(ss -> ss.sscode == sscode, data)
    end
    return sscomposition
end

# Testing module
include("../test/Testing.jl")

end
