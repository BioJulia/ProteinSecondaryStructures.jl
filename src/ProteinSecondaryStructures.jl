module ProteinSecondaryStructures

using TestItems: @testitem
using ProgressMeter: Progress, next!
import Chemfiles
import PDBTools

export SSData

export stride_run
export dssp_run

export is_anyhelix, is_alphahelix, is_pihelix, is_310helix, is_kappahelix
export is_anystrand, is_betastrand, is_betabridge
export is_bend, is_coil, is_turn, is_loop
export ss_composition

# Representation conversion functions
export ss_name
export ss_code
export ss_number

import Base: @kwdef # for 1.6 compatibility

"""
    SSData

A struct to hold secondary structure data for a single residue. The fields are:

* `resname::String` - residue name
* `chain::String` - chain identifier
* `resnum::Int` - residue number
* `sscode::String` - secondary structure code
* `phi::Float64` - phi dihedral angle
* `psi::Float64` - psi dihedral angle
* `area::Float64` - solvent accessible area (stride specific)
* `kappa::Float64` - virtual bond angle (dssp specific)
* `alpha::Float64` - virtual torsion angle (dssp specific)

"""
@kwdef struct SSData
    resname::String
    chain::String
    resnum::Int
    sscode::String # as String because Char is `show`n uglier
    phi::Float64
    psi::Float64
    area::Float64 = 0.0 # stride specific
    kappa::Float64 = 0.0 # dssp specific
    alpha::Float64 = 0.0 # dssp specific
end

# Constructor for SSData from residue data only 
function SSData(resname::String, chain::String, resnum::Int)
    return SSData(resname, chain, resnum, " ", 0.0, 0.0, 0.0, 0.0, 0.0)
end

# Check if two residues are the same given the SSData fields
function residues_match(ss1::SSData, ss2::SSData)
    return ss1.resname == ss2.resname && ss1.chain == ss2.chain && ss1.resnum == ss2.resnum
end

# Stride interface
include("./stride.jl")

# DSSP interface
include("./dssp.jl")

# Common interface 
"""
    ss_classes

Secondary structure classes dictionary. 

"""
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

The classification follows the DSSP standard classes:

https://m3g.github.io/ProteinSecondaryStructures.jl/stable/explanation/#Secondary-structure-classes

# Example

```jldoctest
julia> using ProteinSecondaryStructures

julia> ss_number("H")
2

julia> ss_number('B')
7

julia> ss_number("beta bridge")
7

julia> ss = SSData("ARG", "A", 1, "H", 0.0, 0.0, 0.0, 0.0, 0.0);
SSData("ARG", "A", 1, "H", 0.0, 0.0, 0.0, 0.0, 0.0)

julia> ss_number(ss)
2

```

!!! compat
    This function was added in version 1.5.0

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
    @test ss_number(SSData("ARG", "A", 1, "H", 0.0, 0.0, 0.0, 0.0, 0.0)) == 2
end

"""
    ss_code(code::Union{String,Integer,SSData})

Returns the one-letter secondary structure code. The input may be a
secondary structure `Integer` code, a secondary structure name
(`"310 helix"`, `"alpha helix"`, ..., `"coil"`), or a `SSData` object.

The classification follows the DSSP standard classes:

https://m3g.github.io/ProteinSecondaryStructures.jl/stable/explanation/#Secondary-structure-classes

# Example

```jldoctest
julia> using ProteinSecondaryStructures

julia> ss_code(2)
"H"

julia> ss_code("beta bridge")
"B"

julia> ss = SSData("ARG", "A", 1, "H", 0.0, 0.0, 0.0, 0.0, 0.0)
SSData("ARG", "A", 1, "H", 0.0, 0.0, 0.0, 0.0, 0.0)

julia> ss_code(ss)
"H"

```

!!! compat
    This function was added in version 1.5.0

"""
ss_code(code::Integer) = number_to_code[code]
ss_code(code::AbstractString) = name_to_code[code]
ss_code(code::SSData) = code.sscode 

@testitem "ss_code" begin
    @test ss_code(2) == "H"
    @test ss_code("beta bridge") == "B"
    @test ss_code(SSData("ARG", "A", 1, "H", 0.0, 0.0, 0.0, 0.0, 0.0)) == "H"
end

"""
    ss_name(ss::Union{SSData, Integer, String, Char})

Return the secondary structure name. The input may be a `SSData` object, 
a secondary structure `Integer` code (1-8) or a secondary 
structure code (`G, H, ..., C`).

The secondary structure classes names and codes are:

| SS name             | `ss code`    | `code number`|
|:--------------------|:------------:|:------------:|
| `"310 helix"`       | `"G"`        | `1`          | 
| `"alpha helix"`     | `"H"`        | `2`          |
| `"pi helix"`        | `"I"`        | `3`          |
| `"kappa helix"`     | `"P"`        | `4`          |
| `"turn"`            | `"T"`        | `5`          |
| `"beta strand"`     | `"E"`        | `6`          |
| `"beta bridge"`     | `"B"`        | `7`          |
| `"bend"`            | `"S"`        | `8`          |
| `"coil"`            | `"C"`        | `9`          |
| `"loop"`            | `" "`        | `10`         |

## Example

```jldoctest
julia> using ProteinSecondaryStructures

julia> ss_name("H")
"alpha helix"

julia> ss_name(1)
"310 helix"

julia> ss = SSData("ARG", "A", 1, "H", 0.0, 0.0, 0.0, 0.0, 0.0)
SSData("ARG", "A", 1, "H", 0.0, 0.0, 0.0, 0.0, 0.0)

julia> ss_name(ss)
"alpha helix"
```

"""
ss_name(ss::SSData) = ss_classes[ss.sscode]
ss_name(code_number::Int) = ss_classes[number_to_code[code_number]]
ss_name(sscode::Union{AbstractString,AbstractChar}) = ss_classes[string(sscode)]

"""
    is_anyhelix(ss::SSData)
    is_alphahelix(ss::SSData)
    is_pihelix(ss::SSData)
    is_kappahelix(ss::SSData)
    is_310helix(ss::SSData)
    is_anystrand(ss::SSData)
    is_betastrand(ss::SSData)
    is_betabridge(ss::SSData)
    is_turn(ss::SSData)
    is_bend(ss::SSData)
    is_coil(ss::SSData)
    is_turn(ss::SSData)

Return `true` if the data is of the given secondary structure type.

"""
function _is_class end

@doc (@doc _is_class) is_anyhelix(ss::SSData) = ss.sscode in ("H", "G", "I", "P")
@doc (@doc _is_class) is_alphahelix(ss::SSData) = ss.sscode == "H"
@doc (@doc _is_class) is_pihelix(ss::SSData) = ss.sscode == "I"
@doc (@doc _is_class) is_kappahelix(ss::SSData) = ss.sscode == "P"
@doc (@doc _is_class) is_310helix(ss::SSData) = ss.sscode == "G"
@doc (@doc _is_class) is_anystrand(ss::SSData) = ss.sscode in ("E", "B")
@doc (@doc _is_class) is_betastrand(ss::SSData) = ss.sscode == "E"
@doc (@doc _is_class) is_betabridge(ss::SSData) = ss.sscode == "B"
@doc (@doc _is_class) is_turn(ss::SSData) = ss.sscode == "T"
@doc (@doc _is_class) is_bend(ss::SSData) = ss.sscode == "S"
@doc (@doc _is_class) is_coil(ss::SSData) = ss.sscode == "C"
@doc (@doc _is_class) is_loop(ss::SSData) = ss.sscode == " "

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


#
# The next will be deprecated in 2.0
#

# Trajectory analysis
include("./trajectories.jl")
const class = ss_name 
const ss_code_to_number = ss_number
const ss_number_to_code = ss_code
export ss_classes
export class, ss_code_to_number, ss_number_to_code
export ss_content
export ss_map
@testitem "ss_code_to_number/ss_number_to_code" begin
    @test ss_code_to_number('H') == 2
    @test ss_code_to_number("H") == 2
    @test ss_code_to_number.(split("H B")) == [2, 7]
    @test ss_number_to_code(2) == "H"
    @test ss_number_to_code(Int32(2)) == "H"
end

end # module ProteinSecondaryStructures

