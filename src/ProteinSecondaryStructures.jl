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
export class
export ss_composition
export ss_code_to_number
export ss_number_to_code

export ss_content
export ss_map

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
const classes = Dict{String,String}(
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
const number_to_code = Dict{Int,String}(
    code_to_number[key] => key for key in keys(code_to_number)
)

"""
    ss_code_to_number(code::Union{String,Char})

Converts secondary structure codes to code integer numbers, which 
might be useful for generating secondary structure plots.

The classification follows the DSSP standard classes:

https://m3g.github.io/ProteinSecondaryStructures.jl/stable/explanation/#Secondary-structure-classes

# Example

```jldoctest
julia> using ProteinSecondaryStructures

julia> ss_code_to_number("H")
2

julia> ss_code_to_number('B')
7

julia> class(7)
"beta bridge"
```

!!! compat
    This function was added in version 1.1.0

"""
ss_code_to_number(code::Union{Char,String}) = code_to_number[string(code)]

"""
    ss_number_to_code(code::Integer)

Converts secondary structure integer codes to the corresponding string codes.

The classification follows the DSSP standard classes:

https://m3g.github.io/ProteinSecondaryStructures.jl/stable/explanation/#Secondary-structure-classes

# Example

```jldoctest
julia> using ProteinSecondaryStructures

julia> ss_number_to_code(2)
"H"

julia> ss_number_to_code(7)
"B"
```

!!! compat
    This function was added in version 1.1.0

"""
ss_number_to_code(code::Integer) = number_to_code[code]

"""
    class(ss::Union{SSData, SSData, Integer, String, Char})

Return the secondary structure class. The input may be a `SSData` object, 
a secondary structure `Integer` code (1-8) or a secondary 
structure code (`G, H, ..., C`).

The secondary structure classes are:

| Secondary structure | `ss code`    | `code number`|
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

"""
class(ss::SSData) = classes[ss.sscode]
class(code_number::Int) = classes[number_to_code[code_number]]
class(sscode::Union{Char,String}) = classes[string(sscode)]

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
    for sscode in keys(classes)
        sscomposition[class(sscode)] = count(ss -> ss.sscode == sscode, data)
    end
    return sscomposition
end

# Trajectory analysis
include("./trajectories.jl")

# Testing module
include("../test/Testing.jl")

end # module ProteinSecondaryStructures

