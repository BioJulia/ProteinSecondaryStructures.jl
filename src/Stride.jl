module Stride

using TestItems
import Chemfiles
import PDBTools

export SSData

export stride_run
export dssp_run 

export is_helix, is_alphahelix, is_pihelix, is_310helix 
export is_strand, is_betastrand, is_betabridge
export is_bend, is_coil, is_turn
export class
export ss_composition
export ss_content
export ss_map

import Base: @kwdef # for 1.6 compatibility
@kwdef struct SSData 
    resname::String
    chain::String
    residue::Int
    resnum::Int
    sscode::String
    phi::Float64
    psi::Float64
    area::Float64 = 0.0 # stride only
    kappa::Float64 = 0.0 # dssp only
    alpha::Float64 = 0.0 # dssp only
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
    "T" => "turn",
    "E" => "beta strand",
    "B" => "beta bridge",
    "S" => "bend",
    "C" => "coil",
    " " => "loop",
)

const ssenum = Dict{String,Int}(
    "G" => 1, # 310 helix
    "H" => 2, # alpha helix
    "I" => 3, # pi helix
    "T" => 4, # turn
    "E" => 5, # beta strand
    "B" => 6, # beta bridge
    "S" => 7, # bend
    "C" => 8, # coil
    " " => 9, # loop
)

const code_to_ss = Dict{Int,String}(
   1 => "G", # 310 helix
   2 => "H", # alpha helix
   3 => "I", # pi helix
   4 => "T", # turn
   5 => "E", # beta strand
   6 => "B", # beta bridge
   7 => "S", # bend
   8 => "C", # coil
   9 => " ", # loop
)

"""
    class(ss::Union{SSData, SSData, Int, String})

Return the secondary structure class. The input may be a `SSData` object, 
a secondary structure `Int` code (1-8) or a secondary structure type string (`G, H, ..., C`).

The secondary structure classes are:

| Secondary structure | `sscode`     | `ssenum`    |
|:--------------------|:------------:|:------------:|
| `"310 helix"`       | `"G"`        | `1`          | 
| `"alpha helix"`     | `"H"`        | `2`          |
| `"pi helix"`        | `"I"`        | `3`          |
| `"turn"`            | `"T"`        | `4`          |
| `"beta strand"`     | `"E"`        | `5`          |
| `"beta bridge"`     | `"B"`        | `6`          |
| `"bend"`            | `"S"`        | `7`          |
| `"coil"`            | `"C"`        | `8`          |
| `"loop"`            | `" "`        | `9`          |

"""
class(ss::SSData) = classes[ss.sscode]
class(ssenum::Int) = classes[code_to_ss[ssenum]]
class(sscode::String) = classes[sscode]

@doc """
    is_helix(ss::SSData)
    is_alphahelix(ss::SSData)
    is_pihelix(ss::SSData)
    is_310helix(ss::SSData)
    is_strand(ss::SSData)
    is_betastrand(ss::SSData)
    is_betabridge(ss::SSData)
    is_turn(ss::SSData)
    is_bend(ss::SSData)
    is_coil(ss::SSData)
    is_turn(ss::SSData)

Return `true` if the data is of the given secondary structure type.

"""
function is_function end

@doc (@doc is_function) is_helix(ss::SSData) = ss.sscode in ("H", "G", "I")
@doc (@doc is_function) is_alphahelix(ss::SSData) = ss.sscode == "H"
@doc (@doc is_function) is_pihelix(ss::SSData) = ss.sscode == "I"
@doc (@doc is_function) is_310helix(ss::SSData) = ss.sscode == "G"
@doc (@doc is_function) is_strand(ss::SSData) = ss.sscode in ("E", "B")
@doc (@doc is_function) is_betastrand(ss::SSData) = ss.sscode == "E"
@doc (@doc is_function) is_betabridge(ss::SSData) = ss.sscode == "B"
@doc (@doc is_function) is_turn(ss::SSData) = ss.sscode == "T"
@doc (@doc is_function) is_bend(ss::SSData) = ss.sscode == "S"
@doc (@doc is_function) is_coil(ss::SSData) = ss.sscode == "C"
@doc (@doc is_function) is_loop(ss::SSData) = ss.sscode == " "

"""
    ss_composition(data::AbstractVector{<:SSData})

Calculate the secondary structure composition of the data. Returns a dictionary of
the secondary structure types and their counts.

"""
function ss_composition(data::AbstractVector{<:SSData})
    helix = count(is_helix, data)
    alpha_helix = count(is_alphahelix, data)
    pihelix = count(is_pihelix, data)
    helix310 = count(is_310helix, data)
    strand = count(is_strand, data)
    betastrand = count(is_betastrand, data)
    betabridge = count(is_betabridge, data)
    turn = count(is_turn, data)
    bend = count(is_bend, data)
    coil = count(is_coil, data)
    loop = count(is_loop, data)
    return Dict(
        "helix" => helix,
        "alpha helix" => alpha_helix,
        "pi helix" => pihelix,
        "310 helix" => helix310,
        "strand" => strand,
        "beta strand" => betastrand,
        "beta bridge" => betabridge,
        "turn" => turn,
        "bend" => bend,
        "coil" => coil,
        "loop" => loop,
    )
end

# Trajectory analysis
include("./trajectories.jl")

# Testing module
include("../test/Testing.jl")

end # module Stride

