module Stride

import PDBTools

export StrideData
export secondary_structure
export is_helix, is_alphahelix, is_pihelix, is_310helix 
export is_strand, is_betastrand, is_betabridge
export is_bend, is_coil, is_turn
export class
export ss_composition

# Assume "stride" is in the path
executable = "stride"

import Base: @kwdef # for 1.6 compatibility
@kwdef struct StrideData
    resname::String
    chain::String
    residue::Int
    resnum::Int
    sstype::String
    phi::Float64
    psi::Float64
    area::Float64
end

"""
    secondary_structure(pdb_file::String)
    secondary_structure(atoms::AbstractVector{<:PDBTools.Atom})

Run stride on the pdb file and return a vector containing the `stride` detailed
secondary structure information for each residue.

The `stride` executable must be in the path.

"""
function secondary_structure end

function secondary_structure(pdb_file::String)
    # Run stride on the pdb file
    stride_raw_data = readchomp(pipeline(`$executable $pdb_file`))
    ssvector = StrideData[]
    for line in split(stride_raw_data, "\n")
        if startswith(line, "ASG")
            residue_data = split(line)
            push!(ssvector,
                StrideData(
                    resname = residue_data[2],
                    chain = residue_data[3],
                    residue = parse(Int, residue_data[4]),
                    resnum = parse(Int, residue_data[5]),
                    sstype = residue_data[6],
                    phi = parse(Float64, residue_data[8]),
                    psi = parse(Float64, residue_data[9]),
                    area = parse(Float64, residue_data[10])
                )
            )
        end
    end
    return ssvector
end

function secondary_structure(atoms::AbstractVector{<:PDBTools.Atom})
    tmp_file = tempname()
    PDBTools.writePDB(atoms, tmp_file)
    return secondary_structure(tmp_file)
end

const classes = Dict{String,String}(
    "G" => "3₁₀ helix",
    "H" => "α helix",
    "I" => "π helix",
    "T" => "turn",
    "E" => "β strand",
    "B" => "β bridge",
    "S" => "bend",
    "C" => "coil",
)
class(ss::StrideData) = classes[ss.sstype]

is_helix(ss::StrideData) = ss.sstype in ("H", "G", "I")
is_alphahelix(ss::StrideData) = ss.sstype == "H"
is_pihelix(ss::StrideData) = ss.sstype == "I"
is_310helix(ss::StrideData) = ss.sstype == "G"

is_strand(ss::StrideData) = ss.sstype in ("E", "B")
is_betastrand(ss::StrideData) = ss.sstype == "E"
is_betabridge(ss::StrideData) = ss.sstype == "B"

is_turn(ss::StrideData) = ss.sstype == "T"
is_bend(ss::StrideData) = ss.sstype == "S"
is_coil(ss::StrideData) = ss.sstype == "C"

function ss_composition(data::AbstractVector{<:StrideData})
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
    return Dict(
        "helix" => helix,
        "alpha helix" => alpha_helix,
        "pi helix" => pihelix,
        "3₁₀ helix" => helix310,
        "strand" => strand,
        "beta strand" => betastrand,
        "beta bridge" => betabridge,
        "turn" => turn,
        "bend" => bend,
        "coil" => coil,
    )
end

end # module Stride
