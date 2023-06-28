module Stride

using TestItems

import Chemfiles
import PDBTools

export StrideData
export secondary_structure
export is_helix, is_alphahelix, is_pihelix, is_310helix 
export is_strand, is_betastrand, is_betabridge
export is_bend, is_coil, is_turn
export class
export ss_composition
export ss_content
export ss_map

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

@testitem "from pdb file" begin
    using Stride.Testing
    pdbfile = joinpath(Testing.data_dir,"pdb","pdb1fmc.pdb")
    ss = secondary_structure(pdbfile)
    @test length(ss) == 510
    @test ss_composition(ss) == Dict{String, Int}(
        "3₁₀ helix" => 21, 
        "bend" => 0, 
        "turn" => 70, 
        "helix" => 263, 
        "beta strand" => 77, 
        "alpha helix" => 242, 
        "pi helix" => 0, 
        "beta bridge" => 4, 
        "strand" => 81, 
        "coil" => 96
    )
end

function secondary_structure(atoms::AbstractVector{<:PDBTools.Atom})
    tmp_file = tempname()
    PDBTools.writePDB(atoms, tmp_file)
    return secondary_structure(tmp_file)
end

@testitem "from Vector{<:PDBTools.Atom}" begin
    using Stride.Testing
    using PDBTools
    pdbfile = joinpath(Testing.data_dir,"pdb","pdb1fmc.pdb")
    atoms = readPDB(pdbfile, "protein and chain A")
    ss = secondary_structure(atoms)
    @test length(ss) == 255 
    @test ss_composition(ss) == Dict{String, Int}(
        "3₁₀ helix"   => 11,
        "bend"        => 0,
        "turn"        => 34,
        "helix"       => 132,
        "beta strand" => 39,
        "alpha helix" => 121,
        "pi helix"    => 0,
        "beta bridge" => 2,
        "strand"      => 41,
        "coil"        => 48,
    )
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

const ss_code = Dict{String,Int}(
    "G" => 1,
    "H" => 2,
    "I" => 3,
    "T" => 4,
    "E" => 5,
    "B" => 6,
    "S" => 7,
    "C" => 8,
)

"""
    class(ss::StrideData)

Return the secondary structure class of the data. The classes are:

* `"3₁₀ helix"`
* `"α helix"`
* `"π helix"`
* `"turn"`
* `"β strand"`
* `"β bridge"`
* `"bend"`
* `"coil"`

"""
class(ss::StrideData) = classes[ss.sstype]

@doc """
    is_helix(ss::StrideData)
    is_alphahelix(ss::StrideData)
    is_pihelix(ss::StrideData)
    is_310helix(ss::StrideData)
    is_strand(ss::StrideData)
    is_betastrand(ss::StrideData)
    is_betabridge(ss::StrideData)
    is_turn(ss::StrideData)
    is_bend(ss::StrideData)
    is_coil(ss::StrideData)

Return `true` if the data is of the given secondary structure type.

"""
function is_function end

@doc (@doc is_function) is_helix(ss::StrideData) = ss.sstype in ("H", "G", "I")
@doc (@doc is_function) is_alphahelix(ss::StrideData) = ss.sstype == "H"
@doc (@doc is_function) is_pihelix(ss::StrideData) = ss.sstype == "I"
@doc (@doc is_function) is_310helix(ss::StrideData) = ss.sstype == "G"
@doc (@doc is_function) is_strand(ss::StrideData) = ss.sstype in ("E", "B")
@doc (@doc is_function) is_betastrand(ss::StrideData) = ss.sstype == "E"
@doc (@doc is_function) is_betabridge(ss::StrideData) = ss.sstype == "B"
@doc (@doc is_function) is_turn(ss::StrideData) = ss.sstype == "T"
@doc (@doc is_function) is_bend(ss::StrideData) = ss.sstype == "S"
@doc (@doc is_function) is_coil(ss::StrideData) = ss.sstype == "C"

"""
    ss_composition(data::AbstractVector{<:StrideData})

Calculate the secondary structure composition of the data. Returns a dictionary of
the secondary structure types and their counts.

"""
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

function ss_frame!(
    atoms::AbstractVector{<:PDBTools.Atom},
    atom_indices::Vector{Int},
    frame::Chemfiles.Frame,
)
   coordinates = Chemfiles.positions(frame)
   for (i,col) in enumerate(eachcol(coordinates))
       if i in atom_indices
          x, y, z = col[1], col[2], col[3]
          atoms[i].x = x
          atoms[i].y = y
          atoms[i].z = z
       end
   end
   tmp_file = tempname()
   PDBTools.writePDB(atoms, tmp_file)
   ss = secondary_structure(tmp_file)
   rm(tmp_file) 
   return ss
end

"""
    ss_content(f::F, atoms::AbstractVector{<:PDBTools.Atom}, trajectory::Chemfiles.Trajectory)

Calculate the secondary structure content of the trajectory. `f` is the function that returns,
for each residue, if the secondary structure is of a certain type. For example, to calculate 
the alpha helix content, use `f = is_alphahelix`.

"""
function ss_content(
    f::F, 
    atoms::AbstractVector{<:PDBTools.Atom}, 
    trajectory::Chemfiles.Trajectory, 
) where {F<:Function}
    atom_indices = [atom.index for atom in atoms]
    ss_content = zeros(Float64, length(trajectory))
    for (iframe,frame) in enumerate(trajectory)
        ss = ss_frame!(atoms, atom_indices, frame)
        ss_content[iframe] = count(f, ss) / length(ss)
    end
    return ss_content
end

@testitem "ss_content" begin
    using Stride.Testing
    using PDBTools
    using Chemfiles
    pdbfile = joinpath(Testing.data_dir,"Gromacs","system.pdb")
    trajectory = Trajectory(joinpath(Testing.data_dir,"Gromacs","trajectory.xtc"))
    helical_content = ss_content(is_helix, readPDB(pdbfile, "protein"), trajectory)
    @test length(helical_content) == 26
    @test sum(helical_content)/length(helical_content) ≈ 0.2181174089068826
end

"""
    ss_map(atoms::AbstractVector{<:PDBTools.Atom}, trajectory::Chemfiles.Trajectory)

Calculate the secondary structure map of the trajectory. Returns a matrix of the secondary.

"""
function ss_map(
    atoms::AbstractVector{<:PDBTools.Atom},
    trajectory::Chemfiles.Trajectory,
)
    atom_indices = [atom.index for atom in atoms]
    ss_map = zeros(Int,length(PDBTools.eachresidue(atoms)), length(trajectory))
    for (iframe, frame) in enumerate(trajectory)
        ss = ss_frame!(atoms, atom_indices, frame)
        for (i, ssdata) in pairs(ss)
            ss_map[i,iframe] = ss_code[ssdata.sstype]
        end
    end
    return ss_map
end

@testitem "ss_map" begin
    using Stride.Testing
    using PDBTools
    using Chemfiles
    pdbfile = joinpath(Testing.data_dir,"Gromacs","system.pdb")
    trajectory = Trajectory(joinpath(Testing.data_dir,"Gromacs","trajectory.xtc"))
    ssmap = ss_map(readPDB(pdbfile, "protein"), trajectory)
    @test size(ssmap) == (76, 26)
    @test sum(ssmap) == 9032
end

# Testing module
include("../test/Testing.jl")

end # module Stride

