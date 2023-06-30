#
# Interface for the DSSP secondary stucture determination program
#

# Assume "dssp" is in the path
dssp_executable = "dssp"

"""
    dssp_run(pdb_file::String)
    dssp_run(atoms::AbstractVector{<:PDBTools.Atom})

Run DSSP on the pdb file and return a vector containing the detailed
secondary structure information for each residue.

The `dssp` executable must be in the path or, alternatively, the `Stride.dssp_executable`
variable can be set to the full path of the executable.

"""
function dssp_run end

const dssp_init_line =
    "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA"

function dssp_run(pdb_file::String)
    # Run dssp on the pdb file
    dssp_raw_data = readchomp(pipeline(`$dssp_executable --output-format dssp $pdb_file`))
    ssvector = SSData[]
    data_begin = false
    residue = 0
    for line in split(dssp_raw_data, "\n")
        if line == dssp_init_line
            data_begin = true
            continue
        end
        !data_begin && continue
        line[14] == '!' && continue
        residue += 1
        push!(ssvector,
            SSData(
                resname = PDBTools.threeletter(line[14]),
                chain = "$line[12]",
                residue = residue,
                resnum = parse(Int, line[6:10]),
                sscode = "$line[17]",
                phi = parse(Float64, line[104:109]),
                psi = parse(Float64, line[110:115]),
                kappa = parse(Float64, line[92:97]),
                alpha = parse(Float64, line[98:103]),
            )
        )
    end
    return ssvector
end

@testitem "DSSP: from pdb file" begin
    using Stride.Testing
    pdbfile = joinpath(Testing.data_dir,"pdb","pdb1fmc.pdb")
    ss = dssp_run(pdbfile)
    @test length(ss) == 510
    @test ss_composition(ss) == Dict{String, Int}(
        "310 helix" => 21, 
        "bend" => 0, 
        "turn" => 70, 
        "helix" => 263, 
        "beta strand" => 77, 
        "alpha helix" => 242, 
        "pi helix" => 0, 
        "beta bridge" => 4, 
        "strand" => 81, 
        "coil" => 96,
        "loop" => 96,
    )
end

function dssp_run(atoms::AbstractVector{<:PDBTools.Atom})
    tmp_file = tempname()*".pdb"
    @show tmp_file
    PDBTools.writePDB(atoms, tmp_file)
    ss = dssp_run(tmp_file)
    rm(tmp_file)
    return ss
end

@testitem "DSSP: from Vector{<:PDBTools.Atom}" begin
    using Stride.Testing
    using PDBTools
    pdbfile = joinpath(Testing.data_dir,"pdb","pdb1fmc.pdb")
    atoms = readPDB(pdbfile, "protein and chain A")
    ss = dssp_run(atoms)
    @test length(ss) == 255 
    @test ss_composition(ss) == Dict{String, Int}(
        "310 helix"   => 11,
        "bend"        => 0,
        "turn"        => 34,
        "helix"       => 132,
        "beta strand" => 39,
        "alpha helix" => 121,
        "pi helix"    => 0,
        "beta bridge" => 2,
        "strand"      => 41,
        "coil"        => 48,
        "loop"        => 48,
    )
end