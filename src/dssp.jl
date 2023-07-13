#
# Interface for the DSSP secondary stucture determination program
#

import DSSP_jll
dssp_executable = `$(DSSP_jll.mkdssp()) --mmcif-dictionary $(joinpath(DSSP_jll.artifact_dir, "share", "dssp", "mmcif_pdbx.dic"))`

"""
    dssp_run(pdb_file::String; selection="protein")
    dssp_run(atoms::AbstractVector{<:PDBTools.Atom})

Run DSSP on the pdb file and return a vector containing the detailed
secondary structure information for each residue.

The `dssp` executable must be in the path or, alternatively, the `ProteinSecondaryStructures.dssp_executable`
variable can be set to the full path of the executable.

When passing a PDB file, by default only the atoms belonging to standard protein residues are considered. This can be changed by
setting the `selection` keyword argument to a different selection string, or function, following the `PDBTools.jl` syntax.
Note that `DSSP` will fail if residue or atoms types not recognized. 

"""
function dssp_run end

function dssp_pdb_header()
    strip("""
          HEADER
          TITLE     TITLE
          COMPND    MOL_ID: 1;
          SOURCE    MOL_ID: 1;
          KEYWDS    KEYWORDS
          EXPDTA    X-RAY DIFFRACTION
          AUTHOR    ABC
          CRYST1
          """)
end

const dssp_init_line = "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA"

# From the PDBTools.Atom vector
function dssp_run(atoms::AbstractVector{<:PDBTools.Atom}; fix_header=true)
    tmp_file = tempname() * ".pdb"
    # if the header is not in the correct format, dssp will fail
    if fix_header
        PDBTools.writePDB(atoms, tmp_file; header=dssp_pdb_header(), footer=nothing)
    else
        PDBTools.writePDB(atoms, tmp_file)
    end
    # Run dssp on the pdb file
    dssp_raw_data = try
        readchomp(pipeline(`$dssp_executable --output-format dssp $tmp_file`))
    catch
        "error"
    end
    ssvector = [SSData(r.resname, r.chain, r.resnum) for r in PDBTools.eachresidue(atoms)]
    isnothing(dssp_raw_data) && return ssvector
    data_begin = false
    for line in split(dssp_raw_data, "\n")
        if line == dssp_init_line
            data_begin = true
            continue
        end
        !data_begin && continue
        line[14] == '!' && continue
        ss_residue = SSData(
            resname=PDBTools.threeletter(line[14]),
            chain="$(line[12])",
            resnum=parse(Int, line[6:10]),
            sscode="$(line[17])",
            phi=parse(Float64, line[104:109]),
            psi=parse(Float64, line[110:115]),
            kappa=parse(Float64, line[92:97]),
            alpha=parse(Float64, line[98:103]),
        )
        iss = findfirst(ss -> residue_match(ss_residue, ss), ssvector)
        if !isnothing(iss)
            ssvector[iss] = ss_residue
        end
    end
    rm(tmp_file)
    return ssvector
end

# From the PDB file
function dssp_run(pdb_file::String; selection="protein", fix_header=true)
    atoms = PDBTools.readPDB(pdb_file, selection)
    return dssp_run(atoms; fix_header=fix_header)
end




