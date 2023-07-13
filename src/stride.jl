#
# Interface for the STRIDE secondary structure prediction program
#

import STRIDE_jll
stride_executable = STRIDE_jll.stride_exe()

"""
    stride_run(pdb_file::String; selection="protein")
    stride_run(atoms::AbstractVector{<:PDBTools.Atom})

Run stride on the pdb file and return a vector containing the `stride` detailed
secondary structure information for each residue.

The `stride` executable must be in the path or, alternatively, the `ProteinSecondaryStructures.stride_executable`
variable can be set to the full path of the executable.

When passing a PDB file, by default only the atoms belonging to standard protein residues are considered. This can be changed by
setting the `selection` keyword argument to a different selection string, or function, following the `PDBTools.jl` syntax.
Note that `STRIDE` will fail if residue or atoms types not recognized. 

"""
function stride_run end

function stride_pdb_header()
    strip(
        """  
        HEADER    ABC  PROTEIN                               01-JAN-00   1ABC
        CRYST1
        """)
end

function stride_run(atoms::AbstractVector{<:PDBTools.Atom}; fix_header=true)
    tmp_file = tempname() * ".pdb"
    # If the header is not in the correct format, stride will fail
    if fix_header
        PDBTools.writePDB(atoms, tmp_file; header=stride_pdb_header(), footer=nothing)
    else
        PDBTools.writePDB(atoms, tmp_file)
    end
    # Run stride on the pdb file
    stride_raw_data = try
        readchomp(pipeline(`$stride_executable $tmp_file`))
    catch
        "error"
    end
    ssvector = [SSData(r.resname, r.chain, r.resnum) for r in PDBTools.eachresidue(atoms)]
    #      ASG    Detailed secondary structure assignment
    #      Format:  6-8  Residue name
    #	      10-10 Protein chain identifier
    #	      12-15 PDB	residue	number
    #	      17-20 Ordinal residue number
    #	      25-25 One	letter secondary structure code	**)
    #	      27-39 Full secondary structure name
    #	      43-49 Phi	angle
    #	      53-59 Psi	angle
    #	      65-69 Residue solvent accessible area
    for line in split(stride_raw_data, "\n")
        if startswith(line, "ASG")
            ss_residue = SSData(
                resname=line[6:8],
                chain=line[10:10],
                resnum=parse(Int, line[12:15]),
                sscode=uppercase(line[25:25]), # stride sometimes returns "b" instead of "B"
                phi=parse(Float64, line[43:49]),
                psi=parse(Float64, line[53:59]),
                area=parse(Float64, line[65:69])
            )
            iss = findfirst(ss -> residue_match(ss_residue, ss), ssvector)
            if !isnothing(iss)
                ssvector[iss] = ss_residue
            end
        end
    end
    rm(tmp_file)
    return ssvector
end

# From a PDB file
function stride_run(pdb_file::String; selection="protein", fix_header=true)
    atoms = PDBTools.readPDB(pdb_file, selection)
    return stride_run(atoms; fix_header=fix_header)
end

