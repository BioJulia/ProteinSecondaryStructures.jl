#
# Interface for the STRIDE secondary structure prediction program
#

# Assume "stride" is in the path
stride_executable = "stride"

"""
    stride_run(pdb_file::String)
    stride_run(atoms::AbstractVector{<:PDBTools.Atom})

Run stride on the pdb file and return a vector containing the `stride` detailed
secondary structure information for each residue.

The `stride` executable must be in the path or, alternatively, the `Stride.stride_executable`
variable can be set to the full path of the executable.

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
    # If the header is not in the correct format, stride will fail
    if fix_header
        tmp_file = tempname()*".pdb"
        PDBTools.writePDB(atoms, tmp_file; header=stride_pdb_header())
    else
        tmp_file = pdb_file
    end
    # Run stride on the pdb file
    stride_raw_data = try 
        readchomp(pipeline(`$stride_executable $tmp_file`))
    catch
    end
    ssvector = [ SSData(r.resname, r.chain, r.resnum) for r in PDBTools.eachresidue(atoms) ]
    for line in split(stride_raw_data, "\n")
        if startswith(line, "ASG")
            residue_data = split(line)
            ss_residue = SSData(
                resname = residue_data[2],
                chain = residue_data[3],
                resnum = parse(Int, residue_data[4]),
                sscode = residue_data[6],
                phi = parse(Float64, residue_data[8]),
                psi = parse(Float64, residue_data[9]),
                area = parse(Float64, residue_data[10])
            )
            iss = findfirst(ss -> residue_match(ss_residue, ss), ssvector)
            if !isnothing(iss)
                ssvector[iss] = ss_residue
            end
        end
    end
    return ssvector
end

# From a PDB file
function stride_run(pdb_file::String; fix_header=true) 
    atoms = PDBTools.readPDB(pdb_file, "protein")
    return stride_run(atoms; fix_header=fix_header)
end

