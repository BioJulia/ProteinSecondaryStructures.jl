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

function dssp_run(pdb_file::String; fix_header=true)
    # if the header is not in the correct format, dssp will fail
    if fix_header
        atoms = PDBTools.readPDB(pdb_file, "protein")
        tmp_file = tempname()*".pdb"
        PDBTools.writePDB(atoms, tmp_file; header=dssp_pdb_header(), footer=nothing)
    else
        tmp_file = pdb_file
    end
    # Run dssp on the pdb file
    dssp_init_line =
    "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA"
    dssp_raw_data = readchomp(pipeline(`$dssp_executable --output-format dssp $tmp_file`))
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
                chain = "$(line[12])",
                residue = residue,
                resnum = parse(Int, line[6:10]),
                sscode = "$(line[17])",
                phi = parse(Float64, line[104:109]),
                psi = parse(Float64, line[110:115]),
                kappa = parse(Float64, line[92:97]),
                alpha = parse(Float64, line[98:103]),
            )
        )
    end
    return ssvector
end

function dssp_run(atoms::AbstractVector{<:PDBTools.Atom})
    tmp_file = tempname()*".pdb"
    PDBTools.writePDB(atoms, tmp_file; header=dssp_pdb_header())
    ss = dssp_run(tmp_file; fix_header=false)
    rm(tmp_file)
    return ss
end
