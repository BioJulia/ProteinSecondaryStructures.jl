#
# Interface for the DSSP secondary stucture determination program
#

import DSSP_jll
dssp_executable = `$(DSSP_jll.mkdssp()) --mmcif-dictionary $(DSSP_jll.mmcif_pdbx_dic)`

three_letter_code = Dict{Char,String}(
    'A' => "ALA",
    'R' => "ARG",
    'N' => "ASN",
    'D' => "ASP",
    'C' => "CYS",
    'Q' => "GLN",
    'E' => "GLU",
    'G' => "GLY",
    'H' => "HIS",
    'I' => "ILE",
    'L' => "LEU",
    'K' => "LYS",
    'M' => "MET",
    'F' => "PHE",
    'P' => "PRO",
    'S' => "SER",
    'T' => "THR",
    'W' => "TRP",
    'Y' => "TYR",
    'V' => "VAL",
)

"""
    dssp_run(pdb_file::String; adjust_pdb=true)

Run DSSP on the pdb file and return a vector containing the detailed
secondary structure information for each residue.

The `adjust_pdb` option is used to fix the header of the pdb file before running `dssp`.
In this case, only the ATOM lines are kept in the pdb file. If `adjust_pdb=false`, the pdb file
provided is used as is.

Note that `DSSP` will fail if residue or atoms types not recognized or if the header
of the PDB file does not follow the necessary pattern.

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

function dssp_run(pdbfile::AbstractString; adjust_pdb=true)
    # if the header is not in the correct format, dssp will fail
    pdbfile = if adjust_pdb
        tmp_file = tempname() * ".pdb"
        open(tmp_file, "w") do io
            println(io, dssp_pdb_header())
            for line in readlines(pdbfile)
                if startswith(line, "ATOM") 
                    # DSSP does not support empty chain identifiers
                    if line[22] == ' '
                        line = line[1:21]*'X'*line[23:end]
                    end
                    println(io, line)
                end
            end
        end
        tmp_file
    else
        pdbfile
    end
    ssvector = init_ssvector(pdbfile)
    # Run dssp on the pdb file
    dssp_raw_data = try
        readchomp(pipeline(`$dssp_executable --output-format dssp $pdbfile`))
    catch
        "error executing dssp on $pdbfile"
    end
    isnothing(dssp_raw_data) && return ssvector
    data_begin = false
    for line in eachline(IOBuffer(dssp_raw_data))
        if line == dssp_init_line
            data_begin = true
            continue
        end
        !data_begin && continue
        line[14] == '!' && continue
        ss_residue = SSData(
            resname=three_letter_code[line[14]],
            chain=line[12:12],
            resnum=parse(Int, line[6:10]),
            sscode=line[17:17],
            phi=parse(Float64, line[104:109]),
            psi=parse(Float64, line[110:115]),
            kappa=parse(Float64, line[92:97]),
            alpha=parse(Float64, line[98:103]),
        )
        iss = findfirst(ss -> residues_match(ss_residue, ss), ssvector)
        if !isnothing(iss)
            ssvector[iss] = ss_residue
        end
    end
    return ssvector
end
