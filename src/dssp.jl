#
# Interface for the DSSP secondary stucture determination program
#

import DSSP_jll
dssp_executable = `$(DSSP_jll.mkdssp()) --mmcif-dictionary $(DSSP_jll.mmcif_pdbx_dic)`

"""
    dssp_run(input_file::String; adjust_pdb=false)

Run DSSP on the PDB or CIF file provided and return a vector containing the detailed
secondary structure information for each residue.

The `adjust_pdb` option is used to fix the header of the pdb file before running `dssp`, 
which is a common problem for computing the secondary structure from PDB files.
In this case, only the ATOM lines are kept in the pdb file. If `adjust_pdb=false`, the pdb file
provided is used as is. 

Note that `DSSP` will fail if residue or atoms types not recognized or if the header
of the PDB file does not follow the necessary pattern.

"""
function dssp_run end

const dssp_header = strip("""
          HEADER
          TITLE     TITLE
          COMPND    MOL_ID: 1;
          SOURCE    MOL_ID: 1;
          KEYWDS    KEYWORDS
          EXPDTA    X-RAY DIFFRACTION
          AUTHOR    ABC
          CRYST1
          """)
const dssp_init_line = "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA"

function parse_dssp_output(dssp_output::AbstractString)
    ss_vector = SSData[]
    struct_summary_section = false
    ifield = 0
    fields = Dict{String,Int}()
    for line in eachline(IOBuffer(dssp_output))
        if startswith(line, "_dssp_struct_summary")
            struct_summary_section = true
            ifield += 1
            field_data = strip.(split(line,"."))
            fields[field_data[end]] = ifield    
        end
        if struct_summary_section
            if !startswith(line, "_dssp_struct_summary")
                if line[1] == '#' 
                    break
                end
                data = split(line)
                chain = data[fields["label_asym_id"]] 
                if chain == "."
                    chain = " "
                end
                sscode=data[fields["secondary_structure"]]
                if sscode == "."
                    sscode = " "
                end
                phi = tryparse(Float64, data[fields["phi"]])
                psi = tryparse(Float64, data[fields["psi"]])
                kappa = tryparse(Float64, data[fields["kappa"]])
                alpha = tryparse(Float64, data[fields["alpha"]])
                ss_residue = SSData(
                    resname=data[fields["label_comp_id"]],
                    chain=data[fields["label_asym_id"]],
                    resnum=parse(Int, data[fields["label_seq_id"]]),
                    sscode=sscode,
                    phi=isnothing(phi) ? 0.0 : phi,
                    psi=isnothing(psi) ? 0.0 : psi,
                    kappa=isnothing(kappa) ? 0.0 : kappa,
                    alpha=isnothing(alpha) ? 0.0 : alpha,
                )
                push!(ss_vector, ss_residue)
            end
        end
    end
    return ss_vector
end

function dssp_run(input_file::AbstractString; adjust_pdb=false)
    # if the header is not in the correct format, dssp will fail
    structure_file = if adjust_pdb
        adjust_pdb_file(input_file; header=dssp_header, empty_chain_identifier="X")
    else
        input_file
    end
    # Run dssp
    dssp_output = try
        readchomp(pipeline(`$dssp_executable --output-format mmcif $structure_file`))
    catch
        "error executing dssp on $input_file"
    end
    ss_vector = parse_dssp_output(dssp_output)
    adjust_pdb && rm(structure_file)
    return ss_vector
end
