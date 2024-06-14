#
# Interface for the STRIDE secondary structure prediction program
#

import STRIDE_jll
stride_executable = STRIDE_jll.stride_exe()

"""
    stride_run(pdb_file::AbstractString; adjust_pdb=true)

Run stride on the pdb file and return a vector containing the `stride` detailed
secondary structure information for each residue.

The `adjust_pdb` option is used to fix the header of the pdb file before running `stride`.
In this case, only the ATOM lines are kept in the pdb file. If `adjust_pdb=false`, the pdb file
provided is used as is.

Note that `STRIDE` will fail if residue or atoms types not recognized or if the header
of the PDB file does not follow the necessary pattern.

"""
function stride_run end

const stride_header = strip(
        """  
        HEADER    ABC  PROTEIN                               01-JAN-00   1ABC
        CRYST1
        """)

function stride_run(input_pdb_file::AbstractString; adjust_pdb=true)
    # If the header is not in the correct format, stride will fail
    pdb_file = if adjust_pdb
        adjust_pdb_file(input_pdb_file; header=stride_header, empty_chain_identifier="-")
    else
        input_pdb_file
    end
    ssvector = init_ssvector(pdb_file; empty_chain_identifier="-")
    # Run stride on the pdb file
    stride_raw_data = try
        readchomp(pipeline(`$stride_executable $pdb_file`))
    catch
        "error running stride on $pdb_file"
    end
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
    for line in eachline(IOBuffer(stride_raw_data))
        if startswith(line, "ASG")
            chain = line[10]
            if chain == '-'
                chain = ' '
            end
            ss_residue = SSData(
                resname=line[6:8],
                chain=string(chain),
                resnum=parse(Int, line[12:15]),
                sscode=string(uppercase(line[25])), # stride sometimes returns "b" instead of "B"
                phi=parse(Float64, line[43:49]),
                psi=parse(Float64, line[53:59]),
                area=parse(Float64, line[65:69])
            )
            iss = findfirst(ss -> residues_match(ss_residue, ss), ssvector)
            if !isnothing(iss)
                ssvector[iss] = ss_residue
            end
        end
    end
    adjust_pdb && rm(pdb_file)
    return ssvector
end
