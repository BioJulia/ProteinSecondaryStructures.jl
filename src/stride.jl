#
# Interface for the STRIDE secondary structure prediction program
#

import STRIDE_jll
stride_executable = STRIDE_jll.stride_exe()

"""
    stride_run(pdb_file::AbstractString; adjust_pdb=false)

Run stride on the PDB file and return a vector containing the `stride` detailed
secondary structure information for each residue.

- `adjust_pdb=true` is used to adjust the format of the PDB file before running `stride`,
   which requires a specific header and a specific empty-chain identifier. 
   In this case, only the ATOM lines are kept in the PDB file. 
- If `adjust_pdb=false`, the PDB file provided is used as is.

Note that `STRIDE` will fail if residue or atoms types not recognized or if the header
of the PDB file does not follow the necessary pattern.

!!! note 
    - STRIDE might ignore some residues in the PDB file if they are not recognized, or incomplete. 
    - STRIDE does not support structures in mmCIF format.

"""
function stride_run end

const stride_header = strip(
        """  
        HEADER    ABC  PROTEIN                               01-JAN-00   1ABC
        CRYST1
        """)

function parse_stride_output(stride_output::AbstractString)
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
    ss_vector = SSData[]
    for line in eachline(IOBuffer(stride_output))
        if startswith(line, "ASG")
            chain = line[10]
            chain = chain == '-' ? ' ' : chain
            ss_residue = SSData(
                line[6:8], # resname
                string(chain), # chain
                parse(Int, line[12:15]), # resnum
                string(uppercase(line[25])), # sscode: stride sometimes returns "b" instead of "B"
                parse(Float64, line[43:49]), # phi
                parse(Float64, line[53:59]), # psi
            )
            push!(ss_vector, ss_residue)
        end
    end
    return ss_vector
end

function stride_run(input_pdb_file::AbstractString; adjust_pdb=false)
    # If the header is not in the correct format, stride will fail
    pdb_file = if adjust_pdb
        adjust_pdb_file(input_pdb_file; header=stride_header, empty_chain_identifier="-")
    else
        input_pdb_file
    end
    stride_output = try
        readchomp(pipeline(`$stride_executable $pdb_file`))
    catch
        "error running stride on $pdb_file"
    end
    ss_vector = parse_stride_output(stride_output)
    adjust_pdb && rm(pdb_file)
    return ss_vector
end

@testitem "STRIDE" begin
    using ProteinSecondaryStructures.Testing: examples
    # PDB files
    for i in 1:6
        @test ss_composition(stride_run(examples[i].filename; adjust_pdb=true)) == examples[i].ss_composition[:stride]
    end
end