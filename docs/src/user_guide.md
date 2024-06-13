```@meta
CollapsedDocStrings = true
```

## Compute secondary structures

First, install the package, according to the [Installation](@ref) instructions.

Load the package with:
```julia-repl
julia> using ProteinSecondaryStructures
```

Here, we illustrate the computation of the secondary structure of a PDB file provides
as an example:
```julia-repl
julia> pdbfile = ProteinSecondaryStructures.Testing.examples[1].filename
"/home/user/.julia/dev/ProteinSecondaryStructures/test/data/pdb/pdb1fmc.pdb"
```

Then, to compute the secondary structure using the `STRIDE` algorithm, use:
```julia-repl
julia> ss = stride_run(pdbfile)
510-element Vector{SSData}:
 SSData("MET", "A", 1, "C", 360.0, 150.62, 234.4, 0.0, 0.0)
 SSData("PHE", "A", 2, "C", -69.01, 138.78, 162.9, 0.0, 0.0)
 ⋮
 SSData("ASN", "B", 255, "C", -130.75, 360.0, 114.8, 0.0, 0.0)
```

The output is a vector of `SSData` elements, which contain the residue name, 
chain, residue number, and secondary structure code. The list of codes follow
the DSSP convention, described in [Secondary structure classes](@ref).

!!! note 
    Alternativelly to the `stride_run` function, the `dssp_run` function can be used, to compute the secondary structures as defined by DSSP.

The details of the `SSData` structure that contain the output for each residue are described in the [Data structure](@ref) section.

## Secondary structure composition

Given the `ss` output of the `stride_run` or `dssp_run` functions, an overview of content of the secondary structure can be obtained with the
`ss_composition` function:

```julia-repl
julia> comp = ss_composition(ss)
Dict{String, Int64} with 12 entries:
  "bend"        => 0
  "kappa helix" => 0
  "beta strand" => 77
  "strand"      => 81
  "loop"        => 0
  "310 helix"   => 21
  "turn"        => 70
  "helix"       => 263
  "beta bridge" => 4
  "alpha helix" => 242
  "pi helix"    => 0
  "coil"        => 96

julia> comp["alpha helix"]
242
```

The output is a dictionary containing the number of residues that were classified in each class. As shown above, this number can be retrieved individually.

## Retrieving names, codes, and numeric codes

The name, single-character codes, or numeric codes of the secondary structure of each residue can be retrieved with the `ss_name`, `ss_code`,
and `ss_number` functions. The input of these functions can be an instance of `SSData` or one of the other two secondary structure
classification types (name, code, or number): 

```@docs
ss_name
ss_code
ss_number
```

These functions can be used to obtain arrays of codes, by broadcasting over the
vector of secondary structure data. For example:

```jldoctest
julia> using ProteinSecondaryStructures

julia> using ProteinSecondaryStructures.Testing: examples

julia> ss = stride_run(examples[1].filename);

julia> ss_name.(ss)[1:5]
5-element Vector{String}:
 "coil"
 "coil"
 "coil"
 "310 helix"
 "310 helix"

```



