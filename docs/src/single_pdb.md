# Single PDB files

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

And helper functions are available to obtain boolean vectors to check if each residue belong to some class. For example:

```julia
julia> is_anyhelix(ss[10])
false

julia> is_anyhelix.(ss)
255-element BitVector:
 0
 0
 0
 1
 ⋮
 0
 0
 0
```

where in the last example we used the broadcast operation (the dot) to obtain the result of the application of the function for the complete array of residues.

The following functions are available:
```julia
is_anyhelix
is_alphahelix
is_pihelix
is_310helix
is_kappahelix
is_anystrand
is_betastrand
is_betabridge
is_turn
is_bend
is_coil
```

## Using PDBTools

The [PDBTools](http://m3g.github.io/PDBTools.jl) provides a simple
and lighweight interface to select specific subsets of atoms from the
structures.

Install `PDBTools` with:
```julia-repl
julia> import Pkg; Pkg.add("PDBTools")
```

Then, load both packages and set the PDB file name with, for example:

```julia-repl
julia> using ProteinSecondaryStructures, PDBTools

julia> pdbfile = ProteinSecondaryStructures.Testing.examples[1].filename
"/home/user/.julia/dev/ProteinSecondaryStructures/test/data/pdb/pdb1fmc.pdb"
```

The `PDBTtools.readPDB` function can now read tha atom data from the PDB file, but it is possible to select a subset of atoms from the structure with a simple selection syntax:
```julia-repl
julia> pdb = readPDB(pdbfile, "protein and chain A") # read only protein atoms from chain A
   Array{Atoms,1} with 1876 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     MET     A        1        1   22.518   17.379   31.003  1.00 34.99     1       -         1
       2   CA     MET     A        1        1   23.426   17.764   32.113  1.00 34.03     1       -         2
                                                       ⋮ 
    1876  OXT     ASN     A      255      255   44.876   38.278   -2.437  1.00 54.24     1       -      1876
```

The atoms are loaded in a vector of type `Vector{<:PDBTools.Atom}`. This vector can be fed in place of the PDB file names, and the secondary structures will be computed for the selected subsets of atoms:

```julia-repl
julia> ss = stride_run(pdb)
255-element Vector{SSData}:
 SSData("MET", "A", 1, "C", 360.0, 150.62, 234.4, 0.0, 0.0)
 SSData("PHE", "A", 2, "C", -69.01, 138.78, 162.9, 0.0, 0.0)
 ⋮
 SSData("ASN", "A", 255, "C", -130.97, 360.0, 100.9, 0.0, 0.0)
```




