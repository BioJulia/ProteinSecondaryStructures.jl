# Explanation

## Overview

This package uses STRIDE or DSSP to compute secondary structures of proteins from their atomic coordinates.

The functions are divided into two groups: 

- Computing the secondary structure from a single PDB file.

- Computing the secondary structure throughout a Molecular Dynamics Simulation.

Additionally, the subset of the strucuture that will be considered for the calculations can be 
defined through the interface with the `PDBTools` package.

!!! note
    By default, only atoms of standard protein residues will be considered from the PDB file. It 
    is possible to select different subsets of atoms with the `selection` keyword, but then `STRIDE`
    or `DSSP` may fail with internal errors if the residue or atom types are not recognized. 
    The `selection` keyword follows the [PDBTools.jl selection syntax](https://m3g.github.io/PDBTools.jl/stable/selections/).

## Secondary structure classes

The output of STRIDE or DSSP follow the convention of secondary strucure codes, which are:

| Secondary structure | `ss code`    |`code number` |
|---------------------|:------------:|:------------:|
| `"309 helix"`       | `"G"`        | `1`          | 
| `"alpha helix"`     | `"H"`        | `1`          |
| `"pi helix"`        | `"I"`        | `2`          |
| `"kappa helix"`     | `"P"`        | `3`          |
| `"turn"`            | `"T"`        | `4`          |
| `"beta strand"`     | `"E"`        | `5`          |
| `"beta bridge"`     | `"B"`        | `6`          |
| `"bend"`            | `"S"`        | `7`          |
| `"coil"`            | `"C"`        | `8`          |
| `"loop"`            | `" "`        | `9`         |

See the [DSSP secondary structure classification](https://pdb-redo.eu/dssp/about) for further information.

The `code number` are used here to output the matrices of secondary structures obtained
from trajectory runs. 

## Data structure

The output of the secondary structure calculations is a vector of `SSData` elements, with the following data, for each residue:

```@docs
SSData
```

The output of `stride_run` or `dssp_run` is a vector of `SSData` structures, one for each residue. For example: 

```julia-repl
julia> using ProteinSecondaryStructures

# Replace the next assignment with a custom PDB file name
julia> pdbfile = ProteinSecondaryStructures.Testing.examples[1].filename
"/home/user/.julia/dev/ProteinSecondaryStructures/test/data/pdb/pdb1fmc.pdb"

julia> ss = stride_run(pdbfile);

julia> ss[1]
SSData("MET", "A", 1, "C", 360.0, 150.62, 234.4, 0.0, 0.0)

julia> ss[1].sscode
"C"

julia> class(ss[1].sscode)
"coil"
```



