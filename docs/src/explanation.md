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

| SS name             | `ss code`    |`code number` |
|---------------------|:------------:|:------------:|
| `"310 helix"`       | `"G"`        | `1`          | 
| `"alpha helix"`     | `"H"`        | `2`          |
| `"pi helix"`        | `"I"`        | `3`          |
| `"kappa helix"`     | `"P"`        | `4`          |
| `"turn"`            | `"T"`        | `5`          |
| `"beta strand"`     | `"E"`        | `6`          |
| `"beta bridge"`     | `"B"`        | `7`          |
| `"bend"`            | `"S"`        | `8`          |
| `"coil"`            | `"C"`        | `9`          |
| `"loop"`            | `" "`        | `10`         |

See the [DSSP secondary structure classification](https://pdb-redo.eu/dssp/about) for further information.

The `code number` are used here to output the matrices of secondary structures obtained
from trajectory runs. 

## Data structure

The output of the secondary structure calculations is a vector of `SSData` elements, with the following data, for each residue:

* `resname::String` - residue name
* `chain::String` - chain identifier
* `resnum::Int` - residue number
* `sscode::String` - secondary structure code
* `phi::Float64` - phi dihedral angle
* `psi::Float64` - psi dihedral angle
* `area::Float64` - solvent accessible area (stride specific)
* `kappa::Float64` - virtual bond angle (dssp specific)
* `alpha::Float64` - virtual torsion angle (dssp specific)

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

julia> ss_name(ss[1].sscode)
"coil"
```



