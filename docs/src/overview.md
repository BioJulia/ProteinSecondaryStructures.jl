```@meta
CollapsedDocStrings = true
```

# Overview

This package uses STRIDE or DSSP to compute secondary structures of proteins from their structure files:

- STRIDE supports reading PDB files.
- DSSP supports reading both PDB and mmCIF files.

## Secondary structure classes

The outputs of STRIDE or DSSP follow the convention of secondary strucure codes, which are:

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

## Data structure

The outputs of the secondary structure calculations are vectors of `SSData` elements, with the following data, for each residue:

* `resname::String` - residue name
* `chain::String` - chain identifier
* `resnum::Int` - residue number
* `sscode::String` - secondary structure code
* `phi::Float64` - phi dihedral angle
* `psi::Float64` - psi dihedral angle

For example: 

```jldoctest
julia> using ProteinSecondaryStructures

julia> pdbfile = ProteinSecondaryStructures.Testing.examples[1].filename;

julia> ss = stride_run(pdbfile);

julia> ss[1]
SSData("MET", "A", 1, "C", 360.0, 150.62)

julia> ss[1].sscode
"C"

julia> ss_name(ss[1].sscode)
"coil"
```

### Reference structure

```@docs
SSData
```


