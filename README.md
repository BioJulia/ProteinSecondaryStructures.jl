# Stride.jl

This package parses the [Stride](https://en.wikipedia.org/wiki/STRIDE_(algorithm)) and [DSSP](https://github.com/PDB-REDO/dssp) secondary structure prediction program outputs, to make them convenient to use from Julia, particularly for the analysis of MD simulations. 

1. [Installation](#installation)
2. [Usage](#usage)

   - [Computing from a PDB file](#computing-the-secondary-structure-from-a-pdb-file)
   - [Computing from a PDBTools selection](#computing-from-a-selection-of-atoms-using-pdbtools)
   - [Interpretation of the output](#interpretation-of-the-output)
   - [Computing the secondary structure map along the trajectory](#computing-the-secondary-structure-map)
   - [Computing the secondary structure in a MD trajectory](#computing-the-ammount-of-secondary-structure-classes-in-a-trajectory)


3. [References](#references)

## Installation

Obtain `stride` and/or `dssp` programs from the corresponding repositories, compile it and add it to your path:

For `stride`:

```bash
git clone https://github.com/MDAnalysis/stride
cd stride/src
make
cp ./stride /usr/local/bin # or somewhere in the path
```

For `dssp`:

```
sudo apt install dssp # on Ubuntu-based linux distributions
```

or go to: https://github.com/PDB-REDO/dssp

In Julia, install `Stride` with:

```julia
import Pkg; Pkg.add(url="https://github.com/m3g/Stride.jl")
```

## Usage

### Computing the secondary structure from a PDB file

```julia
julia> using Stride

julia> pdbfile = Stride.Testing.examples[1].filename
"/home/user/.julia/dev/Stride/test/data/pdb/pdb1fmc.pdb"

julia> ss = stride_run(pdbfile)
510-element Vector{SSData}:
 SSData("MET", "A", 1, "C", 360.0, 150.62, 234.4, 0.0, 0.0)
 SSData("PHE", "A", 2, "C", -69.01, 138.78, 162.9, 0.0, 0.0)
 ⋮
 SSData("ASN", "B", 255, "C", -130.75, 360.0, 114.8, 0.0, 0.0)
```

or use `dssp_run` to use the `dssp` method.

### Computing from a selection of atoms using PDBTools

Install `PDBTools` with:
```julia
import Pkg; Pkg.add("PDBTools")
```

Then:

```julia
julia> using Stride, PDBTools

julia> pdbfile = Stride.Testing.examples[1].filename
"/home/user/.julia/dev/Stride/test/data/pdb/pdb1fmc.pdb"

julia> pdb = readPDB(pdbfile, "protein and chain A") # read only protein atoms from chain A
   Array{Atoms,1} with 1876 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     MET     A        1        1   22.518   17.379   31.003  1.00 34.99     1       -         1
       2   CA     MET     A        1        1   23.426   17.764   32.113  1.00 34.03     1       -         2
                                                       ⋮ 
    1876  OXT     ASN     A      255      255   44.876   38.278   -2.437  1.00 54.24     1       -      1876

julia> ss = stride_run(pdb)
255-element Vector{SSData}:
 SSData("MET", "A", 1, "C", 360.0, 150.62, 234.4, 0.0, 0.0)
 SSData("PHE", "A", 2, "C", -69.01, 138.78, 162.9, 0.0, 0.0)
 ⋮
 SSData("ASN", "A", 255, "C", -130.97, 360.0, 100.9, 0.0, 0.0)
```

## Interpretation of the output

The output is a vector of `SSData` elements, with the following data, for each residue:

```julia
struct SSData
    resname::String
    chain::String
    resnum::Int
    sscode::String
    phi::Float64
    psi::Float64
    area::Float64 = 0.0 # stride only
    kappa::Float64 = 0.0 # dssp only
    alpha::Float64 = 0.0 # dssp only
end
```

### Classes of secondary structure

The classes of secondary structure and their codes are:

| Secondary structure | `ss code`    |`code number` |
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

### Secondary structure composition

It is possible to extract the composition of the secondary structure with:

```julia
julia> c = ss_composition(ss)
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

julia> c["alpha helix"]
242
```

The class of secondary structure of each residue can be retrived with the `class` function:

```julia
julia> ss[10]
SSData("ASP", "A", 10, "T", -53.61, 124.03, 78.7, 0.0, 0.0)

julia> class(ss[10])
"turn"
```

And helper functions are available to obtain boolean vectors to check if each residue belong to some class. For example:

```julia
julia> is_helix(ss[10])
false

julia> is_helix.(ss)
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
is_helix
is_alphahelix
is_pihelix
is_310helix
is_kappahelix
is_strand
is_betastrand
is_betabridge
is_turn
is_bend
is_coil
```

### Computing the secondary structure map

```julia
using Stride
import PDBTools
import Chemfiles

pdbfile = Stride.Testing.data_dir*"/Gromacs/system.pdb"
trajectory_file = Stride.Testing.data_dir*"/Gromacs/trajectory.xtc"

atoms = PDBTools.readPDB(pdbfile, "protein")
trajectory = Chemfiles.Trajectory(trajectory_file)

ssmap = ss_map(atoms, trajectory) # returns a matrix
```

Similarly, the method used can be selected with:
```julia
ssmap = ss_map(atoms, trajectory; method=stride_run) # returns a matrix
ssmap = ss_map(atoms, trajectory; method=dssp_run) # returns a matrix
```

This will create a matrix that can be visualized, for instance, with:

```julia
using Plots
heatmap(ssmap,xlabel="step",ylabel="residue")
```

producing the figure:

![heatmap](./test/map.png)

where the colors refer to the `code number` fields of the [classes table](#classes-of-secondary-structure) above.

### Computing the ammount of secondary structure classes in a trajectory 

#### From the secondary structure map

Calling `ss_content` with a class identifier function and a map (as computed above), will return the content
of that class along the trajectory:

```julia
julia> ss_content(is_alphahelix, ssmap)
26-element Vector{Float64}:
 0.21052631578947367
 0.15789473684210525
 ⋮
 0.13157894736842105
```

The composition of classes for a given frame can also be retrieved from the content map:

```julia
julia> ss_composition(ssmap, 6)
Dict{String, Int64} with 10 entries:
  "310 helix"   => 7
  "bend"        => 0
  "turn"        => 17
  "kappa helix" => 0
  "beta strand" => 25
  "beta bridge" => 2
  "alpha helix" => 12
  "pi helix"    => 0
  "loop"        => 0
  "coil"        => 13
```

These functions are useful, because the computation of the secondary structure along the
trajectory (the map) can be costly.

#### For a single class, along the trajectory

```julia
using Stride
import PDBTools
import Chemfiles

pdbfile = Stride.Testing.data_dir*"/Gromacs/system.pdb"
trajectory_file = Stride.Testing.data_dir*"/Gromacs/trajectory.xtc"

atoms = PDBTools.readPDB(pdbfile, "protein")
trajectory = Chemfiles.Trajectory(trajectory_file)

helical_content = ss_content(is_helix, atoms, trajectory)
```

Optionally, the method to compute the secondary structure can be defined, with, for example: 

```julia
helical_content = ss_content(is_helix, atoms, trajectory; method=stride_run)
#or
helical_content = ss_content(is_helix, atoms, trajectory; method=dssp_run)
```

## Note

This is a package under early development. New functionality and breaking changes may occur. 

## References

If you use the `stride` algorithm for secondary structure prediction, please cite:

- Frishman,D & Argos,P. (1995) Knowledge-based secondary structure assignment. Proteins: structure, function and genetics, 23, 566-579.
- Kabsch,W. & Sander,C. (1983) Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers, 22: 2577-2637.

If you use the `dssp` algorithm for secondary structure prediction, please cite:

- Joosten RP, te Beek TAH, Krieger E, Hekkelman ML, Hooft RWW, Schneider R, Sander C, Vriend A series of PDB related databases for everyday needs. Nuc. Acids Res. 2010; 39:D411-D419.
- Kabsch W, Sander C. Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers 1983; 22:2577-2637. 

