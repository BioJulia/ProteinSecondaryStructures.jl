# Molecular Dynamics Trajectories

## Secondary structure map

```julia
using ProteinSecondaryStructures 
using PDBTools: readPDB 
using Chemfiles: Trajectory

pdbfile = ProteinSecondaryStructures.Testing.data_dir*"/Gromacs/system.pdb"
trajectory_file = ProteinSecondaryStructures.Testing.data_dir*"/Gromacs/trajectory.xtc"

atoms = readPDB(pdbfile, "protein")
trajectory = Trajectory(trajectory_file)

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
heatmap(ssmap,
  xlabel="frame",
  ylabel="residue",
  framestyle=:box,
  color=palette(:tab20c,10)
)
```

producing the figure:

![heatmap](./assets/map.svg)

where the colors refer to the `code number` fields of the [classes table](#classes-of-secondary-structure) above.


## Saving and loading a map

The secondary structure map computed is just a matrix of integer codes. Thus, it can be saved or read in any preferred format.
As a suggestion, it is possible to use `writedlm` and `readdlm` function from the `DelimitedFiles` package: 

```julia
using DelimitedFiles
# save data to ssmap.dat
writedlm("ssmap.dat", ssmap)
# load data
ssmat = readdlm("ssmap.dat", Int)
```

## Secondary structure classes

### From the secondary structure map

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

### For a single class, along the trajectory

```julia
using ProteinSecondaryStructures 
using PDBTools: readPDB 
using Chemfiles: Trajectory

pdbfile = ProteinSecondaryStructures.Testing.data_dir*"/Gromacs/system.pdb"
trajectory_file = ProteinSecondaryStructures.Testing.data_dir*"/Gromacs/trajectory.xtc"

atoms = readPDB(pdbfile, "protein")
trajectory = Trajectory(trajectory_file)

helical_content = ss_content(is_anyhelix, atoms, trajectory)
```

Optionally, the method to compute the secondary structure can be defined, with, for example: 

```julia
helical_content = ss_content(is_anyhelix, atoms, trajectory; method=stride_run)
#or
helical_content = ss_content(is_anyhelix, atoms, trajectory; method=dssp_run)
```
