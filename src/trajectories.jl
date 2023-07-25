#
# MD Trajectory analysis
#

function ss_frame!(
    atoms::AbstractVector{<:PDBTools.Atom},
    atom_indices::Vector{Int},
    frame::Chemfiles.Frame;
    method::F=stride_run
) where {F<:Function}
    coordinates = Chemfiles.positions(frame)
    for (i, col) in enumerate(eachcol(coordinates))
        iatom = findfirst(==(i), atom_indices)
        if !isnothing(iatom)
            atoms[iatom].x = col[1]
            atoms[iatom].y = col[2]
            atoms[iatom].z = col[3]
        end
    end
    return method(atoms)
end

"""
    ss_map(pdbfile::AbstractString, trajectory::Chemfiles.Trajectory; selection="protein")
    ss_map(atoms::AbstractVector{<:PDBTools.Atom}, trajectory::Chemfiles.Trajectory)

Calculate the secondary structure map of the trajectory. 
Returns a matrix of secondary structure codes, where each row is a residue and each column is a frame.

The atoms to be considered may be provided directly with an `atoms` vector, or by reading a PDB file,
in which case the `selection` keyword argument is used to select the atoms to be considered.

The `show_progress` keyword argument controls whether a progress bar is shown.

"""
function ss_map(
    atoms::AbstractVector{<:PDBTools.Atom},
    trajectory::Chemfiles.Trajectory;
    method::F=stride_run,
    show_progress=true
) where {F<:Function}
    atom_indices = [atom.index for atom in atoms]
    ss_map = zeros(Int, length(PDBTools.eachresidue(atoms)), length(trajectory))
    p = Progress(length(trajectory); enabled=show_progress)
    for (iframe, frame) in enumerate(trajectory)
        ss = ss_frame!(atoms, atom_indices, frame; method=method)
        for (i, ssdata) in pairs(ss)
            ss_map[i, iframe] = code_to_number[ssdata.sscode]
        end
        next!(p)
    end
    return ss_map
end

function ss_map(
    pdbfile::AbstractString,
    trajectory::Chemfiles.Trajectory;
    method::F=stride_run,
    selection::AbstractString="protein"
) where {F<:Function}
    atoms = PDBTools.readPDB(pdbfile, selection)
    return ss_map(atoms, trajectory; method=method)
end

@testitem "ss_map" begin
    using ProteinSecondaryStructures.Testing
    import PDBTools
    import Chemfiles
    pdbfile = joinpath(Testing.data_dir, "Gromacs", "system.pdb")
    trajectory = Chemfiles.Trajectory(joinpath(Testing.data_dir, "Gromacs", "trajectory.xtc"))
    # With stride
    ssmap = ss_map(PDBTools.readPDB(pdbfile, "protein"), trajectory; method=stride_run)
    @test size(ssmap) == (76, 26)
    @test sum(ssmap) == 10577
    # With DSSP
    ssmap = ss_map(PDBTools.readPDB(pdbfile, "protein"), trajectory; method=dssp_run)
    @test size(ssmap) == (76, 26)
    @test sum(ssmap) == 12069
    # From the pdb file name
    ssmap = ss_map(pdbfile, trajectory; method=stride_run, selection="protein")
    @test size(ssmap) == (76, 26)
    @test sum(ssmap) == 10577
end

"""
    ss_content(f::F, atoms::AbstractVector{<:PDBTools.Atom}, trajectory::Chemfiles.Trajectory)

Calculate the secondary structure content of the trajectory. `f` is the function that returns,
for each residue, if the secondary structure is of a certain type. For example, to calculate 
the alpha helix content, use `f = is_alphahelix`.

The `show_progress` keyword argument controls whether a progress bar is shown.

"""
function ss_content(
    f::F,
    atoms::AbstractVector{<:PDBTools.Atom},
    trajectory::Chemfiles.Trajectory;
    method::G=stride_run,
    show_progress=true
) where {F<:Function,G<:Function}
    atom_indices = [atom.index for atom in atoms]
    ss_content = zeros(Float64, length(trajectory))
    p = Progress(length(trajectory); enabled=show_progress)
    for (iframe, frame) in enumerate(trajectory)
        ss = ss_frame!(atoms, atom_indices, frame; method=method)
        ss_content[iframe] = count(f, ss) / max(1, length(ss))
        next!(p)
    end
    return ss_content
end

is_anyhelix(sscode::Int) = number_to_code[sscode] in ("H", "G", "I", "P")
is_alphahelix(sscode::Int) = number_to_code[sscode] == "H"
is_pihelix(sscode::Int) = number_to_code[sscode] == "I"
is_kappahelix(sscode::Int) = number_to_code[sscode] == "P"
is_310helix(sscode::Int) = number_to_code[sscode] == "G"
is_anystrand(sscode::Int) = number_to_code[sscode] in ("E", "B")
is_betastrand(sscode::Int) = number_to_code[sscode] == "E"
is_betabridge(sscode::Int) = number_to_code[sscode] == "B"
is_turn(sscode::Int) = number_to_code[sscode] == "T"
is_bend(sscode::Int) = number_to_code[sscode] == "S"
is_coil(sscode::Int) = number_to_code[sscode] == "C"
is_loop(sscode::Int) = number_to_code[sscode] == " "

"""
    ss_content(f::F, ssmap::AbstractMatrix{Int})

Calculates the secondary structure content of the trajectory. `f` is the function that returns,
for each residue, if the secondary structure is of a certain type. For example, to calculate
the alpha helix content, use `f = is_alphahelix`.

Here, `ssmap` is the secondary structure map of the trajectory, as returned by the `ss_map` function.

"""
function ss_content(
    f::F,
    ssmap::AbstractMatrix{Int},
) where {F<:Function}
    ss_content = zeros(Float64, size(ssmap, 2))
    for (iframe, sscomposition) in enumerate(eachcol(ssmap))
        ss_content[iframe] = count(f, sscomposition) / max(1, length(sscomposition))
    end
    return ss_content
end

"""
    ss_composition(ssmap::AbstractMatrix{Int}, iframe::Int}

Calculates the secondary structure composition of the trajectory. Returns a dictionary of
the secondary structure types and their counts, for the chosen frame.

"""
function ss_composition(ssmap::AbstractMatrix{Int}, iframe::Int)
    sscomposition = Dict{String,Int}()
    sscodes = @view(ssmap[:, iframe])
    for sscode in keys(number_to_code)
        sscomposition[class(sscode)] = count(==(sscode), sscodes)
    end
    return sscomposition
end

@testitem "ss_content/ss_composition" begin
    using ProteinSecondaryStructures.Testing
    using PDBTools: readPDB
    using Chemfiles: Trajectory
    pdbfile = joinpath(Testing.data_dir, "Gromacs", "system.pdb")
    trajectory = Trajectory(joinpath(Testing.data_dir, "Gromacs", "trajectory.xtc"))
    # With stride
    helical_content = ss_content(is_anyhelix, readPDB(pdbfile, "protein"), trajectory; method=stride_run)
    @test length(helical_content) == 26
    @test sum(helical_content) / length(helical_content) ≈ 0.2181174089068826
    # With DSSP
    helical_content = ss_content(is_anyhelix, readPDB(pdbfile, "protein"), trajectory; method=dssp_run)
    @test length(helical_content) == 26
    @test sum(helical_content) / length(helical_content) ≈ 0.21204453441295545
    # With non-contiguous indexing
    atoms = readPDB(pdbfile, only=at -> (10 <= at.residue < 30) | (40 <= at.residue < 60))
    helical_content = ss_content(is_anyhelix, atoms, trajectory; method=stride_run)
    @test length(helical_content) == 26
    @test sum(helical_content) / length(helical_content) ≈ 0.20288461538461539
    # From the map
    ssmap = ss_map(readPDB(pdbfile, "protein"), trajectory; method=stride_run)
    helical_content = ss_content(is_anyhelix, ssmap)
    @test length(helical_content) == 26
    @test sum(helical_content) / length(helical_content) ≈ 0.2181174089068826
    # Test ss_composition function
    @test ss_composition(ssmap, 5) == Dict("310 helix" => 6, "bend" => 0, "turn" => 17, "kappa helix" => 0, "beta strand" => 25, "beta bridge" => 2, "alpha helix" => 12, "pi helix" => 0, "loop" => 0, "coil" => 14)
end
