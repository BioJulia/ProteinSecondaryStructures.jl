#
# MD Trajectory analysis
#

function ss_frame!(
    atoms::AbstractVector{<:PDBTools.Atom},
    atom_indices::Vector{Int},
    frame::Chemfiles.Frame;
    method::F=stride_run,
) where {F<:Function}
   coordinates = Chemfiles.positions(frame)
   for (i,col) in enumerate(eachcol(coordinates))
       if i in atom_indices
          atoms[i].x = col[1]
          atoms[i].y = col[2]
          atoms[i].z = col[3]
       end
   end
   return method(atoms)
end


"""
    ss_content(f::F, atoms::AbstractVector{<:PDBTools.Atom}, trajectory::Chemfiles.Trajectory)

Calculate the secondary structure content of the trajectory. `f` is the function that returns,
for each residue, if the secondary structure is of a certain type. For example, to calculate 
the alpha helix content, use `f = is_alphahelix`.

"""
function ss_content(
    f::F, 
    atoms::AbstractVector{<:PDBTools.Atom}, 
    trajectory::Chemfiles.Trajectory;
    method::G=stride_run,
) where {F<:Function, G<:Function}
    atom_indices = [atom.index for atom in atoms]
    ss_content = zeros(Float64, length(trajectory))
    for (iframe,frame) in enumerate(trajectory)
        ss = ss_frame!(atoms, atom_indices, frame; method=method)
        ss_content[iframe] = count(f, ss) / max(1,length(ss))
    end
    return ss_content
end

@testitem "ss_content" begin
    using Stride.Testing
    using PDBTools
    using Chemfiles
    pdbfile = joinpath(Testing.data_dir,"Gromacs","system.pdb")
    trajectory = Trajectory(joinpath(Testing.data_dir,"Gromacs","trajectory.xtc"))
    # With stride
    helical_content = ss_content(is_helix, readPDB(pdbfile, "protein"), trajectory; method=stride_run)
    @test length(helical_content) == 26
    @test sum(helical_content)/length(helical_content) ≈ 0.2181174089068826
    # With DSSP
    helical_content = ss_content(is_helix, readPDB(pdbfile, "protein"), trajectory; method=dssp_run)
    @test length(helical_content) == 26
    @test sum(helical_content)/length(helical_content) ≈ 0.2181174089068826
end

"""
    ss_map(atoms::AbstractVector{<:PDBTools.Atom}, trajectory::Chemfiles.Trajectory)

Calculate the secondary structure map of the trajectory. Returns a matrix of the secondary.

"""
function ss_map(
    atoms::AbstractVector{<:PDBTools.Atom},
    trajectory::Chemfiles.Trajectory;
    method::F=stride_run,
) where {F<:Function}
    atom_indices = [atom.index for atom in atoms]
    ss_map = zeros(Int,length(PDBTools.eachresidue(atoms)), length(trajectory))
    for (iframe, frame) in enumerate(trajectory)
        ss = ss_frame!(atoms, atom_indices, frame; method=method)
        for (i, ssdata) in pairs(ss)
            ss_map[i,iframe] = ssenum[ssdata.sscode]
        end
    end
    return ss_map
end

@testitem "ss_map" begin
    using Stride.Testing
    using PDBTools
    using Chemfiles
    pdbfile = joinpath(Testing.data_dir,"Gromacs","system.pdb")
    trajectory = Trajectory(joinpath(Testing.data_dir,"Gromacs","trajectory.xtc"))
    # With stride
    ssmap = ss_map(readPDB(pdbfile, "protein"), trajectory; method=stride_run)
    @test size(ssmap) == (76, 26)
    @test sum(ssmap) == 9032
    # With DSSP
    ssmap = ss_map(readPDB(pdbfile, "protein"), trajectory; method=dssp_run)
    @test size(ssmap) == (76, 26)
    @test sum(ssmap) == 9032
end