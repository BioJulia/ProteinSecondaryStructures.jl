module Testing
    using ..Stride
    using Chemfiles
    using PDBTools

    export data_dir

    data_dir = joinpath(@__DIR__, "data")

    #pdbfile = joinpath(data_dir, "Gromacs/system.pdb")
    #trajectory = Trajectory(joinpath(data_dir, "Gromacs/trajectory.xtc"))

end