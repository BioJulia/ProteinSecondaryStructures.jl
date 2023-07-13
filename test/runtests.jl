using TestItemRunner
@run_package_tests

@testitem "Aqua.test_all" begin
    import Aqua
    Aqua.test_all(ProteinSecondaryStructures)
end

@testitem "from pdb file" begin
    using ProteinSecondaryStructures.Testing
    for i in eachindex(examples)
        @test ss_composition(stride_run(examples[i].filename)) == examples[i].ss_composition[:stride]
        @test ss_composition(dssp_run(examples[i].filename)) == examples[i].ss_composition[:dssp]
    end
end

@testitem "from Vector{<:PDBTools.Atom}" begin
    using ProteinSecondaryStructures.Testing
    using PDBTools
    for i in eachindex(examples)
        atoms = PDBTools.readPDB(examples[i].filename, "protein")
        @test ss_composition(stride_run(atoms)) == examples[i].ss_composition[:stride]
        @test ss_composition(dssp_run(atoms)) == examples[i].ss_composition[:dssp]
    end
end

