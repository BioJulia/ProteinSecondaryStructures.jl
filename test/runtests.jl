using TestItemRunner: @run_package_tests
using TestItems: @testitem
@run_package_tests

@testitem "Aqua.test_all" begin
    import Aqua
    Aqua.test_all(ProteinSecondaryStructures)
end

@testitem "from pdb file" begin
    using ProteinSecondaryStructures.Testing: examples
    for i in eachindex(examples)
        @test ss_composition(stride_run(examples[i].filename)) == examples[i].ss_composition[:stride]
        @test ss_composition(dssp_run(examples[i].filename)) == examples[i].ss_composition[:dssp]
    end
end

@testitem "from Vector{<:PDBTools.Atom}" begin
    using ProteinSecondaryStructures.Testing: examples
    import PDBTools
    for i in eachindex(examples)
        atoms = PDBTools.readPDB(examples[i].filename, "protein")
        @test ss_composition(stride_run(atoms)) == examples[i].ss_composition[:stride]
        @test ss_composition(dssp_run(atoms)) == examples[i].ss_composition[:dssp]
    end
end

@testiem "ss_code_to_number/ss_number_to_code" begin
    @test ss_code_to_number('H') == 2
    @test ss_code_to_number("H") == 2
    @test ss_code_to_number.(split("H B")) == [2, 7]
    @test ss_number_to_code(2) == "H"
    @test ss_number_to_code(Int32(2)) == "H"
end

