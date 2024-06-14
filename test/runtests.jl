using TestItemRunner: @run_package_tests
using TestItems: @testitem
@run_package_tests

@testitem "Aqua.test_all" begin
    import Aqua
    Aqua.test_all(ProteinSecondaryStructures)
end

@testitem "secondary structure" begin
    using ProteinSecondaryStructures.Testing: examples
    for i in 1:6
        @test ss_composition(stride_run(examples[i].filename; adjust_pdb=true)) == examples[i].ss_composition[:stride]
        @test ss_composition(dssp_run(examples[i].filename; adjust_pdb=true)) == examples[i].ss_composition[:dssp]
    end
    @test ss_composition(dssp_run(examples[7].filename)) == examples[7].ss_composition[:dssp]
end
