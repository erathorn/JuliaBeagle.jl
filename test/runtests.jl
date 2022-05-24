using JuliaBeagle
using Test

@testset "JuliaBeagle.jl" begin
    # Write your tests here.
    tree = ParseNewick(
        "(((0:0.110833,1:0.0137979)10:0.146124,(2:0.197891,(3:0.132967,(4:0.0378759,5:0.089252)11:0.101833)12:0.184301)
         13:0.0450774)14:0.335725,6:0.153197,(7:0.0216218,(8:0.0781687,9:0.120419)15:0.0209114)16:0.0209771);",
    )
    ntax, nchar, gap, miss, symbols, df, langs =
        MCPhylo.ParseNexus("simudata.nex")
    df = MCPhylo.datafortree(df, langs, tree, symbols, gap, miss, log_space=false)
    pden = ones(4) / 4
    pd = PhyloDist(tree, pden, [1.0], [1.0], JC)
end
