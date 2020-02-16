using OceanDistributions
using Test

@testset "OceanDistributions.jl" begin
    Tcensus=readoceandistribution("../examples/Tcensus.txt")
    tmp=sum(Tcensus)[1]
    @test isapprox(tmp,1.3349978769392906e18; atol=1e16)
end
