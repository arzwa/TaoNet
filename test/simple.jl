using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using Test, TaoNet, NewickTree, Distributions, LightGraphs

t = readnw("(((atr:2.47,(osa:1.82,vvi:1.82):0.65):0.91,(gbi:2.90,"*
    "(gmo:1.77,wmi:1.77):1.13):0.48):0.80,(afi:0.89,scu:0.89):3.31);")
n = length(getleaves(t))

@testset "SimpleTaoNet" begin
    @testset "No events will happen" begin
        tn = SimpleTaoNet(λ=0., μ=0., ν=0., root=DiscreteUniform(1,1))
        G = rand(tn, t)
        @test nv(G) == n
        @test ne(G) == n*(n-1)÷2
    end

    @testset "μ = 0, one component ~ species tree" begin
        tn = SimpleTaoNet(λ=0.4, μ=0., ν=0., root=DiscreteUniform(1,1))
        for i=1:10
            G = rand(tn, t)
            cp = cluster_profile(G)
            @test length(cp[1]) == n
        end
    end
end
