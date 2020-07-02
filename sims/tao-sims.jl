using Pkg; Pkg.activate("/home/arzwa/dev/TaoNet")
using NewickTree, TaoNet, CSV, Distributions, StatsBase, UUIDs, DataFrames

# configuration ----------------------------------------------------------------
treefile   = joinpath(@__DIR__, "data/monocots2.nw")  # reference species tree
outputfile = joinpath(@__DIR__, "results.csv")  # output file path
tmpdir     = "/tmp/taonet-test"; mkpath(tmpdir)  # working directory
nrep       = 10    # nr. of replicates to simulate
nfamilies  = 1000  # nr. of families to simulate per replicate (!= # clusters)
plambda    = Uniform(0., 3.)   # distribution from which to sample λ
pnu        = Uniform(0., 3.)   # distribution from which to sample ν
pbeta      = Uniform(1, 10)    # distribution from which to sample β
peta       = Beta(6, 3)      # distribution from which to sample η
# NOTE: if you want to set a constant value for λ, μ, ... just use a Normal
# distribution with zero variance, e.g. if you want to fix λ at 0.1 use
# `pλ = Normal(0.1, 0.)`
# NOTE: I don't think it's worthwhile ding bootstrap analysis/LRT on the
# simulated data
# NOTE: `@__DIR__` expands to a string with the path to the directory where
# the file in which `@__DIR__` is called (i.e. this file) is located.
# ------------------------------------------------------------------------------

# some helper setup, we use a four-parameter model (i.e. a critical BDP)
struct SimulationPriors
    λ; ν; β; η
end

function Base.rand(p::SimulationPriors)
    β = Beta(1, rand(p.β))
    root = TaoNet.ShiftedGeometric(rand(p.η))
    λ = rand(p.λ)
    SimpleTaoNet(λ=λ, μ=λ, ν=rand(p.ν), pr=β, pd=β, root=root)
end

StatsBase.params(p::SimpleTaoNet) = (λ=p.λ, μ=p.μ, ν=p.ν, αr=p.pr.α,
    βr=p.pr.β, αd=p.pd.α, βd=p.pd.β, η=p.root.p, q=p.q)

const iqtfiles = ["ckp.gz", "log", "mldist", "treefile", "iqtree",
    "parstree", "bionj", "treefile.rfdist", "treefile.log"]

function unroot(tree)
    t = deepcopy(tree)
    bifchild = first(filter(!isleaf, children(t)))
    push!(t, bifchild[1])
    push!(t, bifchild[2])
    bifchild[1].parent = t
    bifchild[2].parent = t
    delete!(t, bifchild)
    return t
end

# function for doing the simulations
function dosims(speciestree, prior, nfamilies, nrep)
    tmptrefile = joinpath(tmpdir, "sptree.nw")
    writenw(tmptrefile, unroot(speciestree))
    map(1:nrep) do i
        # draw a random parameter set (model parameterization)
        model = rand(prior)
        @info i params(model)

        # simulate the synteny network and obtain the matrix
        graphs   = rand(model, speciestree, nfamilies)
        clusters = cluster_profile(graphs, speciestree)
        binarydf = clusters .> 0
        tmpfile  = joinpath(tmpdir, string(uuid4()) * ".phy")

        to_phylip(tmpfile, binarydf)
        cmd = `iqtree -s $tmpfile -m MK+R+FO -st MORPH`
        run(cmd)
        nwtree = readline("$tmpfile.treefile")
        run(`iqtree -t $tmpfile.treefile -rf $tmptrefile`)
        rfdist = "$tmpfile.treefile.rfdist"
        rf = parse(Int, split(readlines(rfdist)[2])[2])
        @info "Robinson-Foulds distance = $rf"
        rm.(["$tmpfile.$suffix" for suffix in iqtfiles])
        rm(tmpfile)
        (θ=params(model), clusters=clusters, res=(rf=rf, tree=nwtree))
    end
end

sptre = readnw(readline(treefile))
prior = SimulationPriors(plambda, pnu, pbeta, peta)
out   = dosims(sptre, prior, nfamilies, nrep)

df = hcat(DataFrame(first.(out)), DataFrame(last.(out)))
CSV.write(joinpath(@__DIR__, outputfile), df)

# ps = map([:λ, :ν, :η]) do x
#     scatter(df[!,x], df[!,:rf], ylabel="RF distance", xlabel=string(x),
#         color=:black, markersize=3, alpha=0.5, yticks=0:2:12)
# end
# plot(ps..., grid=false, legend=false, size=(600,250))
# savefig(joinpath(@__DIR__, "results/test.pdf"))
