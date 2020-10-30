using Pkg; Pkg.activate("/home/arzwa/dev/TaoNet")
using NewickTree, TaoNet, CSV, Distributions, StatsBase, UUIDs, DataFrames

# configuration ----------------------------------------------------------------
treefile   = joinpath(@__DIR__, "data/monocots3.nw")  # reference species tree
outputfile = joinpath(@__DIR__, "results.csv")  # output file path
tmpdir     = "/tmp/taonet-test"; mkpath(tmpdir)  # working directory
nrep       = 250  # nr. of replicates to simulate
nfamilies  = 500  # nr. of families to simulate per replicate (!= # clusters)

# Now I implement the simulations slightly different, the configuration is 
# basically the following function, which generates an array of `N` model 
# objects, representing a model for each family.
function randmodel(N)
    # fixed parameters
    λmean = 0.38
    λstd1 = 0.5  # σ of log-normal distribution of λ values across *simulations*
    λstd2 = 0.5  # σ of log-normal distribution of λ values across *families*
    pbeta = Beta(2,2)  
    root  = TaoNet.ShiftedGeometric(0.66)

    # draw mean λ and mean ν from same distribution
    λ, ν = rand(LogNormal(log(λmean), λstd1), 2)

    # draw family specific rates
    λs = rand(LogNormal(log(λ), λstd2), N)
    νs = rand(LogNormal(log(ν), λstd2), N)
    
    # draw β parameters
    b  = rand(pbeta)
    β  = Beta(b, b)

    # generate N models with family-specific rates
    M = map(i->SimpleTaoNet(λ=λs[i], μ=λs[i], ν=νs[i], pr=β, pd=β, root=root), 1:N)
    θ = (λ=λ, ν=ν, b=b)
    M, θ
end
# -----------------------------------------------------------------------------

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
function dosims(speciestree, modelfun, nfamilies, nrep)
    tmptrefile = joinpath(tmpdir, "sptree.nw")
    writenw(tmptrefile, unroot(speciestree))
    map(1:nrep) do i
        # draw a random parameter set (model parameterization)
        model, θ = modelfun(nfamilies)
        @info "simulation" i θ

        # simulate the synteny network and obtain the matrix
        graphs   = simulate(model, speciestree)
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
        (θ=θ, clusters=clusters, res=(rf=rf, tree=nwtree))
    end
end

sptre = readnw(readline(treefile))
out   = dosims(sptre, randmodel, nfamilies, nrep)

df = hcat(DataFrame(first.(out)), DataFrame(last.(out)))
CSV.write(joinpath(@__DIR__, outputfile), df)

# Example of the kind of model we're simulating from
M, x = randmodel(100);

p1 = plot(LogNormal(log(0.38), 0.5), grid=false, label="prior")
vline!([x.λ], label="average λ across families")
vline!([x.ν], label="average μ across families")
p2 = scatter([tn.λ for tn in M], [tn.ν for tn in M])
p3 = histogram([tn.λ for tn in M], bins=50)
p4 = histogram([tn.ν for tn in M], bins=50)
plot(p1, p2, p3, p4, grid=false)
