module TaoNet

using NewickTree, Parameters, Distributions
using LightGraphs, MetaGraphs, GraphRecipes
using StatsBase, ColorSchemes, DataFrames
export SimpleTaoNet, cluster_profile, profile, taonetplot

# A simpler implementation of the TaoNet, without branch-specific rates etc.
# (more involved but outdated version in `TaoNets.jl`)
randexp(λ) = -log(rand())/λ

iswgd(n) = startswith(name(n), "wgd")
iswgt(n) = startswith(name(n), "wgt")
wgdid(n) = parse(Int, split(name(n), "_")[2])

struct ShiftedGeometric{T} <: DiscreteUnivariateDistribution
    p::T
end
Base.rand(d::ShiftedGeometric) = rand(Geometric(d.p)) + 1

@with_kw struct SimpleTaoNet{T,V<:DiscreteUnivariateDistribution}
    λ::T       # duplication rate
    μ::T       # loss rate
    ν::T       # rearrangement rate
    q::Vector{T} = Float64[]  # wgd retention rates
    pr::Beta{Float64} = Beta(1,3) # edge retention probability rearrangement
    pd::Beta{Float64} = Beta(1,3) # edge retention probability upon duplication
    root::V = ShiftedGeometric(0.66)
end

getevents(tn::SimpleTaoNet) = [dup!, loss!, rearrange!]
randrootstate(tn::SimpleTaoNet) = rand(tn.root)
totalrate(tn::SimpleTaoNet, G) = sum(nv(G) .* [tn.λ, tn.μ, tn.ν])

function randevent!(G, tn::SimpleTaoNet)
    ev! = sample(getevents(tn), Weights([tn.λ, tn.μ, tn.ν]))
    ev!(G, tn)
end

function dup!(G, tn::SimpleTaoNet)
    @unpack pd = tn
    parent = rand(1:nv(G))
    add_vertex!(G, :sp, get_prop(G, parent, :sp))
    lv = nv(G)
    r = rand(pd)
    for n in neighbors(G, parent)
        rand() < r ? add_edge!(G, lv, n) : continue
    end
end

function loss!(G, tn::SimpleTaoNet)
    v = rand(1:nv(G))
    rem_vertex!(G, v)
end

function rearrange!(G, tn::SimpleTaoNet)
    @unpack pr = tn
    n = rand(1:nv(G))
    nn = length(neighbors(G, n))
    nlost = nn - floor(Int64, rand(Binomial(nn, rand(pr))))
    if nlost != 0
        for v in rand(neighbors(G, n), nlost)
            rem_edge!(G, n, v)
        end
    end
end

function Base.rand(tn::SimpleTaoNet, tree, N::Integer;
        condition=(G)->nv(G) > 1 && ne(G) > 0)
    G = rand(tn, tree)
    Gs = condition(G) ? typeof(G)[G] : typeof(G)[]
    while length(Gs) < N
        G = rand(tn, tree)
        condition(G) && push!(Gs, G)
    end
    return Gs
end

function Base.rand(tn::SimpleTaoNet, tree)
    # initialize the graph
    X₀ = randrootstate(tn)
    G = MetaGraph(complete_graph(X₀))
    root = getroot(tree)
    for i=1:X₀
        set_prop!(G, i, :sp, id(root))
    end

    # iterate over 'time slices'
    t = 0.0
    h = labeled_history(tree)
    for (tend, node) in h
        t = sim_edges!(G, t, tend, tn)
        sim_node!(G, node, tn)
    end
    sim_edges!(G, t, treeheight(tree), tn)
    return G
end

function sim_edges!(G, t, tend, tn::SimpleTaoNet)
    t += randexp(totalrate(tn, G))
    while t < tend  # simulate
        randevent!(G, tn)
        t += randexp(totalrate(tn, G))
    end
    return tend
end

function sim_node!(G, node, tn::SimpleTaoNet)
    iswgd(node) && return sim_wgd!(G, node, tn)
    iswgt(node) && return sim_wgt!(G, node, tn)
    for n in filter_vertices(G, :sp, id(node))  # speciation
        set_prop!(G, n, :sp, id(node[1]))
        add_vertex!(G, :sp, id(node[2]))
        lv = nv(G)
        for nn in neighbors(G, n)
            add_edge!(G, nn, lv)
        end
        add_edge!(G, n, lv)
    end
end

function sim_wgd!(G, node, tn::SimpleTaoNet)
    q = tn.q[wgdid(node)]
    for n in filter_vertices(G, :sp, id(node))
        set_prop!(G, n, :sp, id(node[1]))
        if rand() < q
            add_vertex!(G, :sp, id(node[1]))
            lv = nv(G)
            for nn in neighbors(G, n)
                add_edge!(G, nn, lv)
            end
            add_edge!(G, n, lv)
        end
    end
end

function sim_wgt!(G, node, tn::SimpleTaoNet)
    q = tn.q[wgdid(node)]
    for n in filter_vertices(G, :sp, id(node))
        set_prop!(G, n, :sp, id(node[1]))
        retained = sample(1:3, Weights([(1.0-q)^2, 2q*(1.0-q), q^2]))
        for i=1:retained-1
            add_vertex!(G, :sp, x.i)
            lv = nv(G)
            for nn in neighbors(G, n)
                add_edge!(G, nn, lv)
            end
            add_edge!(G, n, lv)
        end
    end
end

height(n) = sum(distance.(NewickTree.getpath(n)[1:end-1]))
treeheight(n) = height(getleaves(n)[1])

function labeled_history(tree)
    h = [(height(n), n) for n in prewalk(tree) if !isleaf(n)]
    sort!(h, by=x->x[1])
    return h
end

function profile(Gs::AbstractVector, tree)
    nd = Dict(name(n)=>id(n) for n in getleaves(tree))
    DataFrame(profile.(Gs, Ref(nd)))
end

function profile(G::MetaGraph)
    cmap = countmap([v[:sp] for (k,v) in G.vprops])
end

function profile(G::MetaGraph, ndict::Dict)
    p = profile(G)
    (; [Symbol(k)=>haskey(p, v) ? p[v] : 0 for (k,v) in ndict]...)
end

function cluster_profile(Gs::AbstractVector, tree)
    nd = Dict(name(n)=>id(n) for n in getleaves(tree))
    DataFrame(vcat(cluster_profile.(Gs, Ref(nd))...))
end

function cluster_profile(G::MetaGraph)
    cc = connected_components(G)
    [countmap([get_prop(G, v, :sp) for v in c]) for c in cc if length(c) > 1]
end

function cluster_profile(G::MetaGraph, ndict::Dict)
    cs = cluster_profile(G)
    [(; [Symbol(k)=>haskey(c, v) ? c[v] : 0 for (k,v) in ndict]...) for c in cs]
end

function taonetplot(G, tree; cs = ColorSchemes.Spectral_11,
        linecolor=:darkgrey, linealpha=0.5, dim=2,
        markersize=0.2, nodeshape=:circle, kwargs...)
    nd = Dict(id(n)=>name(n) for n in getleaves(tree))
    nl = length(nd)
    cs = Dict(id(n)=>get(cs, i/nl) for (i,n) in enumerate(getleaves(tree)))
    mc = [cs[get_prop(G, i, :sp)] for i=1:nv(G)]
    ns = [nd[get_prop(G, i, :sp)] for i=1:nv(G)]
    graphplot(G, markercolor=mc, names=ns, linecolor=linecolor,
        linealpha=linealpha, dim=dim, nodeshape=nodeshape;
        markersize=markersize, kwargs...)
end

end # module
