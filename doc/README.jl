# # TaoNet

# This module currently implements a simulation algorithm for a simple Markov
# model of the evolution of *gene family syntenic networks*. It simulates the
# synteny graph (sensu Zhao & Schranz 2018) in stages across a given a dated
# species tree using a Gillespie-like algorithm. The continuous time Markov
# model includes the following events:

# 1. gene duplication at *per-gene* duplication rate `λ`
# 2. gene loss at *per gene* loss rate `μ`
# 3. rearrangement at *per gene* rearrangement rate `ν`

# Two more parameters are defined in the simple model:

# 1. the probability an edge of a node is retained upon rearrangement `pr`
# 2. the probability an edge of a parent node is copied to a duplicated node `pd`

# In the present implementation, the latter two parameters are not set
# directly, but assigned a Beta prior distribution from which `pr` and `pd` are
# randomly sampled independently for each event.

# ## Usage

# This is a julia package, and was developed using julia 1.5. To install the
# package open a julia REPL (type `julia` in a terminal window), then enter the
# package manager by typing `]` and execute the following lines
#
# ```
# add https://github.com/arzwa/NewickTree.jl#master
# add https://github.com/arzwa/TaoNet.jl#master
# ```

# You might also want to `add Plots`

# ## Example
# Note that `Plots` can take a long time to load...
using NewickTree, TaoNet, Plots, Random, Distributions
Random.seed!(195);

# Get the species tree
t = readnw("(((atr:2.47,(osa:1.82,vvi:1.82):0.65):0.91,(gbi:2.90,"*
    "(gmo:1.77,wmi:1.77):1.13):0.48):0.80,(afi:0.89,scu:0.89):3.31);")

# Now define the model
tn = SimpleTaoNet(λ=0.2, μ=0.15, ν=0.05, pr=Beta(1,5), pd=Beta(1,5))

# and simulate a random syntenic network
G = simulate(tn, t)

# and we can plot it, this might take some time to load
taonetplot(G, t; curves=false)

# Now the interesting stuff will mostly be to simulate genome-scale profiles
Gs = simulate(tn, t, 10)
profile_df = profile(Gs, t)

# This profile is however not directly related to the kind of profiles Tao
# makes. Tao's profiles are based on clustering the MCScanX based network. I
# think the equivalent of Tao's clusters in the simulated networks would be to
# simply take all connected components from the simulated gene family networks.
# This is implemented with the following function
clusters_df = cluster_profile(Gs, t)

# which can be easily 'binarized'
binary_df = clusters_df .> 0

# which can be used for phylogenetic analysis. To do phylogenetic analysis, I
# guess it's convenient to have the matrix in more 'phylogenetic' output
# formats? Not sure, but it's implemented: use `to_fasta("output_file_name",
# binary_df)` or `to_phylip("output_file_name", binary_df)`.

# Note that a condition can be given to `simulate` to only retain simulated graphs
# that satisfy the condition, for instance having more than one edge. The
# default condition is having at least one edge and more than one node, which
# looks like this
using LightGraphs
Gs = simulate(tn, t, 10, cond=(G)->nv(G) > 1 && ne(G) > 0)

using Literate  #src
Literate.markdown(@__FILE__, joinpath(@__DIR__, ".."), execute=true, documenter=false)  #src
