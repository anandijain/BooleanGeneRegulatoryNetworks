using GeneRegulatoryNetworks, Graphs, CSV, DataFrames
using SimpleWeightedGraphs
using LinearAlgebra, Arpack, Plots
using GLMakie, GraphMakie
using SymbolicUtils: Sym
using SymbolicUtils

to_df(fn; kwargs...) = CSV.read(fn, DataFrame; kwargs...)
df = to_df("data/HumanNet-FN.tsv"; header=false)
c1, c2 = df.Column1, df.Column2
x = sort(unique(c1))
y = sort(unique(c2))
v = union(x, y)
vid = eachindex(v)
d = Dict(vid .=> v)
d2 = Dict(v .=> vid)
e = c1 .=> c2
ed = Dict(e)

# es = map(x->getindex(d2, x)=>, c1)
# g = SimpleGraph(length(vid))
g = SimpleDiGraph(length(vid))
# todo check that there isn't an edge from a-> b and an edge b-> a

# 0.952395 seconds (12.61 M allocations: 370.791 MiB, 38.39% gc time)
@time for r in eachrow(df)
    add_edge!(g, d2[r[1]] => d2[r[2]])
end

# am = adjacency_matrix(g)
# E = eigs(am)
# maximum_adjacency_visit(g)
# bc = betweenness_centrality(g)

xs = map(x -> d[x] => degree(g, x), vertices(g));
ys = map(x -> x => degree(g, x), vertices(g));
sort!(xs; by=last, rev=true)
sort!(ys; by=last, rev=true)
V = xs[1][1]
V = ys[1][1]
hood = neighborhood(g, xs[1][1], 1)
hood = neighborhood(g, V, 1)
ig, vs = induced_subgraph(g, hood)
bc = betweenness_centrality(ig)
dc = degree_centrality(ig)
id = d[V]
V2 = ys[2][1]
V2 in hood
id2 = d[V2]
hood2 = neighborhood(g, V2, 1)
ig2, vs = induced_subgraph(g, hood)
bc = betweenness_centrality(ig2)
dc = degree_centrality(ig2)


# g2 = SimpleWeightedGraph(length(vid))
# E = edgetype(g2)[]
# todo check that there isn't an edge from a-> b and an edge b-> a
# for r in eachrow(df)
# push!(E, edgetype(g2)(d2[r[1]], d2[r[2]], r[3]))

# end

# 179.321221 seconds (3.10 M allocations: 102.237 MiB, 0.03% compilation time)
# @time for e in E
#     add_edge!(g2, e)
# end
# add_e
# @time add_edge!(g, E)
# gdict = Dict("unweighted" => g, "weighed" => g2)
# open("HumanNetGraphs.lgz", "w") do io
#     Graphs.savelg_mult(io, gdict)
# end
# savegraph("HumanNetFCWeighted.lgz", g2, d, format=LGFormat)


# gs = loadgraphs("HumanNetGraphs.lgz")
# g = gs["unweighted"]

xs = map(x -> d[x] => degree(g, x), vertices(g));
ys = map(x -> x => degree(g, x), vertices(g));
sort!(xs; by=last, rev=true)
sort!(ys; by=last, rev=true)
V = xs[1][1]
V = ys[1][1]
hood = neighborhood(g, xs[1][1], 1)
hood = neighborhood(g, V, 1)
igg, vs = induced_subgraph(g, hood)
ig = copy(igg)
rem_vertex!(ig, 1)

function foo(g)
    v = vertices(g)
    d = degree(g)
    v .=> d, v, d
end
#deletes 0 degree edges

function prune(g)
    ig = copy(g)
    vd, v, d = foo(ig)
    rem_vertices!(ig, first.(filter(x -> last(x) == 0, vd)))
    vd, v, d = foo(ig)

    sort!(vd; by=last, rev=true)
    for (vv, dd) in vd
        # should be sorting by weighting here, not by least vid
        ns = neighbors(ig, vv)
        isempty(ns) && continue
        # n = rand(ns)
        n = first(ns)
        for e in setdiff(ns, [n])
            rem_edge!(ig, vv, e)
        end
    end
    vd, v, d = foo(ig)
    rem_vertices!(ig, first.(filter(x -> last(x) == 0, vd)))
    ig
end


"crossing point subgraph "
function crossing_subgraph(g)

    df2 = DataFrame(v=vertices(g), d=degree(g))
    sort!(df2, :d; rev=true)
    NV = 1
    NE = 0
    n = 1
    vids = nothing
    ig = nothing
    while NV > NE
        @info n, NV, NE
        ig, vids = induced_subgraph(g, df2.v[1:n])
        NV = nv(ig)
        NE = ne(ig)
        n += 1
    end
    ig
end

Graphs.getindex(G::AbstractGraph, bv::AbstractVector{Bool}) = G[findall(bv)]
dget(d, idxs) = map(x -> getindex(d, x), idxs)
# dget(d, vids)

# D = degree(G)
# ncc = []
# for i in 1:100
#     tmp = G[D.>=i]
#     push!(ncc, length(connected_components(tmp)))
# end

"a"
graphplot(ig; nlabels=repr.(dget(d, vids)), arrowsize=20)

const GENE_URL = "https://www.ncbi.nlm.nih.gov/gene/"

function get_gene_name(id)
    Cascadia.text(eachmatch(Selector("#summaryDl > dd.noline"), get_html(joinpath(GENE_URL, "?term=$id")).root)[1].children[1])
end


"Nest"
function nest_apply(f, x, n)
    for i in 1:n
        x = f(x)
    end
    x
end

"NestList"
function nest_apply_save(f, x, n)
    V = Vector{typeof(x)}(undef, n)
    for i in 1:n
        x = f(x)
        V[i] = x
    end
    V
end



ig = ig[degree(ig).>0]
ccs = connected_components(ig)
sort!(ccs; by=length, rev=true)
is_connected(ig[ccs[1]])

p = prune(g)
gccs = connected_components(p)
g[first(sort!(gccs; by=length, rev=true))]
gs = nest_apply_save(connected_prune, g, 40)
gs = largest_connected_subgraph.(gs)
gs = map(add_self_loops!, gs)
G = gs[end-4]
arity = round(Int, log2(nv(G)), RoundUp)
# to_func(boolean_function(G))

using GeneRegulatoryNetworks: to_func, get_vals
fex = boolean_function(G)
f_ = to_func(fex)
# vs = get_vals(f_, arity);
vs = get_vals(G)
# allequal(eachrow(reduce(hcat, vs)'))

fex = boolean_function(G)
n = get_n(G)
# @syms x::Bool[1:n]
x = Sym{Bool}.(Symbol.(:x, string.(1:n)))
fex2 = boolean_function(G, x);
# vs = get_vals(f_2, round(Int, log2(nv(G)), RoundUp));
ex = fex2[3]
str = eq_str_to_wl(string(ex))
sfex = simplify.(fex2);
tree = Tree(ex)

x = Sym{Bool}.(Symbol.(:x, string.(1:n)))
x1, x2, x3 = x[1], x[2], x[3]

using AbstractTrees
using AbstractTrees: printnode
struct SymNode
    ex
    children::Vector{SymNode}
end
ex1 = (x1 & x2)
ex2 = (x1 & (x2 | x3))
ex2 = (x1 & (x2 & x3))
a = arguments(ex2)
SymNode(ex1, children(ex1))
MyNode(ex, children=SymNode[arguments]) where {T} = MyNode{T}(op, children)


t2 = [1, [2, 3, [4, 5]]]
print_tree(t2)
T = typeof(ex2)
# AbstractTrees.children(node::T) = arguments(node)
AbstractTrees.children(node::T) = istree(node) ? arguments(node) : T[]
# AbstractTrees.children(node::T) = istree(node) ? arguments(node) : T[]

node_id(node) = istree(node) ? string(operation(node)) : nameof(node)
AbstractTrees.printnode(io::IO, node::T) = print(io, node_id(node))
printnode(stdout, ex2)

# f_2 = to_func(fex)
# allequal(eachrow(reduce(hcat, vs)'))

# allequal(eachrow(reduce(hcat, get_f(gg).(bs))))


# f2(x1, x2) = (x1 & x2, x1 & !x2)


# get_vals(f, n) = map(x -> f(x...), GeneRegulatoryNetworks.get_bools(n))
# get_vals(f2, 2)
# get_vals(to_func(boolean_function(ga)), 2)

# vs = get_vals(f_, arity);

# f_2(x1, x2) = (x1&x2, x1&x2)
# f_2(xs) = f_2(xs...)

g = SimpleGraph()

# function tree_to_graph(g, node)
#     # for (i, node) in PreOrderDFS(tree)
#     i = 0
#     id = node_id(ex2)
#     add_vertex!(g)
#     i += 1
#     parent_v = i

#     for k in children(ex2)
#         add_vertex!(g)
#         i += 1
#         add_edge!(g, parent_v, i)
#     end

#     # kids = 

# end

es = []
for (i, node) in PreOrderDFS(tree)
    push!(es, (i, node_id(node)), children(node))
end


using Graphs
g = SimpleGraph(5)
combinations(1:5, 2)
collect(combinations(1:5, 2))
e_labels = collect(combinations(1:5, 2))
1:binomial(5, 2)
e_ = 1:binomial(5, 2) .=> e_labels
e_labels
possible_gs = collect(combinations(collect(1:binomial(5, 2)), 5))
e_labels[possible_gs[1]]
for ev in e_labels[possible_gs[1]]
    add_edge!(g, ev...)
end
