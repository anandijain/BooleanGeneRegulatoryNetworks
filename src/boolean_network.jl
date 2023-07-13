
# sync_state_transition_digraph(f) = sync_state_transition_digraph(f, union(Symbolics.get_variables.(f)...))

function sync_state_transition_digraph(f, x; showplot=false)
    n = length(x)
    nvertices = 2^n
    g = SimpleDiGraph(nvertices) # convention v1 = 0000 v2 = 0001 etc v16 = 1111
    sync_state_transition_digraph!(g, f, x, n; showplot=showplot)
    g
end

function sync_state_transition_digraph!(g, f, x, n; showplot=false)
    bs = bool_itr(n)
    fx = falses(length(x))
    built_f = build_function(f, x...; expression=Val{false})[2]
    for b in bs
        b = collect(b)
        built_f(fx, b...)
        vin_idx = from_bools(b) + 1
        vout_idx = from_bools(fx) + 1
        if b != fx # maybe self loops are ok
            add_edge!(g, vin_idx, vout_idx)
        end
    end
    showplot && display(gplot(g; nodelabel=bstrs))
    nothing
end

# async_state_transition_digraph(f) = async_state_transition_digraph(f, union(Symbolics.get_variables.(f)...))
function async_state_transition_digraph(f, x; showplot=false)
    n = length(x)
    nvertices = 2^n
    g = SimpleDiGraph(nvertices) # convention v1 = 0000 v2 = 0001 etc v16 = 1111
    async_state_transition_digraph!(g, f, x, n; showplot=showplot)
    g
end

function async_state_transition_digraph!(g, f, x, n; showplot=false)
    nvertices = nv(g)
    bs = bool_itr(n)
    fx = falses(length(x))
    built_f = build_function(f, x...; expression=Val{false})[2]
    for b in bs
        b = collect(b)
        built_f(fx, b...)
        bds = bitsdiff(b, fx)
        vin_idx = from_bools(b) + 1
        for bd in bds
            newv = copy(b)
            newv[bd] = fx[bd]
            vout_idx = from_bools(newv) + 1
            add_edge!(g, vin_idx, vout_idx)
        end
    end
    showplot && display(gplot(g; nodelabel=bstrs))
    nothing
end

boolean_f_from_sstg(g) = boolean_f_from_sstg(g, make_boolean_variables(round(Int, log2(nv(g)), RoundUp)))
function boolean_f_from_sstg(g, x)
    # @assert all(<=(1), outdegree(g))
    N = nv(g)
    nd = round(Int, log2(N), RoundUp)
    V = vertices(g)
    D = Boolin.binary_decomposition.(V .- 1, nd)
    @assert all(==(length(D[1])), length.(D))
    VD = Dict(V .=> D)
    dsts = [VD[e.dst] for e in collect(edges(g))]

    tups = Boolin._get_tups(x)
    M = Bool.(reduce(hcat, dsts)')

    map(x -> Boolin.boolean_function([], tups, x), eachcol(M))
end

local_interaction_graph(f, v) = SimpleWeightedDiGraph(GeneRegulatoryNetworks.jac(f, v))
global_interaction_graph(f, n) = union(map(b -> local_interaction_graph(f, b), bools(n)))
# global_interaction_graph(f) = global_interaction_graph(f, length(union(Symbolics.get_variables.(f)...)))

# intentional piracy and it's a little sus, because you lose some info eg union([A 1-> B, A -1 -> B]) = A B
union(gs::Vector{<:SimpleWeightedDiGraph}) = SimpleWeightedDiGraph(sum(weights.(gs)))

"""
```julia
g = SimpleDiGraph(sparse([4, 1, 1, 2, 3, 4, 4], [1, 2, 3, 4, 4, 5, 6], [1, 1, 1, 1, 1, 1, 1], 6, 6))
t = [5, 6]

GeneRegulatoryNetworks.istrapset(g, t) # true

add_edge!(g, 5, 4)

GeneRegulatoryNetworks.istrapset(g, t) # false

```
"""
function istrapset(g, t)
    V = vertices(g)
    V_t = setdiff(V, t)
    for v in t
        for Vt in V_t
            has_path(g, v, Vt) && return false
        end
    end
    true
end

"""
```julia
g = SimpleDiGraph(sparse([4, 1, 1, 2, 3, 4, 4], [1, 2, 3, 4, 4, 5, 6], [1, 1, 1, 1, 1, 1, 1], 6, 6))
gplot(g, nodelabel=1:nv(g))
t = [5, 6]

GeneRegulatoryNetworks.is_noreturnset(g, t) # true

rem_edge!(g, 4, 5)
rem_edge!(g, 4, 6)

gplot(g, nodelabel=1:nv(g))
GeneRegulatoryNetworks.is_noreturnset(g, t) # false

```
"""
function is_noreturnset(g, t)
    V = vertices(g)
    V_t = setdiff(V, t)
    for v in t
        for Vt in V_t
            has_path(g, Vt, v) && return true
        end
    end
    false
end

"
for f: {0, 1}^n -> {0, 1}^n, the sstg has nv(g) == 2^n.

the number of labeled digraphs in n nodes is 2^(2^n) / 2^n
"
function all_sstgs(n)
    x = make_boolean_variables(n)
    fs_itr = boolean_functions(x,n)
    S = Set{SimpleDiGraph}()
    # Gs = SimpleDiGraph[]
    for (i, f) in enumerate(fs_itr)
    # for f in fs_itr
        stg = sync_state_transition_digraph(collect(f), x)
        # push!(Gs, stg)
        i % 32000 == 0 && @info i, stg
        push_graph!(S, stg)
        # @info i
        length(S) == 100 && break
    end
    S#, Gs
    # Gs
end
