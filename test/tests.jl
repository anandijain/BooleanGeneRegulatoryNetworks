using GeneRegulatoryNetworks, Test
using Graphs, SimpleWeightedGraphs#, GraphPlot
using Catalyst
using Boolin
using OrdinaryDiffEq

@variables x1 x2 x3 x4
x = [x1, x2, x3, x4]
# example 2.3
f = [x1 | x2, x1 & x4, !x1 & x4, !x3]

g = sync_state_transition_digraph(f, x)
@test !Graphs.is_connected(g)

# example 1.1
f = [x1 | !x2, !x1 | (x3 & x4), !x2, x1]

g = async_state_transition_digraph(f, x)
T = [10, 12, 16, 14]
@test GeneRegulatoryNetworks.istrapset(g, T)

# TODO: add test/assert that hamming distance for every edge is 1 in async 
bf = BooleanFunction(f)
v = falses(4)
@test GeneRegulatoryNetworks.partial(bf, v, 3, 2) == -1

g = GeneRegulatoryNetworks.global_interaction_graph(bf, 4)
@test (nv(g) == 4) && (ne(g) == 7)
@test g isa SimpleWeightedDiGraph

# example 2.2
f = [x1 | !x2, !x1 | x2]
x = [x1, x2]
u = collect(only(@variables u[1:2]))
ts = GeneRegulatoryNetworks.interpolate_boolean_function(f)
@test isequal(ts, [(1 - u[1]) * (1 - u[2]) + (1 - u[2]) * u[1] + u[1] * u[2], (1 - u[1]) * (1 - u[2]) + (1 - u[1]) * u[2] + u[1] * u[2]])

@named sys = GeneRegulatoryNetworks.boolean_f_to_ode(f, x)
@test sys isa ODESystem

x = make_boolean_variables(3)
g = cycle_digraph(8)
g2 = sync_state_transition_digraph(boolean_f_from_sstg(g, x), x)
@test g == g2

f = boolean_function.([30, 40, 50, 60], 4)
x = make_boolean_variables(4)
g = sync_state_transition_digraph(f, x)
@named sys = boolean_f_to_ode(f, x)
prob = ODEProblem(sys, ones(4), (0, 100.), rand(36))
prob = ODEProblem(sys, ones(4), (0, 100.), ones(36))
sol = solve(prob, Rosenbrock23())
plot(sol)
