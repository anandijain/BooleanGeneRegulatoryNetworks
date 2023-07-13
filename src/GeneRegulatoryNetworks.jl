module GeneRegulatoryNetworks
using SparseArrays
using Graphs, SimpleWeightedGraphs#, GraphPlot
using Symbolics, Catalyst
using Boolin, GraphHelpers
using Base.Iterators
import Base: union

include("boolean_network.jl")
include("lower.jl")
include("utils.jl")

export bitsdiff, bitscomm
export sync_state_transition_digraph, async_state_transition_digraph
export boolean_f_to_ode, boolean_f_from_sstg

end # module
