# GeneRegulatoryNetworks.jl

[![Build status](https://badge.buildkite.com/ba257503f8b1524ac709b1bcb2bedbeaedf5e5349247d01f4b.svg)](https://buildkite.com/julia-computing-1/generegulatorynetworks-dot-jl)

## Current Functionality

* Constructing directed graphs from boolean functions.
    
```julia 
f(x1, x2, x3, x4) = (x1 | x2, x1 & x4, !x1 & x4, !x3)
f(xs) = f(xs...)
g = sync_state_transition_digraph(f, 4)
# or 
f(x1, x2, x3, x4) = (x1 | !x2, !x1 | (x3 & x4), !x2, x1)
f(xs) = f(xs...)
g = async_state_transition_digraph(f, 4)

```
* Lowering boolean functions into ODESystem
```julia 
f(x1, x2) = (x1 | !x2, !x1 | x2)
f(xs) = f(xs...)
@named sys = GeneRegulatoryNetworks.boolean_f_to_ode(f, 2)
```
## Next Steps

1. Use a database of TFs to construct big networks for use in our automated model builder.

2. Hook it into the GUI


References:

http://www.mi.fu-berlin.de/en/math/groups/dibimath/PhD/schwieger/index.html

