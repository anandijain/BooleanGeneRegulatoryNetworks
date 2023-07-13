using GeneRegulatoryNetworks, Boolin, Symbolics, Graphs, GraphHelpers
using Base.Iterators
using Symbolics: substituter, symtype, value

nb(n) = 2^(2^n)
nb(n, m) = nb(n)^m
nbn(n) = nb(n,n)

# trying to figure out if the number of unique state transition graphs is
nbu(n, m) = binomial(nb(n), m)
# or
# nbu(n,m) = nb(n,m)/factorial(m)
nbu(n) = nbu(n,n)


