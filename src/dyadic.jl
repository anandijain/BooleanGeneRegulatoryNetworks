rademacher(i, w) = 2*w[i] - 1
rademacher(i) = w->rademacher(i, w)
# rademacher() = w->rademacher(i, w)

sn(w) = mapreduce(i->rademacher(i, w), +, 1:length(w))
random_walk_positions(n) = map(sn, GeneRegulatoryNetworks.gen_bools(n))

xs = []
nzs = []
for i in 1:10
    ps = random_walk_positions(i)
    nz = count(isequal(0), ps[1:length(ps) ÷ 2])
    
    push!(xs, ps)
    push!(nzs, nz)
end


nzss(n) = binomial(2*n+1, n+1)