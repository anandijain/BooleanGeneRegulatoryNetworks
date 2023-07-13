@register_symbolic (Base.:&)(x, y)::Bool
# @register_symbolic (Base.:&)(xs)::Bool
@register_symbolic (Base.:|)(x, y)::Bool
# @register_symbolic (Base.:|)(x...)::Bool
# @register_symbolic (Base.:|)(xs)::Bool
@register_symbolic (Base.:!)(x)::Bool

# fix for set of BitVectors
bitsdiff(a, b) = findall(a .!= b)
bitscomm(a, b) = findall(a .== b)

hamming(k, l) = sum(k .!= l)

# defn 1.24
function flipsome(v, idxs)
    for idx in idxs
        v[idx] = !v[idx]
    end
    v
end

function partial(f, v, i, j)
    fv = f(v)
    # fv = map(f -> substitute(f, Dict(x .=> v)), f)
    v_j = flipsome(copy(v), j)
    fv_j = f(v_j)
    # fv_j = map(f -> substitute(f, Dict(x .=> v_j)), f)
    (fv_j[i] - fv[i]) / (v_j[j] - v[j])
end

gradient(f, v, i) = map(j -> partial(f, v, i, j), 1:length(v))
jac(f, v) = mapreduce(i -> gradient(f, v, i), hcat, 1:length(v))
