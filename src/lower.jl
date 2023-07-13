function interpolate_boolean_function(f)
    x = boolean_variables(f)
    n = length(x)
    interpolate_boolean_function(f, x, n, collect(first(@variables u[1:n])))
end

function interpolate_boolean_function(f, x, n, u)
    bs = bools(n)
    # fx = map(b -> map(f -> substitute(f, Dict(x .=> b)), f), bs)
    bf = BooleanFunction(f)
    fx = map(b->bf(x, b), bs)
    _tt = Dict(bs .=> fx)
    ts = Num[]
    for i in 1:n
        t = Num(0)
        for b in bs
            fi_x = _tt[b][i]
            t += fi_x * itp_term(b, u)
        end
        push!(ts, t)
    end
    ts
end

"""
come up w better name
"""
function itp_term(x, u)
    @assert length(x) == length(u)
    prod(x[j] * u[j] + (1 - x[j]) * (1 - u[j]) for j in 1:length(x))
end

"""

might want to allow option for normalized hill or no-hill functions

also maybe move the @parameters and @variables to be args
"""
function boolean_f_to_ode(f, x; kw...)
    n = length(x)
    @parameters t d[1:n] K[1:n, 1:n] hilln[1:n, 1:n]
    @variables u(t)[1:n]
    ts = interpolate_boolean_function(f, x, n, u)
    D = Differential(t)
    eqs = Equation[]
    for i in 1:n
        for j in 1:n
            h = Catalyst.hill(u[j], 1, K[i, j], hilln[i, j])
            ts[i] = substitute(ts[i], Dict(u[j] => h))
        end
        push!(eqs, D(u[i]) ~ d[i] * (ts[i] - u[i]))
        # push!(eqs, D(u[i]) ~ ts[i]) # no lifetime or degradation parameters
    end
    ODESystem(eqs; kw...)
end
