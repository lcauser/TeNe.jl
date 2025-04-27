#=
    Calculate Traces of StateOperators
=#

### Trace of StateOperators 
export trace
"""
    trace(Os::StateOperator...)

Compute the trace of a string of StateOperators.

# Examples 

```jldoctest
julia> O = productso(10, [0 1; 1 0]);
julia> trace(O, O)
1024.0 + 0.0im
```
"""
function trace(Os::StateOperator...)
    _op_trace_validation(Os...)

    # Deal with length one seperately; this needs work later...
    if length(Os) == 1
        ten = tensor(Os[1])
        for _ = 1:length(Os[1])
            ten = trace(ten, 1, 2)
        end
        return ten[]
    end

    # Create a tensor from the cache
    T = _promote_tensor_eltype(Os...)
    perms = tuple(map(j->isodd(j) ? j + 1 : j - 1, Base.OneTo(ndims(tensor(Os[1]))))...)
    if istranspose(Os[1])
        dims = map(j->size(tensor(Os[1]), j), perms)
        ten = cache(T, dims, 2, 1)
        permutedims!(ten, tensor(Os[1]), perms)
    else
        dims = size(tensor(Os[1]))
        ten = cache(T, dims, 2, 1) .= tensor(Os[1])
    end
    if isconj(Os[1])
        ten .= conj.(ten)
    end

    # Contract with center operators 
    perms2 = _so_so_product_perms(length(Os[1]))
    for i in range(2, length(Os)-1)
        ten = contract(
            ten,
            tensor(Os[i]),
            collect(2:2:ndims(ten)),
            collect((istranspose(Os[i]) ? 2 : 1):2:ndims(tensor(Os[i]))),
            false,
            isconj(Os[i]),
        )
        ten = permutedims(ten, perms2)
    end

    # Contract with the last 
    dims = Tuple(Base.OneTo(ndims(ten)))
    if istranspose(Os[end])
        return contract(ten, tensor(Os[end]), dims, dims, false, isconj(Os[end]))[]
    else
        return contract(ten, tensor(Os[end]), dims, perms, false, isconj(Os[end]))[]
    end
end
