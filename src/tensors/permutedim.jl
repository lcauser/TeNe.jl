#=
    Permuting a single dimension within a tensor.
=#

export permutedim, permutedim!

"""
    permutedim(x, i, j; kwargs...)

Permute dimension i to j for tensor x.

# Optional Keyword Arguments
    
    - 'tocache::Bool=false': store the result in the second level of the cache?
    - 'sublevel::Int=1': if stored in cache, at which sublevel?
"""
function permutedim(x, i::Int, j::Int; tocache::Bool=false, sublevel::Int=1)
    i == j && return copy(x)
    if j==-1 j=ndims(x) end
    _permutedim_checkbounds(x, i, j)
    order = _permutedim_ordering(x, i, j)
    dims = (size(x, k) for k in order)
    if tocache
        z = cache(eltype(x), dims, 2, sublevel; backend=typeof(get_backend(x)))
    else
        z = zeros(eltype(x), dims...)
    end
    permutedims!(z, x, order)
    return z
end

"""
    permutedim!(z, x, i, j; kwargs...)

Permute dimension i to j for tensor x. Store the result in z.
In place version of permutedim.
"""
function permutedim!(z, x, i::Int, j::Int)
    i == j && return nothing
    if j==-1 j=ndims(x) end
    _permutedim_checkbounds(x, i, j)
    order = _permutedim_ordering(x, i, j)
    _permutedim_checkresult(z, x, order)
    permutedims!(z, x, order)
end


### New ordering 
function _permutedim_ordering(x, i, j)
    # Any room for optimising so that a vector is not created?
    order = collect(1:ndims(x))
    dir = j > i ? 1 : -1
    for k = 1:abs(i-j)
        temp = order[i + dir*k]
        order[i + dir*k] = order[i + dir*(k-1)]
        order[i + dir*(k-1)] = temp
    end

    return tuple(order...)
end


### Checks 
function _permutedim_checkbounds(x, i, j)
    if i < 1 || i > ndims(x)
        throw(ArgumentError("Dimension $(i) is out of bounds for number of dimensions $(ndims(x))."))
    end
    if j < 1 || j > ndims(x)
        throw(ArgumentError("Dimension $(j) is out of bounds for number of dimensions $(ndims(x))."))
    end
end

function _permutedim_checkresult(z, x, dims)
    if !all(i->(size(z, i)==size(x, dims[i])), 1:ndims(x))
        throw(ArgumentError("Destination tensor has the wrong dimensions."))
    end
end


