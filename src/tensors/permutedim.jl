#=
    Permuting a single dimension within a tensor.
=#

export permutedim, permutedim!

"""
    permutedim(x, i, j; kwargs...)

Permute dimension with position `i` to position `j` for tensor `x`.

# Optional Keyword Arguments
    
    - `tocache::Bool=true`: store the result in the second level of the cache?
    - `sublevel=:auto`: if stored in cache, at which sublevel? :auto finds non-aliased memory

# Examples 

```jldoctest
julia> x = randn(ComplexF64, 2, 3, 4, 5);
julia> x = permutedim(x, 2, 4);
julia> size(x)
(2, 4, 5, 3)
```
"""
function permutedim(x, i::Int, j::Int; tocache::Bool=true, sublevel=:auto)
    i == j && return copy(x)
    if j==-1 j=ndims(x) end
    _permutedim_checkbounds(x, i, j)
    order = _permutedim_ordering(x, i, j)
    dims = map(k->size(x, k), order)
    if tocache
        if sublevel == :auto 
            sublevel = 1
            z = cache(eltype(x), dims, 2, sublevel; backend=typeof(get_backend(x)))
            if prod(dims) == length(x)
                while Base.mightalias(x, z)
                    sublevel += 1
                    z = cache(eltype(x), dims, 2, sublevel; backend=typeof(get_backend(x)))
                end
            end
        else
            z = cache(eltype(x), dims, 2, sublevel; backend=typeof(get_backend(x)))
        end
    else
        z = zeros(eltype(x), dims...)
    end
    permutedims!(z, x, order)
    return z
end

"""
    permutedim!(z, x, i, j; kwargs...)

Permute dimension `i` to `j` for tensor `x`. Store the result in `z`.
In place version of permutedim.

# Examples 

```jldoctest
julia> x = randn(ComplexF64, 2, 3, 4, 5);
julia> z = similar(x, (2, 4, 5, 3));
julia> permutedims!(z, x, 2, 4);
```
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


### permutedims(...) with automatic allocation
function _permutedims(x, order; sublevel=:auto)
    dims = map(k->size(x, k), order)
    if sublevel == :auto 
        sublevel = 1
        z = cache(eltype(x), dims, 2, sublevel; backend=typeof(get_backend(x)))
        if prod(dims) == length(x)
            while Base.mightalias(x, z)
                sublevel += 1
                z = cache(eltype(x), dims, 2, sublevel; backend=typeof(get_backend(x)))
            end
        end
    else
        z = cache(eltype(x), dims, 2, sublevel; backend=typeof(get_backend(x)))
    end
    permutedims!(z, x, order)
    return z
end