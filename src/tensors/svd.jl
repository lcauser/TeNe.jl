#=
    Singular value decomposition of tensors.

    For future: this implementantion could use some work. For example, using cached memory
    and an implementation of SVD which provides U, S and V (that are preallocated).
=#

export tsvd

### Perform a singular value decomposition 
"""
    svd(x, dims; kwargs...)
    svd(x, dim::Int; kwargs...)

Computer a singular value decomposition of tensor `x`. Seperates the dimensions
`dims` from the remainder.

# Key arguments

    - `cutoff::Float64=0.0`: Truncation criteria to reduce the bond dimension.
      Good values range from 1e-8 to 1e-14.
    - `mindim::Int=1`: Mininum dimension for truncated.
    - `maxdim::Int=0`: Maximum bond dimension for truncation. Set to 0 to have
      no limit.
    
"""
function tsvd(x, dims; cutoff::Float64=0.0, mindim::Int=1, maxdim::Int=0)
    _svd_checkdims(x, dims)
    return _svd(x, dims; cutoff, mindim, maxdim)
end

function tsvd(x, dim::Int; kwargs...)
    if dim == -1 dim = ndims(x) end 
    return tsvd(x, Tuple(dim); kwargs...)
end

### Unsafe SVD 
function _svd(x, dims; cutoff::Float64=0.0, mindim::Int=1, maxdim::Int=0)
    # Permute and reshape 
    sinners, souters, pixs = _svd_permute_dims(x, dims)
    if !all(i->pixs[i]==i, 1:ndims(x))
        y = cache(x, (sinners..., souters...))
        permutedims!(y, x, pixs)
    else 
        y = x
    end
    y = reshape(y, (prod(sinners), prod(souters)))
    
    # Do the SVD using LinearAlgebra 
    local t
    try
        t = LinearAlgebra.svd(y, alg=LinearAlgebra.DivideAndConquer())
    catch e
        t = LinearAlgebra.svd(y, alg=LinearAlgebra.QRIteration())
    end

    # Truncatation criteria 
    idx = findfirst(t.S .== 0)
    maxdim = max(min(maxdim, isnothing(idx) ? length(t.S) : idx-1), 1) 
    mindim = min(mindim, maxdim)
    if cutoff != 0.0
        S2 = t.S .^ 2
        S2cum = reverse(cumsum(reverse(S2))) ./ sum(S2)
        idxs = findlast([x > cutoff for x = S2cum])
        idxs = isnothing(idxs) ? 1 : idxs
        maxdim = min(maxdim, idxs)
    end
    vals = max(maxdim, mindim)

    # Truncate 
    U = t.U[:, 1:vals]
    S = Diagonal(t.S[1:vals])
    V = t.Vt[1:vals, :]

    # Ungroup the dimensions
    U = reshape(U, (sinners..., vals))
    V = reshape(V, (vals, souters...))

    return U, S, V
end

### Permuting dimensions 
function _svd_permute_dims(x, dims)
    inner_dims = tuple(setdiff(1:ndims(x), dims)...)
    sinners = map(i->size(x, i), inner_dims)
    souters = map(i->size(x, i), dims)
    return sinners, souters, tuple(inner_dims..., dims...)
end

### Checks 
function _svd_checkdims(x, dims)
    if ndims(x) < length(dims)
        throw(ArgumentError("The number of dimensions given exceed that of the tensor."))
    end
    if !all(i -> isinteger(i) , dims)
        throw(ArgumentError("Dimensions $(dims) are not integers."))
    end
    if !all(i -> ( i > 0 && i <= ndims(x)), dims)
        throw(ArgumentError("Dimensions $(dims) exceed the size of the tensor: $(ndims(x))."))
    end
end