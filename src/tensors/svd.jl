#=
    Singular value decomposition of tensors.
=#

### Perform a singular value decomposition 
"""
    svd(x, dims; kwargs...)
    svd(x, dim::Int; kwargs...)

Computer a singular value decomposition of tensor x. Seperates the dimensions
dims from the remainder.

# Key arguments

    - `cutoff::Float64=0.0`: Truncation criteria to reduce the bond dimension.
      Good values range from 1e-8 to 1e-14.
    - `mindim::Int=1`: Mininum dimension for truncated.
    - `maxdim::Int=0`: Maximum bond dimension for truncation. Set to 0 to have
      no limit.
    
"""
function svd(x, dims; cutoff::Float64=0.0, mindim::Int=1, maxdim::Int=0)
end

function svd(x, dim::Int; kwargs...)
    if dim == -1 dim = ndims(x) end 
    return svd(x, (dims); kwargs...)
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
    y = reshape(y, (prod(sinners), prod(sinners)))
    
    # Do the SVD using LinearAlgebra 
    local t
    try
        t = LinearAlgebra.svd(y, alg=LinearAlgebra.DivideAndConquer())
    catch e
        t = LinearAlgebra.svd(y, alg=LinearAlgebra.QRIteration())
    end
    U, S, V = t.U, t.S, t.Vt

    # Truncatation criteria 
    idx = findfirst(t.S .== 0)
    maxdim = max(min(maxdim, isnothing(idx) ? length(t.S) : idx-1), 1) 
    mindim = min(mindim, maxdim)
    if cutoff != 0.0

    end
end

### Permuting dimensions 
function _svd_permute_dims(x, dims)
    inner_dims = tuple(setdiff(1:ndims(x), dims)...)
    sinners = map(i->size(x, i), inner_dims)
    souters = map(i->size(x, i), dims)
    return sinners, souters, tuple(inner_dims..., outerdims...)
end

### Checks 
function _svd_checkdims(x, dims)
    if ndims(x) >= length(dims)
        throw(ArgumentError("The number of dimensions given exceed that of the tensor."))
    end
    if !all(i -> isinteger(i) , dims)
        throw(ArgumentError("Dimensions $(dims) are not integers."))
    end
    if !all(i -> ( i > 0 && i <= ndims(x)), dims)
        throw(ArgumentError("Dimensions $(dims) exceed the size of the tensor: $(ndims(x))."))
    end
end