#=
    Taking the exponential of tensors; works by reshaping into a square matrix,
    exponentiating, and returning it back to the original form.
=#

export exp

### Functions to exponentiate 
"""
    exp(x, innerdims, outerdims; kwargs...)
    exp(x, outerdims; kwargs...)

Calculate the exponential over a set of outer dimensions. 
If unspecified, the innerdims will be assumed to be the remaining dimensions in order.

# Key arguments

    - `prefactor`: multiply the matrix by a prefactor before exponetiating.
"""
function exp(x, outerdims; prefactor=1)
    innerdims = setdiff(Base.OneTo(ndims(x)), outerdims)
    pxs, pxs_return, sinners, souters = _exp_permuted_dimensions(x, innerdims, outerdims)
    _exp_check_dims(sinners, souters)
    return _exp(x, pxs, pxs_return, souters, sinners, prefactor)
end
exp(x, outerdim::Int; kwargs...) = exp(x, (outerdim); kwargs...)


function exp(x, innerdims, outerdims; prefactor=1)
    _exp_check_inds(x, innerdims, outerdims)
    pxs, pxs_return, sinners, souters = _exp_permuted_dimensions(x, innerdims, outerdims)
    _exp_check_dims(sinners, souters)
    return _exp(x, pxs, pxs_return, souters, sinners, prefactor)
end
exp(x, innerdim::Int, outerdim::Int; kwargs...) = exp(x, (innerdim), (outerdim); kwargs...)

### Unsafe exponential 
function _exp(x, pxs, pxs_return, souters, sinners, prefactor)
    # Permute the matrix 
    y = cache(x, map(i->size(x, i), pxs))
    permutedims!(y, x, pxs)

    # Exponentiate
    y = reshape(y, (prod(souters), prod(sinners)))
    if prefactor != 1
        y .*= prefactor
    end
    y = LinearAlgebra.exp(y)
    y = reshape(y, (souters..., sinners...))

    # Permute back 
    z = cache(x)
    permutedims!(z, y, pxs_return)
    y = reshape(y, size(z))
    y .= z
    return y
end

### Find the remaining dimensions
function _exp_permuted_dimensions(x, innerdims, outerdims)
    sinners = map(i -> size(x, i), outerdims)
    souters = map(i -> size(x, i), innerdims)
    pxs = (innerdims..., outerdims...)
    pxs_return = map(i -> findfirst(i .== pxs), Base.OneTo(ndims(x)))
    return pxs, pxs_return, sinners, souters
end

### Checks 
function _exp_check_dims(sinners, souters)
    if prod(sinners) != prod(souters)
        throw(ArgumentError("Dimensions $(sinners) and $(souters) do not give a square matrix."))
    end
end

function _exp_check_inds(x, innerdims, outerdims)
    if !all(i -> !(i in outerdims), innerdims)
        throw(ArgumentError("The inner and outer dimensions have shared indices."))
    end
    if (length(innerdims) != length(outerdims)) || (length(innerdims) != (ndims(x) / 2))
        throw(ArgumentError("The number of dimensions do not match."))
    end
end