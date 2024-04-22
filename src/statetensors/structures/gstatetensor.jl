#=
    Generalised state tensors describe a full-rank description of a tensor
    for some many-body system. The `generalised` here simply means that any rank 
    tensor, such as a vector or matrix, is described by this object.
    Functionality for specific structures, such as state vectors or state matrices
    can be found seperately.
=#

mutable struct GStateTensor{r, T<:AbstractArray} <: AbstractStateTensor
    dim::Int
    tensor::T
end
export GStateTensor


### State tensor properties
export rank 
"""
    Base.eltype(::GStateTensor)

Returns the type of parameters within a state tensor.
"""
Base.eltype(ψ::GStateTensor) = Base.eltype(ψ.tensor)

"""
    rank(::GStateTensor)

Returns the rank of a state tensor.
"""
TeNe.rank(::GStateTensor{r, T}) where {r, T} = r

"""
    Base.length(::GStateTensor)

The length of a state tensor.
"""
Base.length(ψ::GStateTensor) = fld(ndims(ψ.tensor), rank(ψ))

### Norms 
export norm, normalize!
"""
    norm(ψ::GStateTensor)

Calculate the norm of a state tensor.
"""
function TeNe.norm(ψ::GStateTensor)
    return LinearAlgebra.norm(ψ.tensor)
end

"""
    normalize!(ψ::GStateTensor)

Normalize a state tensor.
"""
function TeNe.normalize!(ψ::GStateTensor)
    ψ.tensor .*= TeNe.norm(ψ)^-1
end