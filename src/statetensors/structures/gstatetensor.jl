#=
    Generalised state tensors describe a full-rank description of a tensor
    for some many-body system. The `generalised` here simply means that any rank 
    tensor, such as a vector or matrix, is described by this object.
    Functionality for specific structures, such as state vectors or state matrices
    can be found seperately.
=#

mutable struct GStateTensor{r,T<:AbstractArray} <: AbstractStateTensor
    tensor::T
end
export GStateTensor

### Indexing 
Base.getindex(ψ::GStateTensor, i...) = Base.getindex(ψ.tensor, i...)
Base.firstindex(ψ::GStateTensor) = Base.firstindex(ψ.tensor)
Base.lastindex(ψ::GStateTensor) = Base.lastindex(ψ.tensor)
Base.eachindex(ψ::GStateTensor) = Base.eachindex(ψ.tensor)

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
TeNe.rank(::GStateTensor{r}) where {r} = r

"""
    Base.length(::GStateTensor)

The length of a state tensor.
"""
Base.length(ψ::GStateTensor) = fld(ndims(ψ.tensor), rank(ψ))

tensor(ψ::GStateTensor) = ψ.tensor


### Dimensions
export ind, inds, dim, dims

"""
    ind(::GStateTensor, which::Int, site::Int)

Return the index of the tensor for dimension `which` at position `site`.
"""
function ind(::GStateTensor{r}, which::Int, site::Int) where {r}
    return r*(site-1)+which
end

"""
    inds(::GStateTensor, which::Int)

Return the indices of the tensor for dimension `which` at each site.
"""
function inds(ψ::GStateTensor{r}, which::Int) where {r}
    return Tuple(which:r:(r*length(ψ)))
end

"""
    dim(ψ::GStateTensor, [which::Int, site::Int])

Return the physical dimension of a StateTensor. The axis can be specified using
`which`, and furthermore the `site`. Returns `0` for heterogeneous systems.
"""
function dim(ψ::GStateTensor, which::Int, site::Int)
    return size(tensor(ψ), ind(ψ, which, site))
end

function dim(ψ::GStateTensor, which::Int)
    ds = dims(ψ, which)
    if all(map(j->j==ds[1], ds))
        return ds[1]
    else
        return 0
    end
end

function dim(ψ::GStateTensor)
    ds = size(tensor(ψ))
    if all(map(j->j==ds[1], ds))
        return ds[1]
    else
        return 0
    end
end

"""
    dims(ψ::GStateTensor, which::Int)

Return the size of the tensor for dimensions `which` at each site.
"""
function dims(ψ::GStateTensor, which::Int)
    return map(j->size(tensor(ψ), j), inds(ψ, which))
end

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


### Base operations
function *(ψ::GStateTensor, a::Number)
    ten = ψ.tensor .* a
    return GStateTensor(rank(ψ), dim(ψ), ten)
end
*(a::Number, ψ::GStateTensor) = *(ψ, a)
/(ψ::GStateTensor, a::Number) = *(ψ, 1/a)

function /(a::Number, ψ::GStateTensor)
    ten = a ./ ψ.tensor
    return GStateTensor(rank(ψ), dim(ψ), ten)
end

function +(ψ::GStateTensor, ϕ::GStateTensor)
    ten = ψ.tensor .+ ϕ.tensor
    return GStateTensor(rank(ψ), dim(ψ), ten)
end

function -(ψ::GStateTensor, ϕ::GStateTensor)
    ten = ψ.tensor .- ϕ.tensor
    return GStateTensor(rank(ψ), dim(ψ), ten)
end

function Base.show(io::IO, ψ::GStateTensor)
    println(io, "$(typeof(ψ))")
    println(io, Base.size(ψ.tensor))
end

Base.copy(ψ::GStateTensor) = typeof(ψ)(ψ.tensor)
Base.deepcopy(ψ::GStateTensor) = typeof(ψ)(Base.copy(ψ.tensor))

### Initalising 
function GStateTensor(rank::Int, dim::Int, length::Int; T::Type = ComplexF64)
    tensor = zeros(T, map(j->dim, Base.OneTo(rank*length))...)
    return GStateTensor{rank,typeof(tensor)}(tensor)
end

function GStateTensor(rank::Int, tensor::Q) where {Q<:AbstractArray}
    return GStateTensor{rank,typeof(tensor)}(tensor)
end

function randomgst(rank::Int, d::Int, N::Int; T::Type = ComplexF64)
    ψ = GStateTensor(rank, randn(T, map(j->d, Base.OneTo(rank*N))...))
    normalize!(ψ)
    return ψ
end

function productgst(N::Int, A::AbstractArray; T::Type = ComplexF64)
    tensor = ones(T)
    for i in Base.OneTo(N)
        tensor = tensorproduct(tensor, A; tocache = i!=N)
    end
    return GStateTensor(ndims(A), tensor)
end


### Save and write 
function HDF5.write(
    parent::Union{HDF5.File,HDF5.Group},
    name::AbstractString,
    M::GStateTensor{r},
) where {r}
    g = create_group(parent, name)
    attributes(g)["type"] = "StateTensor"
    attributes(g)["version"] = 1
    write(g, "rank", r)
    write(g, "tensor", M.tensor)
end


function HDF5.read(
    parent::Union{HDF5.File,HDF5.Group},
    name::AbstractString,
    ::Type{GStateTensor},
)
    g = open_group(parent, name)
    if read(attributes(g)["type"]) != "StateTensor"
        error("HDF5 group of file does not contain State Tensor data.")
    end
    tensor = read(g, "tensor")
    rank = read(g, "rank")
    return GStateTensor(rank, tensor)
end


### Conjugation 
export conj, isconj
struct ConjGStateTensor{r,T} <: GStateTensorTrait where {r,T}
    StateTensor::GStateTensor{r,T}
end
TeNe.conj(ψ::GStateTensor) = ConjGStateTensor(ψ)
TeNe.conj(ψ::ConjGStateTensor) = ψ.StateTensor
isconj(ψ::ConjGStateTensor) = true
