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
    dim(::GStateTensor)

Returns the physical dimension of a state tensor. Returns `0` for heterogeneous 
systems (i.e. an invariant physical dimension).
"""
TeNe.dim(ψ::GStateTensor) = ψ.dim

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


### Base operations
import Base.+, Base.-, Base.*, Base./

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

Base.copy(ψ::GStateTensor) = typeof(ψ)(dim(ψ), ψ.tensor)
Base.deepcopy(ψ::GStateTensor) = typeof(ψ)(Base.copy(dim(ψ)), Base.copy(ψ.tensor))

### Initalising 
function GStateTensor(rank::Int, dim::Int, length::Int; T::Type=ComplexF64)
    tensor = zeros(T, map(j->dim, Base.OneTo(rank*length))...)
    return GStateTensor{rank, typeof(tensor)}(dim, tensor)
end

function GStateTensor(rank::Int, dim::Int, tensor)
    return GStateTensor{rank, typeof(tensor)}(dim, tensor)
end

export randomgst
function randomgst(rank::Int, d::Int, N::Int; T::Type=ComplexF64)
    ψ = GStateTensor(rank, d, rand(T, map(j->d, Base.OneTo(rank*N))...))
    normalize!(ψ)
    return ψ
end


### Save and write 
function HDF5.write(parent::Union{HDF5.File, HDF5.Group}, name::AbstractString,
                    M::GStateTensor{r, T}) where {r, T}
    g = create_group(parent, name)
    attributes(g)["type"] = "StateTensor"
    attributes(g)["version"] = 1
    write(g, "rank", r)
    write(g, "dim", dim(M))
    write(g, "tensor", M.tensor)
end


function HDF5.read(parent::Union{HDF5.File, HDF5.Group}, name::AbstractString,
                   ::Type{GStateTensor})
    g = open_group(parent, name)
    if read(attributes(g)["type"]) != "StateTensor"
        error("HDF5 group of file does not contain State Tensor data.")
    end
    dim = read(g, "dim")
    tensor = read(g, "tensor")
    rank = read(g, "rank")
    return GStateTensor(rank, dim, tensor)
end


### Conjugation 
export conj, isconj 
struct ConjGStateTensor{r, T} <: GStateTensorTrait where {r, T}
    StateTensor::GStateTensor{r, T}
end
TeNe.conj(ψ::GStateTensor) = ConjGMPS(ψ)
TeNe.conj(ψ::ConjGStateTensor) = ψ.StateTensor
isconj(ψ::ConjGStateTensor) = true