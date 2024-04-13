"""
    GMPS(rank::Int, dim::Int, tensors::Vector{Array{Complex{Float64}}},
         center::Int)
    GMPS(rank::Int, dim::Int, length::Int)

Create a generalised matrix product state with physical dimension dim.
"""

mutable struct GMPS{r} <: AbstractMPS
    dim::Int
    tensors::Vector{<:AbstractArray}
    center::Int
end
export GMPS

function GMPS(rank::Int, d::Int, N::Int; T::Type=ComplexF64)
    tensors = [zeros(T, (1, [d for __ = Base.OneTo(rank)]..., 1)) for _ = Base.OneTo(N)]
    return GMPS{rank}(d, tensors, 0)
end

### Indexing an MPS 
"""
    Base.getindex(ψ::GMPS, i::Int)

Returns the tensor of GMPS `ψ` at site `i`.
"""
Base.getindex(ψ::GMPS, i::Int) = ψ.tensors[i]

"""
    Base.firstindex(ψ::GMPS)

Returns the first index of GMPS `ψ`.
"""
Base.firstindex(ψ::GMPS) = Base.firstindex(ψ.tensors)

"""
    Base.lastindex(ψ::GMPS)

Returns the last index of GMPS `ψ`.
"""
Base.lastindex(ψ::GMPS) = Base.lastindex(ψ.tensors)

"""
    Base.eachindex(ψ::GMPS)

Returns each index within an MPS `ψ`.
"""
Base.eachindex(ψ::GMPS) = Base.eachindex(ψ.tensors)

"""
    Base.setindex!(ψ::GMPS, x, i::Int)

Set the tensor at site `i` of GMPS `ψ` to `x`.
"""
function Base.setindex!(ψ::GMPS, x, i::Int)
    ψ.tensors[i...] = x
end

### MPS Properties 
export rank, dim, center, bonddim, maxbonddim

"""
    Base.eltype(::GMPS)

Returns the type of parameters within a GMPS.
"""
Base.eltype(ψ::GMPS) = Base.eltype(ψ[begin])

"""
    LinearAlgebra.rank(::GMPS)

Returns the rank of an MPS object.
"""
LinearAlgebra.rank(::GMPS{r}) where {r} = r

"""
    dim(::GMPS)

The size of the physical dimensions in a GMPS. Returns `0` for heterogeneous 
systems (i.e. an invariant physical dimension).
"""
dim(ψ::GMPS) = ψ.dim

"""
    Base.length(::GMPS)

The length of a GMPS.
"""
Base.length(ψ::GMPS) = Base.length(ψ.tensors)

"""
    center(::GMPS)

The orthogonal center of an MPS. Returns `0` if nothing is set.
"""
center(ψ::GMPS) = ψ.center

"""
    bonddim(::GMPS, idx::Int)

Return the bond dimension size between idx and idx + 1. Returns nothing if
out of range.
"""
function bonddim(ψ::GMPS, site::Int)
    (site < 1 || site > Base.length(ψ)) && return nothing
    return size(ψ[site+1])[1]
end


"""
    maxbonddim(::GMPS)

Calculate the maximum bond dimension within an GMPS.
"""
function maxbonddim(ψ::GMPS)
    D = 0
    for i = Base.OneTo(length(ψ)-1)
        D = max(D, bonddim(ψ, i))
    end
    return D
end


### Norms
"""
    LinearAlgebra.norm(ψ::GMPS)

Calculate the norm of an GMPS.
"""
function LinearAlgebra.norm(ψ::GMPS)
    if center(ψ) == 0
        movecenter!(ψ, 1)
    end
    return LinearAlgebra.norm(ψ[center(ψ)])
end

"""
    normalize!(ψ::GMPS)

Normalize a GMPS.
"""
function LinearAlgebra.normalize!(ψ::GMPS)
    if center(ψ) == 0
        movecenter!(ψ, 1)
    end
    ψ[center(ψ)] .*= LinearAlgebra.norm(ψ)^-1
end


### Moving the orthogonal center
export moveleft!, moveright!, movecenter!

"""
    moveleft!(ψ::GMPS, idx::Int; kwargs...)

Move the gauge from a tensor within the GMPS to the left.
"""
function moveleft!(ψ::GMPS, idx::Int; kwargs...)
    if 1 < idx && idx <= length(ψ)
        U, S, V = svd(ψ[idx], 1; kwargs...)
        V = contract(S, V, 2, 1)
        ψ[idx] = U
        ψ[idx-1] = contract(ψ[idx-1], V, 2+rank(ψ), 2)
    end
end

"""
    moveright!(ψ::GMPS, idx::Int; kwargs...)

Move the gauge from a tensor within the GMPS to the right.
"""
function moveright!(ψ::GMPS, idx; kwargs...)
    if 0 < idx && idx < length(ψ)
        U, S, V = svd(ψ[idx], 2+rank(ψ); kwargs...)
        V = contract(S, V, 2, 1)
        ψ[idx] = U
        ψ[idx+1] = contract(V, ψ[idx+1], 2, 1)
    end
end

"""
    movecenter!(ψ::GMPS, idx::Int; kwargs...)

Move the orthogonal center of an GMPS.
"""
function movecenter!(ψ::GMPS, idx::Int; kwargs...)
    (idx < 1 || idx > length(ψ)) && error("The index is out of range.")
    if center(ψ) == 0
        for i = 1:idx-1
            moveright!(ψ, i; kwargs...)
        end
        N = length(ψ)
        for i = 1:N-idx
            moveleft!(ψ, N+1-i; kwargs...)
        end
    else
        if idx > center(ψ)
            for i = center(ψ):idx-1
                moveright!(ψ, i; kwargs...)
            end
        elseif idx < center(ψ)
            for i = 1:center(ψ)-idx
                moveleft!(ψ, center(ψ)+1-i; kwargs...)
            end
        end
    end
    ψ.center = idx
end


### Information operators 
function Base.show(io::IO, ψ::GMPS)
    println(io, "$(typeof(ψ))")
    for i = Base.OneTo(Base.length(ψ))
        println(io, "[$(i)] $(Base.size(ψ[i]))")
    end
end

### Creating copies
Base.copy(ψ::GMPS) = typeof(ψ)(LinearAlgebra.rank(ψ), dim(ψ), ψ.tensors, center(ψ))
Base.deepcopy(ψ::GMPS) = typeof(ψ)(Base.copy(LinearAlgebra.rank(ψ)), Base.copy(dim(ψ)), Base.copy(ψ.tensors),
                                        Base.copy(center(ψ)))
