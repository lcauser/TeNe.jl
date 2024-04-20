#=
    Generalised matrix product states describe a decomposition of one-dimensional
    space into a product of tensors, where each tensor has the physical dimensions
    for a local subsystem. The generalised here means that any rank tensor can be
    decomposed into MPS representation: e.g., a vector is an MPS, operator is an
    MPO... Functionality for specific structures, such as an MPS, can be found
    seperately.
=#

"""
    GMPS(rank::Int, dim::Int, tensors::Vector{Array{Complex{Float64}}},
         center::Int)
    GMPS(rank::Int, dim::Int, length::Int)

Create a generalised matrix product state with physical dimension dim.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensors.
"""

mutable struct GMPS{r} <: AbstractMPS
    dim::Int
    tensors::Vector{<:AbstractArray}
    center::Int
end
export GMPS

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

#Base.collect(ψ::GMPS) = ψ

### MPS Properties 
export rank, dim, center, bonddim, maxbonddim

"""
    Base.eltype(::GMPS)

Returns the type of parameters within a GMPS.
"""
Base.eltype(ψ::GMPS) = Base.eltype(ψ[begin])

"""
    rank(::GMPS)

Returns the rank of an MPS object.
"""
TeNe.rank(::GMPS{r}) where {r} = r

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
export norm, normalize!
"""
    norm(ψ::GMPS)

Calculate the norm of an GMPS.
"""
function TeNe.norm(ψ::GMPS)
    if center(ψ) == 0
        movecenter!(ψ, 1)
    end
    return LinearAlgebra.norm(ψ[center(ψ)])
end

"""
    normalize!(ψ::GMPS)

Normalize a GMPS.
"""
function TeNe.normalize!(ψ::GMPS)
    if center(ψ) == 0
        movecenter!(ψ, 1)
    end
    ψ[center(ψ)] .*= TeNe.norm(ψ)^-1
end


### Moving the orthogonal center
export movecenter!

function _moveleft!(ψ::GMPS, idx::Int; kwargs...)
    # Maybe TODO later: re-use memory
    if 1 < idx && idx <= length(ψ)
        U, S, V = tsvd(ψ[idx], Base.range(2, ndims(ψ[idx])); kwargs...)
        U = contract(U, S, 2, 1)
        ψ[idx] = V
        ψ[idx-1] = contract(ψ[idx-1], U, 2+rank(ψ), 1; tocache=false)
    end
end

function _moveright!(ψ::GMPS, idx; kwargs...)
    if 0 < idx && idx < length(ψ)
        # Maybe TODO later: re-use memory
        U, S, V = tsvd(ψ[idx], 2+rank(ψ); kwargs...)
        V = contract(S, V, 2, 1)
        ψ[idx] = U
        ψ[idx+1] = contract(V, ψ[idx+1], 2, 1; tocache=false)
    end
end

"""
    movecenter!(ψ::GMPS, idx::Int; kwargs...)

Move the orthogonal center of an GMPS `ψ` to center `idx`.

# Optional Keyword Arguments

    - `cutoff::Float64=0.0`: Truncation criteria to reduce the bond dimension.
      Good values range from 1e-8 to 1e-14.
    - `mindim::Int=1`: Mininum dimension for truncated.
    - `maxdim::Int=0`: Maximum bond dimension for truncation. Set to 0 to have
      no limit.
"""
function movecenter!(ψ::GMPS, idx::Int; kwargs...)
    # Checks bounds 
    if idx < 1 || idx > length(ψ)
        throw(BoundsError("Center index is out of bounds."))
    end
    
    # Move center
    if center(ψ) == 0
        for i = Base.OneTo(idx-1)
            _moveright!(ψ, i; kwargs...)
        end
        for i = length(ψ):-1:idx+1
            _moveleft!(ψ, i; kwargs...)
        end
    else
        if idx > center(ψ)
            for i = center(ψ):idx-1
                _moveright!(ψ, i; kwargs...)
            end
        elseif idx < center(ψ)
            for i = center(ψ):-1:idx+1
                _moveleft!(ψ, i; kwargs...)
            end
        end
    end
    ψ.center = idx
end

### Replace the sites within a MPS with a contraction of multiple sites
"""
    replacesites!(ψ::GMPS, A, site::Int, direction::Bool=false, normalize::Bool=false; kwargs...)

Replace the tensors over a small range of sites in a GMPS.

# Arguments 

    - `ψ`: The GMPS.
    - `A': The updated composite tensor.
    - `site::Int`: The first site in the range which is replaced.
    - `direction::Bool=false`: The sweeping direction; `true` is sweeping left, `false` sweeping right.
    - `normalize::Bool=false`: Normalize the GMPS after replacing the sites?

# Optional Keyword Arguments

    - `cutoff::Float64=0.0`: Truncation criteria to reduce the bond dimension.
      Good values range from 1e-8 to 1e-14.
    - `mindim::Int=1`: Mininum dimension for truncated.
    - `maxdim::Int=0`: Maximum bond dimension for truncation. Set to 0 to have
      no limit.
"""
function replacesites!(ψ::GMPS, A, site::Int, direction::Bool=false, normalize::Bool=false; kwargs...)
    # Determine the number of sites
    nsites = (ndims(A) - 2) / rank(ψ)

    # Deal with case of just one site
    if nsites == 1
        ψ[site] = A
        ctr = site + 1 - 2*direction
        if 0 < ctr && ctr <= length(ψ)
            movecenter!(ψ, ctr)
        end
        if normalize
            normalize!(ψ)
        end
    end

    # Repeatidly apply SVD to split the tensors
    U = A
    for i = 1:nsites-1
        if direction
            ### Sweeping Left
            # Group together the last indices
            idxs = collect(length(size(U))-rank(ψ):length(size(U)))
            U, cmb = combineidxs(U, idxs)

            # Find the next site to update
            site1 = site+nsites-i

            # Apply SVD and determine the tensors
            U, S, V = tsvd(U, ndims(U); kwargs...)
            U = contract(U, S, length(size(U)), 1)
            D = site1 == length(ψ) ? 1 : size(ψ[site1+1])[1]
            dims = (size(S)[2], [dim(ψ) for i=1:rank(ψ)]..., D)
            V = reshape(V, dims)

            # Update tensors
            ψ[site1] = V
        else
            ### Sweeping right
            # Combine first two indexs together
            idxs = collect(1:1+rank(ψ))
            U, cmb = combineidxs(U, idxs)

            # Find the next site to update
            site1 = site + i - 1

            # Apply SVD and determine the tensors
            U, S, V = svd(U, ndims(U); kwargs...)
            U = contract(U, S, length(size(U)), 1)
            U = moveidx(U, length(size(U)), 1)
            D = site1 == 1 ? 1 : size(ψ[site1-1])[2+rank(ψ)]
            dims = (size(S)[2], D, [dim(ψ) for i=1:rank(ψ)]...)
            V = reshape(V, dims)
            V = moveidx(V, 1, -1)

            # Update tensors
            ψ[site1] = V
        end
    end

    # Update the final site
    site1 = direction ? site : site + nsites - 1
    ψ[site1] = U
    ψ.center = site1
    if normalize
        normalize!(ψ)
    end
end

### Products with numbers
import Base.*, Base./
function *(ψ::GMPS, a::Number)
    ϕ = deepcopy(ψ)
    if center(ψ) != 0
        ϕ.tensors[center(ϕ)] .*= a
    else
        ϕ.tensors[1] .*= a
    end
    return ϕ
end
*(a::Number, ψ::GMPS) = *(ψ, a)
/(ψ::GMPS, a::Number) = *(ψ, 1/a)

### Addition and subtraction 
# TODO later: add variational option
import Base.+, Base.-

function +(ψ::GMPS{r}, ϕ::GMPS{r}) where {r}
    # Checks 
    if !issimilar(ψ, ϕ)
        throw(ArgumentError("GMPS have differing properties."))
    end
    movecenter!(ψ, 1)
    movecenter!(ϕ, 1)
    return _add_exact(ψ, ϕ, false; cutoff=_TeNe_cutoff)
end

function -(ψ::GMPS{r}, ϕ::GMPS{r}) where {r}
    # Checks 
    if !issimilar(ψ, ϕ)
        throw(ArgumentError("GMPS have differing properties."))
    end
    movecenter!(ψ, 1)
    movecenter!(ϕ, 1)
    return _add_exact(ψ, ϕ, true; cutoff=_TeNe_cutoff)
end

function _add_exact(ψ::GMPS{r}, ϕ::GMPS{r}, subtract::Bool=false; kwargs...) where {r}
    Ψ = GMPS(r, dim(ψ), length(ψ); T=Base.promote_op(*, eltype(ψ), eltype(ϕ)))
    Ψ.center = 1
    for i in eachindex(Ψ)
        # Create the tensor 
        dim1 = i == 1 ? 1 : size(ψ[i], 1)+size(ϕ[i], 1)
        dim2 = i == length(Ψ) ? 1 : size(ψ[i], 2+r)+size(ϕ[i], 2+r)
        dims_phys = map(j->size(ψ[i], 1+j), Base.OneTo(rank(ψ)))
        tensor = zeros(eltype(Ψ), dim1, dims_phys..., dim2)
        
        # Find the dimensions of the tensors 
        dims_phys = map(j->1:size(tensor, 1+j), Base.OneTo(rank(ψ)))
        dims1 = (i == 1 ? (1:1) : 1:size(ψ[i], 1), dims_phys..., 
                 i == length(ψ) ? (1:1) : 1:size(ψ[i], 2+r))
        dims2 = (i == 1 ? (1:1) : size(ψ[i], 1)+1:size(tensor, 1), dims_phys..., 
                 i == length(ψ) ? (1:1) : size(ψ[i], 2+r)+1:size(tensor, 2+r))

        # Create the tensor
        tensor[dims1...] .= ψ[i]
        tensor[dims2...] .= (subtract && i == 1) ? -1*ϕ[i] : ϕ[i]
        Ψ[i] = tensor
        movecenter!(Ψ, i; kwargs...)
    end
    movecenter!(Ψ, 1; kwargs...)
    return Ψ
end



### Adjusting the bond dimensions 
export truncate!, expand!

"""
    truncate!(ψ::GMPS; kwargs...)

Truncate the bond dimension of a GMPS `ψ`.

# Optional Keyword Arguments
    
    - `cutoff::Float64=0.0`: Truncation criteria to reduce the bond dimension.
    Good values range from 1e-8 to 1e-14.
    - `mindim::Int=1`: Mininum dimension for truncated.
    - `maxdim::Int=0`: Maximum bond dimension for truncation. Set to 0 to have
    no limit.
"""
function truncate!(ψ::GMPS; kwargs...)
    if ψ.center != 1 && ψ.center != length(ψ)
        movecenter!(ψ, 1)
    end
    ctr = center(ψ) == 1 ? length(ψ) : 1
    movecenter!(ψ, ctr; kwargs...)
end

"""
    expand!(ψ::GMPS, bonddim::Int, noise=0.0)

Increase the bond dimension of a GMPS `ψ` to `bonddim`. 
Optionally, make the new parameters noisy, with strength `noise`.
"""
function expand!(ψ::GMPS, bonddim::Int, noise=0)
    movecenter!(ψ, 1)
    for i = 1:length(ψ)
        D1 = i == 1 ? 1 : bonddim 
        D2 = i == length(ψ) ? 1 : bonddim
        tensor = noise .* randn(eltype(ψ[i]), D1, map(j -> size(ψ[i], j), 2:ndims(ψ[i])-1)..., D2)
        tensor[map(j->Base.OneTo(size(ψ[i], j)), Base.OneTo(ndims(ψ[i])))...] .= ψ[i]
        ψ[i] = tensor
        if i > 1
            _moveright!(ψ, i-1)
        end
    end
    ψ.center = length(ψ)
    movecenter!(ψ, 1)
end


### Initialising GMPS 
export randomgmps
function GMPS(rank::Int, d::Int, N::Int; T::Type=ComplexF64)
    tensors = [zeros(T, (1, [d for __ = Base.OneTo(rank)]..., 1)) for _ = Base.OneTo(N)]
    return GMPS{rank}(d, tensors, 0)
end

"""
    randomgmps(rank::Int, dim::Int, length::Int, bonddim::Int; kwargs...)

Create a GMPS with random entries.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensors.
"""
function randomgmps(rank::Int, dim::Int, length::Int, bonddim::Int; T::Type=ComplexF64)
    # Create the GMPS
    ψ = GMPS(rank, dim, length)
    for i = 1:length
        D1 = i == 1 ? 1 : bonddim
        D2 = i == length ? 1 : bonddim
        idxs = (D1, ntuple(i->dim, rank)..., D2)
        ψ[i] = randn(T, idxs)
        if i > 1
            _moveright!(ψ, i-1)
            ψ.tensors[i] ./= norm(ψ[i])
        end
    end
    ψ.center = length
    movecenter!(ψ, 1)
    return ψ
end

### Information
function Base.show(io::IO, ψ::GMPS)
    println(io, "$(typeof(ψ))")
    for i = Base.OneTo(Base.length(ψ))
        println(io, "[$(i)] $(Base.size(ψ[i]))")
    end
end

### Creating copies
Base.copy(ψ::GMPS) = typeof(ψ)(dim(ψ), ψ.tensors, center(ψ))
Base.deepcopy(ψ::GMPS) = typeof(ψ)(Base.copy(dim(ψ)), Base.copy(ψ.tensors),
                                        Base.copy(center(ψ)))


### Save and write
function HDF5.write(parent::Union{HDF5.File, HDF5.Group}, name::AbstractString,
                    M::GMPS{r}) where {r}
    g = create_group(parent, name)
    attributes(g)["type"] = "MPS"
    attributes(g)["version"] = 1
    write(g, "length", length(M))
    write(g, "center", center(M))
    write(g, "rank", r)
    for i = 1:length(M)
        write(g, "MPS[$(i)]", M[i])
    end
end


function HDF5.read(parent::Union{HDF5.File, HDF5.Group}, name::AbstractString,
                    ::Type{GMPS})
    g = open_group(parent, name)
    if read(attributes(g)["type"]) != "MPS"
        error("HDF5 group of file does not contain MPS data.")
    end
    N = read(g, "length")
    center = read(g, "center")
    tensors = [read(g, "MPS[$(i)]") for i=Base.OneTo(N)]
    rank = read(g, "rank")
    return GMPS{rank}(size(tensors[1])[2], tensors, center)
end

### Conjugation of GMPS 
export conj, isconj
struct ConjGMPS{r} <: GMPSTrait where {r}
    MPS::GMPS{r}
end
TeNe.conj(ψ::GMPS) = ConjGMPS(ψ)
TeNe.conj(ψ::ConjGMPS) = ψ.MPS
isconj(ψ::ConjGMPS) = true