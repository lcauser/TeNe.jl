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

### MPS Properties 
export rank, dim, center, bonddim, maxbonddim

"""
    Base.eltype(::GMPS)

Returns the type of parameters within a GMPS.
"""
Base.eltype(ψ::GMPS) = _promote_tensor_eltype(ψ.tensors...)

"""
    rank(::GMPS)

Returns the rank of an MPS object.
"""
TeNe.rank(::GMPS{r}) where {r} = r


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

### Physical dimenons

"""
    dim(ψ::GMPS, [which::Int, site::Int])

The size of the physical dimensions in a GMPS. Returns `0` for heterogeneous 
systems (i.e. an invariant physical dimension). The axis and the lattice site
can also be specified.
"""
function dim(ψ::GMPS)
    ds = Tuple(map(j->dim(ψ, j), Base.OneTo(rank(ψ))))
    if all(map(j->j==ds[1], ds))
        return ds[1]
    else
        return 0
    end
end
function dim(ψ::GMPS, which::Int)
    ds = dims(ψ, which)
    if all(map(j->j==ds[1], ds))
        return ds[1]
    else
        return 0
    end
end
dim(ψ::GMPS, which::Int, site::Int) = size(ψ[site], 1+which)

"""
    dim(ψ::GMPS, which::Int)

The size of the physical dimensions across an axis in a GMPS.
"""
function dims(ψ::GMPS, which::Int)
    return Tuple(map(j->dim(ψ, which, j), eachindex(ψ)))
end

### Bond dimensions
"""
    bonddim(::GMPS, idx::Int)

Return the bond dimension size between idx and idx + 1. Returns nothing if
out of range.
"""
function bonddim(ψ::GMPS, site::Int)
    (site < firstindex(ψ)-1 || site > lastindex(ψ)) && return nothing
    site == firstindex(ψ)-1 && return 1
    return size(ψ[site])[2+rank(ψ)]
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
        movecenter!(ψ, firstindex(ψ))
    end
    return LinearAlgebra.norm(ψ[center(ψ)])
end

"""
    normalize!(ψ::GMPS)

Normalize a GMPS.
"""
function TeNe.normalize!(ψ::GMPS)
    if center(ψ) == 0
        movecenter!(ψ, firstindex(ψ))
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

Move the orthogonal center of an GMPS `ψ` to `idx`.

# Optional Keyword Arguments

    - `cutoff::Float64=0.0`: Truncation criteria to reduce the bond dimension.
      Good values range from 1e-8 to 1e-14.
    - `mindim::Int=1`: Minimum dimension for the truncation.
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
        for i = Base.range(firstindex(ψ), idx-1)
            _moveright!(ψ, i; kwargs...)
        end
        for i = Base.range(lastindex(ψ), idx+1, step=-1)
            _moveleft!(ψ, i; kwargs...)
        end
    else
        if idx > center(ψ)
            for i = Base.range(center(ψ), idx-1)
                _moveright!(ψ, i; kwargs...)
            end
        elseif idx < center(ψ)
            for i = Base.range(center(ψ), idx+1, step=-1)
                _moveleft!(ψ, i; kwargs...)
            end
        end
    end
    ψ.center = idx
end

### Replace the sites within a MPS with a contraction of multiple sites
"""
    replacesites!(ψ::GMPS, A::AbstractArray, site::Int, direction::Bool=false, normalize::Bool=false; kwargs...)

Replace the tensors over a small range of sites in a GMPS.

# Arguments 

    - `ψ::GMPS`: The GMPS.
    - `A::AbstractArray': The updated composite tensor.
    - `site::Int`: The first site in the range which is replaced.
    - `direction::Bool=false`: The sweeping direction; `true` is sweeping left, `false` sweeping right.

# Optional Keyword Arguments

    - `normalize::Bool=false`: Normalize the GMPS after replacing the sites?
    - `cutoff::Float64=0.0`: Truncation criteria to reduce the bond dimension.
      Good values range from 1e-8 to 1e-14.
    - `mindim::Int=1`: Mininum dimension for truncated.
    - `maxdim::Int=0`: Maximum bond dimension for truncation. Set to 0 to have
      no limit.
"""
function replacesites!(ψ::GMPS, A::AbstractArray, site::Int, direction::Bool=false;
        normalize::Bool=false, kwargs...)
    # Determine the number of sites
    nsites = fld((ndims(A) - 2), rank(ψ))

    if nsites == 1
        # Deal with case of just one site
        ψ[site] .= A
        if normalize
            normalize!(ψ)
        end
        return nothing
    else
        # Have to restore MPS form using SVDs
        U = A
        for i = Base.OneTo(nsites-1)
            if direction 
                # Sweeping left
                U, S, V = tsvd(U, Tuple(Base.range(ndims(U)-rank(ψ), ndims(U))); kwargs...)
                ψ[site - i + 1] = V
                U = contract(U, S, ndims(U), 1)
            else
                # Sweeping right 
                V, S, U = tsvd(U, Tuple(Base.range(2+rank(ψ), ndims(U))); kwargs...)
                ψ[site + i - 1] = V 
                U = contract(S, U, 2, 1) 
            end
        end
        lastsite = direction ? site + 1 - nsites : site - 1 + nsites 
        if size(ψ[lastsite]) == size(U)
            ψ[lastsite] .= U
        else
            ψ[lastsite] = copy(U)
        end
        ψ.center = lastsite
    end

    if normalize
        normalize!(ψ)
    end
end

### Products with numbers
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

function +(ψ::GMPS{r}, ϕ::GMPS{r}) where {r}
    _vec_vec_validation(ψ, ϕ)
    movecenter!(ψ, firstindex(ψ))
    movecenter!(ϕ, firstindex(ϕ))
    return _add_exact(ψ, ϕ, false; cutoff=_TeNe_cutoff)
end

function -(ψ::GMPS{r}, ϕ::GMPS{r}) where {r}
    _vec_vec_validation(ψ, ϕ)
    movecenter!(ψ, firstindex(ψ))
    movecenter!(ϕ, firstindex(ϕ))
    return _add_exact(ψ, ϕ, true; cutoff=_TeNe_cutoff)
end

function _add_exact(ψ::GMPS{r}, ϕ::GMPS{r}, subtract::Bool=false; kwargs...) where {r}
    Ψ = GMPS(r, dim(ψ), length(ψ); T=_promote_tensor_eltype(ψ, ϕ))
    Ψ.center = firstindex(Ψ)
    for i in eachindex(Ψ)
        # Create the tensor 
        dim1 = i == 1 ? 1 : size(ψ[i], 1)+size(ϕ[i], 1)
        dim2 = i == length(Ψ) ? 1 : size(ψ[i], 2+r)+size(ϕ[i], 2+r)
        dims_phys = map(j->size(ψ[i], 1+j), Base.OneTo(rank(ψ)))
        tensor = zeros(eltype(Ψ), dim1, dims_phys..., dim2)
        
        # Find the dimensions of the tensors 
        dims_phys = map(j->1:size(tensor, 1+j), Base.OneTo(rank(ψ)))
        dims1 = (i == 1 ? Base.OneTo(1) : axes(ψ[i], 1), dims_phys..., 
                 i == length(ψ) ? Base.OneTo(1) : axes(ψ[i], 2+r))
        dims2 = (i == 1 ? Base.OneTo(1) : range(size(ψ[i], 1)+1, size(tensor, 1)), dims_phys..., 
                 i == length(ψ) ? Base.OneTo(1) :  range(size(ψ[i], 2+r)+1, size(tensor, 2+r)))

        # Create the tensor
        tensor[dims1...] .= ψ[i]
        tensor[dims2...] .= (subtract && i == 1) ? -1*ϕ[i] : ϕ[i]
        Ψ[i] = tensor
        movecenter!(Ψ, i; kwargs...)
    end
    movecenter!(Ψ, firstindex(Ψ); kwargs...)
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
    - `mindim::Int=1`: Minimum dimension for the truncation.
    - `maxdim::Int=0`: Maximum bond dimension for truncation. Set to 0 to have
    no limit.
"""
function truncate!(ψ::GMPS; kwargs...)
    if ψ.center != firstindex(ψ) && ψ.center != lastindex(ψ)
        movecenter!(ψ, firstindex(ψ))
    end
    ctr = center(ψ) == firstindex(ψ) ? lastindex(ψ) : firstindex(ψ)
    movecenter!(ψ, ctr; kwargs...)
end

"""
    expand!(ψ::GMPS, bonddim::Int, noise=0.0)

Increase the bond dimension of a GMPS `ψ` to `bonddim`. 
Optionally, make the new parameters noisy, with strength `noise`.
"""
function expand!(ψ::GMPS, bonddim::Int, noise=0)
    movecenter!(ψ, firstindex(ψ))
    for i in eachindex(ψ)
        D1 = i == firstindex(ψ) ? 1 : bonddim 
        D2 = i == lastindex(ψ) ? 1 : bonddim
        tensor = noise .* randn(eltype(ψ[i]), D1, map(j -> size(ψ[i], j), 2:ndims(ψ[i])-1)..., D2)
        tensor[map(j->Base.OneTo(size(ψ[i], j)), Base.OneTo(ndims(ψ[i])))...] .= ψ[i]
        ψ[i] = tensor
        if i > firstindex(ψ)
            _moveright!(ψ, i-1)
        end
    end
    ψ.center = lastindex(ψ)
    movecenter!(ψ, firstindex(ψ))
end


### Initialising GMPS 
export randomgmps
function GMPS(rank::Int, d::Int, N::Int; T::Type=ComplexF64)
    tensors = [zeros(T, (1, [d for __ = Base.OneTo(rank)]..., 1)) for _ = Base.OneTo(N)]
    return GMPS{rank}(tensors, 0)
end

"""
    randomgmps(rank::Int, dim::Int, length::Int, bonddim::Int; kwargs...)

Create a GMPS with random entries.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensors.
"""
function randomgmps(rank::Int, dim::Int, length::Int, bonddim::Int; T::Type=ComplexF64)
    # Create the GMPS
    ψ = GMPS(rank, dim, length; T=T)
    for i = Base.OneTo(length)
        D1 = i == firstindex(ψ) ? 1 : bonddim
        D2 = i == lastindex(ψ) ? 1 : bonddim
        idxs = (D1, ntuple(i->dim, rank)..., D2)
        ψ[i] = randn(T, idxs)
        if i > firstindex(ψ)
            _moveright!(ψ, i-1)
            ψ.tensors[i] ./= norm(ψ[i])
        end
    end
    ψ.center = length
    movecenter!(ψ, firstindex(ψ))
    return ψ
end

"""
    GMPS(ψ::Union{GStateTensor, GStateTensorTrait}; kwargs...)

Write a GStateTensor as a GMPS.

# Optional Keyword Arguments
    
    - `cutoff::Float64=0.0`: Truncation criteria to reduce the bond dimension.
    Good values range from 1e-8 to 1e-14.
    - `mindim::Int=1`: Minimum dimension for the truncation.
    - `maxdim::Int=0`: Maximum bond dimension for truncation. Set to 0 to have
    no limit.

# Examples
```julia-repl
julia> ψ = randomsv(2, 10);
julia> ψ = GMPS(ψ; cutoff=1e-12);
```
"""
function GMPS(ψ::GStateTensor; kwargs...)
    ψmps = GMPS(rank(ψ), dim(ψ), length(ψ); T=eltype(ψ))
    ten = reshape(tensor(ψ), 1, size(tensor(ψ))..., 1)
    for i = Base.OneTo(length(ψ)-1)
        U, S, ten = tsvd(ten, Tuple(2+rank(ψ):ndims(ten)); kwargs...)
        ψmps[i] = U 
        ten = contract(S, ten, 2, 1)
    end
    ψmps[end] = ten
    ψmps.center = length(ψmps)
    return ψmps
end

### Information
function Base.show(io::IO, ψ::GMPS)
    println(io, "$(typeof(ψ))")
    for i = Base.OneTo(Base.length(ψ))
        println(io, "[$(i)] $(Base.size(ψ[i]))")
    end
end

### Creating copies
Base.copy(ψ::GMPS) = typeof(ψ)(ψ.tensors, center(ψ))
Base.deepcopy(ψ::GMPS) = typeof(ψ)(Base.deepcopy(ψ.tensors), Base.deepcopy(center(ψ)))


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
    return GMPS{rank}(tensors, center)
end

### Conjugation of GMPS 
export conj, isconj
struct ConjGMPS{r} <: GMPSTrait where {r}
    MPS::GMPS{r}
end
TeNe.conj(ψ::GMPS) = ConjGMPS(ψ)
TeNe.conj(ψ::ConjGMPS) = ψ.MPS
isconj(ψ::ConjGMPS) = true