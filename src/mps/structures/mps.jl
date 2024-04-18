#=
    Matrix product states are typically used to represent wavefunctions of 1D
    quantum many-body systems, but can generally be used to describe N discrete
    degrees of freedom.
=#

const MPS = Union{GMPS{1}, ConjGMPS{1}}
export MPS 

export ismps
"""
    ismps(ψ)

Check to see if an object is an MPS.
"""
function ismps(ψ)
    return typeof(ψ) <: MPS
end

### Initalising MPSs 
export randommps, productmps

"""
    MPS(dim::Int, length::Int)

Create an MPS with physical dimension `dim` and `length` sites.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensors.
"""
function MPS(dim::Int, length::Int; kwargs...)
    return GMPS(1, dim, length; kwargs...)
end

"""
    randommps(dim::Int, length::Int, bonddim::Int)

Create an MPS with dimensions `dim`, size `length` and bond dimension `bonddim`,
with random tensors.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensors.
"""
function randommps(dim::Int, length::Int, bonddim::Int; kwargs...)
    return randomgmps(1, dim, length, bonddim; kwargs...)
end

"""
    productmps(N::Int, A<:AbstractArray; kwargs...)

Create a product MPS of size `N`, composed of array `A`. 
`A` can be a vector or rank-3 tensor.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensors.
    - `normalise::Bool=false`: Normalise the MPS after creating it?
"""
function productmps(N::Int, A::Q; T::Type=ComplexF64, normalize::Bool=false) where {Q<:AbstractArray}
    if ndims(A) == 1
        ψ = MPS(length(A), N; T=T)
        for i = Base.OneTo(N)
            ψ[i] = reshape(copy(A), (1, length(A), 1))
        end
    elseif ndims(A) == 3
        ψ = MPS(size(A, 2), N; T=T)
        ψ[1] = A[1:1, :, :]
        for i = Base.range(2, N)
            ψ[i] = copy(A)
        end
        ψ[N] = A[:, :, end:end]
    else
        throw(ArgumentError("You must provide an array with just one or three dimensions."))
    end
    if normalize 
        movecenter!(ψ, 1)
        normalize!(ψ)
    end
    return ψ
end


### Inner products 
export inner, dot

"""
    inner(ψ::MPS, ϕ::MPS)

Calculate the inner product of two MPSs `ψ` and `ϕ`.
"""
function inner(ψ::MPS, ϕ::MPS)
    # Checks 
    if !issimilar(ψ, ϕ)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    return _mps_mps_product(ψ, ϕ)
end

"""
    dot(ψ::MPS, ϕ::MPS)

Calculate the inner product of two MPSs `ψ` and `ϕ`.
"""
dot(ψ::MPS, ϕ::MPS) = inner(ψ, ϕ)

import Base.*
"""
    *(ψ::MPS, ϕ::MPS)

Calculate the inner product of two MPSs `ψ` and `ϕ`.
"""
*(ψ::MPS, ϕ::MPS) = inner(ψ, ϕ)


function _mps_mps_product(ψ::MPS, ϕ::MPS)
    # Type info...
    T = Base.promote_op(*, eltype(ψ), eltype(ϕ))
    conjψ = !isconj(ψ) # Not because inner product has conj on bra by default...
    conjϕ = isconj(ϕ)

    # Contract the network...
    dims_prev = (size(ψ[begin], 1), size(ϕ[begin], 1))
    block = cache(T, dims_prev, 2, 1)
    block .= 1
    for i = 1:length(ψ)
        # Caching...
        dims1 = (dims_prev[2], size(ψ[i], 2), size(ψ[i], 3))
        dims2 = (size(ψ[i], 3), size(ϕ[i], 3))
        sub_level = (prod(dims_prev) == prod(dims1)) || (prod(dims1) == prod(dims2)) ? 2 : 1

        # Contract the new block 
        block_new = contract(block, ψ[i], 1, 1, false, conjψ; sublevel=sub_level)
        block = contract(block_new, ϕ[i], (1, 2), (1, 2), false, conjϕ)
        dims_prev = dims2
    end
    return block[]
end