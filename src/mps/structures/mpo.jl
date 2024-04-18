#=
    Matrix product operators are typically used to represent Hamiltonians or operators
    of 1D quantum many-body systems, but can more generally be used to describe a rank-2
    tensor of N discrete degrees of freedom.
=#

### Traits for MPOs 
export conj, transpose, adjoint, istranspose
# Transpose
struct TransposeMPO <: GMPSTrait
    MPS::GMPS{2}
end
TeNe.transpose(O::GMPS{2}) = TransposeMPO(O)
TeNe.transpose(O::TransposeMPO) = O.MPS
istranspose(ψ::TransposeMPO) = true

# Adjoint
struct AdjointMPO <: GMPSTrait
    MPS::GMPS{2}
end
TeNe.adjoint(O::GMPS{2}) = AdjointMPO(O)
TeNe.adjoint(O::AdjointMPO) = O.MPS
isconj(ψ::AdjointMPO) = true
istranspose(ψ::AdjointMPO) = true

# Trait rules
TeNe.conj(O::TransposeMPO) = AdjointMPO(O.MPS)
TeNe.conj(O::AdjointMPO) = TransposeMPO(O.MPS)
TeNe.transpose(O::ConjGMPS{2}) = AdjointMPO(O.MPS)
TeNe.transpose(O::AdjointMPO) = ConjGMPS(O.MPS)
TeNe.adjoint(O::ConjGMPS{2}) = TransposeMPO(O.MPS)
TeNe.adjoint(O::TransposeMPO) = ConjGMPS(O.MPS)


const MPO = Union{GMPS{2}, ConjGMPS{2}, TransposeMPO, AdjointMPO}
export MPO 


export ismpo
"""
    ismpo(O)

Check to see if an object is an MPO.
"""
function ismpo(O)
    return typeof(O) <: MPO
end


### Initalising MPOs 
export randommpo, productmpo
"""
    MPO(dim::Int, length::Int)

Create an MPO with physical dimension `dim` and `length` sites.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensors.
"""
function MPO(dim::Int, length::Int; kwargs...)
    return GMPS(2, dim, length; kwargs...)
end


"""
    randommpo(dim::Int, length::Int, bonddim::Int)

Create an MPO with dimensions `dim`, size `length` and bond dimension `bonddim`,
with random tensors.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensors.
"""
function randommpo(dim::Int, length::Int, bonddim::Int; kwargs...)
    return randomgmps(2, dim, length, bonddim; kwargs...)
end


"""
    productmpo(N::Int, A<:AbstractArray; kwargs...)

Create a product MPO of size `N`, composed of array `A`. 
`A` can be a matrix or rank-4 tensor.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensors.
"""
function productmpo(N::Int, A::Q; T::Type=ComplexF64) where {Q<:AbstractArray}
    if ndims(A) == 2
        if size(A, 1) != size(A, 2) 
            throw(ArgumentError("A must be a square matrix!"))
        end
        O = MPO(size(A, 1), N; T=T)
        for i = Base.OneTo(N)
            O[i] = reshape(copy(A), (1, size(A)..., 1))
        end
    elseif ndims(A) == 4
        if size(A, 2) != size(A, 3) 
            throw(ArgumentError("Array must have physical dimensions of the same length!"))
        end
        if size(A, 1) != size(A, 4) 
            throw(ArgumentError("Array must have equal bond dimensions!"))
        end
        O = MPO(size(A, 2), N; T=T)
        O[1] = A[1:1, :, :, :]
        for i = Base.range(2, N)
            O[i] = copy(A)
        end
        O[N] = A[:, :, :, end:end]
    else
        throw(ArgumentError("You must provide an array with just two or four dimensions."))
    end
    return O
end

### Inner products 
#=
function _mps_mpo_mps_product(ψ::MPS, Os::MPO..., ϕ::MPS)
    # Type info...
    T = Base.promote_op(*, eltype(ψ), eltype(ϕ), eltype.(Os)...)

    # Contraction 
    dims_prev = (size(ψ[begin], 1), map(O->size(O[begin], 1), Os)..., size(ϕ[begin], 1))
    sub_level_prev = 1
    block = cache(T, dims_prev, 2, sub_level_prev)
    block .= 1
    for i = 1:length(ψ)
        # Contract with ψ
        dims1 = (dims_prev[2:end]..., size(ψ[i], 2), size(ψ[i], 3))

        # Caching...
        dims1 = (dims_prev[2], size(ψ[i], 2), size(ψ[i], 3))
        dims2 = (size(ψ[i], 3), size(ϕ[i], 3))
        sub_level = (prod(dims_prev) == prod(dims1)) || (prod(dims1) == prod(dims2)) ? 2 : 1

        # Contract the new block 
        block_new = contract(block, ψ[i], 1, 1, false, conjψ; sublevel=sub_level)
        block = contract(block_new, ϕ[i], (1, 2), (1, 2), false, conjϕ)
        dims_prev = dims2
    end
end
=#