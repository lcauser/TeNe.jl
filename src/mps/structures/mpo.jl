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
"""
    inner(ψ::MPS, Os::MPO..., ϕ::MPS)

Calculate the expectation of a string of operators `Os` with respect to MPSs `ψ` and `ϕ`.
"""
function inner(ψ::MPS, ϕs::Union{GMPS, GMPSTrait}...)
    # Checks 
    if !issimilar(ψ, ϕs...)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    for i = 1:length(ϕs)-1
        if !ismpo(ϕs[i])
            throw(ArgumentError("The inner terms in the braket must be MPOs."))
        end
    end
    if !ismps(ϕs[end])
        throw(ArgumentError("The last term in the MPS must be an MPO."))
    end
    return _mps_mpo_mps_product(ψ, ϕs[end], ϕs[1:end-1]...)
end

function _mps_mpo_mps_product(ψ::MPS, ϕ::MPS, Os::MPO...)
    # Type info...
    T = Base.promote_op(*, eltype(ψ), eltype(ϕ), eltype.(Os)...)
    conjψ = !isconj(ψ) # Not because inner product has conj on bra by default...
    conjϕ = isconj(ϕ)
    conjOs = map(O->isconj(O), Os)
    transOs = map(O->istranspose(O), Os)

    # Contraction 
    dims = (size(ψ[begin], 1), map(O->size(O[begin], 1), Os)..., size(ϕ[begin], 1))
    block = cache(T, dims, 2, 1) .= 1
    for i = 1:length(ψ)
        block = contract(block, ψ[i], 1, 1, false, conjψ)
        for j = 1:length(Os)
            block = contract(block, Os[j][i], (1, ndims(block)-1), (1, transOs[j] ? 3 : 2), false, conjOs[j])
        end
        block = contract(block, ϕ[i], (1, ndims(block)-1), (1, 2), false, conjϕ)
    end

    return block[]
end