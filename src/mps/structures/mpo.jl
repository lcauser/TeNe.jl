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


### Computing traces of MPO (products)
export trace

"""
    trace(Os::MPO...)

Compute the trace of the product of many MPOs.
"""
function TeNe.trace(Os::MPO...)
    # Checks 
    issimilar(Os...)

    if length(Os) == 1
        return _mpo_trace(Os...)
    else
        return _mpo_mpo_trace(Os...)
    end
end

function _mpo_trace(O::MPO)
     # Type info...
     T = eltype(O)
     conjO = isconj(O)

     # Do the contraction 
    block = cache(T, size(O[begin], 1), 2, 1) .= 1
    for i in eachindex(O)
        block = contract(block, trace(O[i], 2, 3), 1, 1, false, conjO)
    end
    return block[]
end

function _mpo_mpo_trace(Os::MPO...)
    # Type info...
    T = Base.promote_op(*, eltype.(Os)...)
    conjOs = map(O->isconj(O), Os)
    transOs = map(O->istranspose(O), Os)

    # Do the contraction 
    block = cache(T, map(O->size(O[begin], 1), Os), 2, 1) .= 1
    for i in eachindex(Os[1])
        # Contract with the first MPO in the term
        block = contract(block, Os[1][i], 1, 1, false, conjOs[1])
        if transOs[1]
            block = permutedim(block, ndims(block)-1, ndims(block)-2)
        end

        # Contract with the central MPOs
        for j = 2:length(Os)-1
            block = contract(block, Os[j][i], (1, ndims(block)-1),
                            (1, transOs[j] ? 3 : 2), false, conjOs[j])
        end

        # Contract with final MPO 
        block = contract(block, Os[end][i], (1, ndims(block)-1, 2),
                         (1, transOs[end] ? 3 : 2, transOs[end] ? 2 : 3), false, conjOs[end])
    end

    return block[]
end

### Applying an MPO 

function _mpo_mps_zipup!(ϕ::MPS, O::MPO, ψ::MPS; kwargs...)
    # Move canonical centre of both to the first site 
    movecenter!(O, 1)
    movecenter!(ψ, 1)

    # Type info...
    T = Base.promote_op(*, eltype(O), eltype(ψ))
    conjO = isconj(O)
    conjψ = isconj(ψ)
    transO = istranspose(O)

    # Do the contraction
    block = cache(T, (1, size(O[begin], 1), size(ψ[begin], 1)), 2, 1) .= 1
    for i = 1:length(O)
        block = contract(block, O[i], 2, 1, false, conjO)
        block = contract(block, ψ[i], (2, transO ? 3 : 4), (1, 2), false, conjψ)
        if i < length(O)
            U, S, block = tsvd(block, (3, 4); kwargs...)
            block = contract(S, block, 2, 1)
            ϕ[i] = U
        else
            block = reshape(block, (size(block, 1), size(block, 2), 1))
            ϕ[i] = copy(block)
        end
    end
    ϕ.center = length(ϕ)
    movecenter!(ϕ, 1; kwargs...)
end