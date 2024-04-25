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
    randommpo(dim::Int, length::Int, bonddim::Int; kwargs...)

Create an MPO with dimensions `dim`, size `length` and bond dimension `bonddim`,
with random tensors.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensors.
"""
function randommpo(dim::Int, length::Int, bonddim::Int; kwargs...)
    return randomgmps(2, dim, length, bonddim; kwargs...)
end


"""
    productmpo(N::Int, A::AbstractArray; kwargs...)

Create a product MPO of size `N`, composed of array `A`. 
`A` can be a matrix or rank-4 tensor.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensors.
"""
function productmpo(N::Int, A::AbstractArray; T::Type=ComplexF64)
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
        for i = Base.range(1+firstindex(O), lastindex(O))
            O[i] = copy(A)
        end
        O[lastindex(O)] = A[:, :, :, end:end]
    else
        throw(ArgumentError("You must provide an array with just two or four dimensions."))
    end
    return O
end


### Applying an MPO 
# TODO: add densitymatrix method and variational method...
export applympo, applympo!

"""
    applympo(O::MPO, ψ::MPS; kwargs...)
    applympo(ψ::MPS, O::MPO; kwargs...)
    *(O::MPO, ψ::MPS)
    *(ψ::MPS, O::MPO)
    

Apply MPO `O` to MPS `ψ`.

# Optional Keyword Arguments

    - `alg=:naive`: The algorithm to carry out the multiplication. Use :naive 
      for an exact application, followed by SVD truncation. Use :zipup for quicker
      but less precise applications of the MPO. Use :densitymatrix for a one-shot method 
      with high accuracy. Use :variational for the slowest, but must optimal application.
    - `cutoff::Float64=0.0`: Truncation criteria to reduce the bond dimension.
      Good values range from 1e-8 to 1e-14.
    - `mindim::Int=1`: Minimum dimension for the truncation.
    - `maxdim::Int=0`: Maximum bond dimension for truncation. Set to 0 to have
      no limit.
"""
function applympo(O::MPO, ψ::MPS; alg=:naive, kwargs...)
    if !issimilar(O, ψ)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    ϕ = MPS(dim(ψ), length(ψ); T=Base.promote_op(*, eltype(O), eltype(ψ)))
    if alg==:naive
        _mpo_mps_naive!(ϕ, O, ψ; kwargs...)
    elseif alg==:zipup 
        _mpo_mps_zipup!(ϕ, O, ψ; kwargs...)
    else
        throw(ArgumentError("The algorithm $(alg) is unknown."))
    end
    return ϕ
end
applympo(ψ::MPS, O::MPO; kwargs...) = applympo(transpose(O), ψ; kwargs...)
*(O::MPO, ψ::MPS) = applympo(O, ψ; cutoff=1e-12)
*(ψ::MPS, O::MPO) = applympo(ψ, O; cutoff=1e-12)


"""
    applympo(O1::MPO, O2::MPO; kwargs...)

Apply MPO `O1` to MPO `O2`.

# Optional Keyword Arguments

    - `alg=:naive`: The algorithm to carry out the multiplication. Use :naive 
      for an exact application, followed by SVD truncation. Use :zipup for quicker
      but less precise applications of the MPO. Use :densitymatrix for a one-shot method 
      with high accuracy. Use :variational for the slowest, but must optimal application.
    - `cutoff::Float64=0.0`: Truncation criteria to reduce the bond dimension.
      Good values range from 1e-8 to 1e-14.
    - `mindim::Int=1`: Minimum dimension for the truncation.
    - `maxdim::Int=0`: Maximum bond dimension for truncation. Set to 0 to have
      no limit.
"""
function applympo(O1::MPO, O2::MPO; alg=:naive, kwargs...)
    if !issimilar(O1, O2)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    O = MPO(dim(O1), length(O2); T=Base.promote_op(*, eltype(O1), eltype(O2)))
    if alg==:naive
        _mpo_mpo_naive!(O, O1, O2; kwargs...)
    elseif alg==:zipup 
        _mpo_mpo_zipup!(O, O1, O2; kwargs...)
    else
        throw(ArgumentError("The algorithm $(alg) is unknown."))
    end
    return O
end
*(O1::MPO, O2::MPO) = applympo(O1, O2; cutoff=_TeNe_cutoff)

# Naive method; do the contraction exactly and then truncate
function _mpo_mps_naive!(ϕ::MPS, O::MPO, ψ::MPS; kwargs...)
    # Move canonical centre of both to the first site 
    movecenter!(O, 1)
    movecenter!(ψ, 1)
    
    # Type info...
    conjO = isconj(O)
    conjψ = isconj(ψ)
    transO = istranspose(O)

    # Do the contraction 
    ϕ.center = firstindex(ϕ)
    for i in eachindex(ϕ)
        tensor = contract(O[i], ψ[i], transO ? 2 : 3, 2, conjO, conjψ)
        tensor = _permutedims(tensor, (1, 4, 2, 3, 5))
        tensor = reshape(tensor, (size(tensor, 1)*size(tensor, 2), size(tensor, 3), size(tensor, 4)*size(tensor, 5)))
        ϕ[i] = copy(tensor)
        movecenter!(ϕ, i; kwargs...)
    end
    movecenter!(ϕ, firstindex(ϕ); kwargs...)
end

function _mpo_mpo_naive!(O::MPO, O1::MPO, O2::MPO; kwargs...)
    # Move canonical centre of both to the first site 
    movecenter!(O1, 1)
    movecenter!(O2, 1)
    
    # Type info...
    conjO1 = isconj(O1)
    conjO2 = isconj(O2)
    transO1 = istranspose(O1)
    transO2 = istranspose(O2)

    # Do the contraction 
    O.center = firstindex(O)
    for i in eachindex(O)
        tensor = contract(O1[i], O2[i], transO1 ? 2 : 3, transO2 ? 3 : 2, conjO1, conjO2)
        tensor = _permutedims(tensor, (1, 4, 2, 5, 3, 6))
        tensor = reshape(tensor, (size(tensor, 1)*size(tensor, 2), size(tensor, 3), size(tensor, 4),
                                  size(tensor, 5)*size(tensor, 6)))
        O[i] = copy(tensor)
        movecenter!(O, i; kwargs...)
    end
    movecenter!(O, firstindex(O); kwargs...)
end

# Zip-up method: see ``E.M. Stoudenmire, Steven R. White, New J. Phys. 12, 055026 (2010)``
function _mpo_mps_zipup!(ϕ::MPS, O::MPO, ψ::MPS; kwargs...)
    # Move canonical centre of both to the first site 
    movecenter!(O, firstindex(O))
    movecenter!(ψ, firstindex(ψ))

    # Type info...
    conjO = isconj(O)
    conjψ = isconj(ψ)
    transO = istranspose(O)

    # Do the contraction
    block = cache(eltype(ϕ), (1, size(O[begin], 1), size(ψ[begin], 1)), 2, 1) .= 1
    for i in eachindex(O)
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
    movecenter!(ϕ, firstindex(ϕ); kwargs...)
end

function _mpo_mpo_zipup!(O::MPO, O1::MPO, O2::MPO; kwargs...)
    # Move canonical centre of both to the first site 
    movecenter!(O1, firstindex(O1))
    movecenter!(O2, firstindex(O2))

    # Type info...
    conjO1 = isconj(O1)
    conjO2 = isconj(O2)
    transO1 = istranspose(O1)
    transO2 = istranspose(O2)

    # Do the contraction
    block = cache(eltype(O), (1, size(O1[begin], 1), size(O2[begin], 1)), 2, 1) .= 1
    for i in eachindex(O)
        block = contract(block, O1[i], 2, 1, false, conjO1)
        block = contract(block, O2[i], (2, transO1 ? 3 : 4), (1, transO2 ? 3 : 2), false, conjO2)
        if i < length(O)
            U, S, block = tsvd(block, (3, 5); kwargs...)
            block = contract(S, block, 2, 1)
            O[i] = U
        else
            block = reshape(block, (size(block, 1), size(block, 2), size(block, 4), 1))
            O[i] = copy(block)
        end
    end
    O.center = length(O)
    movecenter!(O, firstindex(O); kwargs...)
end

"""
    applympo(O::MPO, ψ::StateVector)
    applympo(ψ::StateVector, O::MPO)
    *(O::MPO, ψ::StateVector)
    *(ψ::StateVector, O::MPO)

Apply MPO `O` to StateVector `ψ`.

# Examples 

```julia-repl
julia> ψ = randomsv(2, 10);
julia> O = productmpo(10, [0 1; 1 0]);
julia> ϕ = O * ψ;
"""
function applympo(O::MPO, ψ::StateVector)
    ϕ = GStateTensor(1, dim(ψ), length(ψ))
    applympo!(ϕ, O, ψ)
    return ϕ
end
applympo(ψ::StateVector, O::MPO) = applympo(transpose(O), ψ)
*(O::MPO, ψ::StateVector) = applympo(O, ψ)
*(ψ::StateVector, O::MPO) = applympo(ψ, O)


"""
    applympo!(ϕ::StateVector, O::MPO, ψ::StateVector)
    applympo!(ϕ::StateVector, ψ::StateVector, O::MPO)

Apply MPO `O` to StateVector `ψ`, and store the result in `ϕ`.

# Examples 

```julia-repl
julia> ψ = randomsv(2, 10);
julia> O = productmpo(10, [0 1; 1 0]);
julia> ϕ = StateVector(2, 10);
julia> applympo!(ϕ, O, ψ);
"""
function applympo!(ϕ::StateVector, O::MPO, ψ::StateVector)
    if !issimilar(ϕ, O, ψ)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    # Apply the first tensor 
    ten = reshape(tensor(ψ), (size(tensor(ψ))..., 1))
    ten = contract(ten, O[1], (ndims(ten), 1), (1, istranspose(O) ? 2 : 3), isconj(ψ), isconj(O))

    # Loop through the remaining tensors 
    for i in range(firstindex(O)+1, lastindex(O))
        ten = contract(ten, O[i], (ndims(ten), 1), (1, istranspose(O) ? 2 : 3), false, isconj(O))
    end
    ten = reshape(ten, size(ten)[begin:end-1])
    tensor(ϕ) .= isconj(ϕ) ? conj.(ten) : ten
end
applympo!(ϕ::StateVector, ψ::StateVector, O::MPO) = applyMPO(ϕ, O, ψ)


### Inner products 
"""
    inner(ψ::MPS, O::MPO, ϕ::MPS)
    inner(ψ::MPS, O1::MPO, O2::MPO, ϕ::MPS)

Calculate the expectation of a string of operators `Os` with respect to MPSs `ψ` and `ϕ`.
"""
function inner(ψ::MPS, ϕs::Union{MPS, MPO}...)
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
    for i in eachindex(ψ)
        block = contract(block, ψ[i], 1, 1, false, conjψ)
        for j in eachindex(Os)
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
    for i in eachindex(Os[begin])
        # Contract with the first MPO in the term
        block = contract(block, Os[1][i], 1, 1, false, conjOs[1])
        if transOs[1]
            block = permutedim(block, ndims(block)-1, ndims(block)-2)
        end

        # Contract with the central MPOs
        for j in range(firstindex(Os)+1, lastindex(Os)-1)
            block = contract(block, Os[j][i], (1, ndims(block)-1),
                            (1, transOs[j] ? 3 : 2), false, conjOs[j])
        end

        # Contract with final MPO 
        block = contract(block, Os[end][i], (1, ndims(block)-1, 2),
                         (1, transOs[end] ? 3 : 2, transOs[end] ? 2 : 3), false, conjOs[end])
    end

    return block[]
end