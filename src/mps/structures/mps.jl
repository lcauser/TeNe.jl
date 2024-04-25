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
    MPS(dim::Int, length::Int; kwargs...)

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
        for i = Base.range(1+firstindex(ψ), lastindex(ψ))
            ψ[i] = copy(A)
        end
        ψ[N] = A[:, :, end:end]
    else
        throw(ArgumentError("You must provide an array with just one or three dimensions."))
    end
    if normalize 
        movecenter!(ψ, firstindex(ψ))
        normalize!(ψ)
    end
    return ψ
end


### Entanglement entropy 
export entropy
"""
    entropy(ψ::MPS, site::Int)

Calculate the entanglement entropy of an MPS `ψ` between sites `i` and `i+1`.
"""
function entropy(ψ::MPS, site::Int)
    movecenter!(ψ, site) # Move center to the site where entropy is calculated
    _, S, _ = tsvd(ψ[site], 3)
    S ./= norm(ψ)
    S2 = diag(S) .^ 2
    return -1 * sum(S2 .* log.(S2))
end


### Sampling a configuration from an MPS 
export sample
"""
    sample(ψ::MPS)

Sample the MPS `ψ` with the interpretation that it is a wavefunction (or Born ansatz).
"""
function sample(ψ::MPS)
    movecenter!(ψ, firstindex(ψ)) # move centre to begin for reduced sampling cost
    config = zeros(Int, length(ψ))
    block = cache(eltype(ψ), (size(ψ[begin], 1)), 2, 1) .= 1
    block = cache_ones((size(ψ[begin], 1)))
    P = 1.0
    for i in eachindex(ψ)
        block = contract(block, ψ[i], 1, 1)
        probs = contract(block, block, 2, 2, true, false)
        cum_prob = 0.0
        cum_probs = zeros(Float64, size(ψ[i], 2))
        for j in axes(ψ[i], 2)
            cum_prob += real(probs[j, j])
            cum_probs[j] = cum_prob
        end
        config[i] = findfirst(cum_probs .> cum_prob*rand())
        block = cache(size(block, 2), block, probs) .= (@view block[config[i], :])
        P *= real(probs[config[i], config[i]] / cum_prob)
    end
    return config, P
end
