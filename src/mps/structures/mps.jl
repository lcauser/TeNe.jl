#=
    Matrix product states are typically used to represent wavefunctions of 1D
    quantum many-body systems, but can generally be used to describe N discrete
    degrees of freedom.
=#

const MPS = Union{GMPS{1},ConjGMPS{1}}
export MPS

export ismps
"""
    ismps(ψ)

Check to see if an object is an MPS.
"""
function ismps(ψ)
    return typeof(ψ) <: MPS
end

export dim, dims
dim(ψ::MPS, site::Int) = size(ψ[site], 2)
function dim(ψ::MPS)
    ds = dims(ψ)
    if all(map(j->j==ds[1], ds))
        return ds[1]
    else
        return 0
    end
end
dims(ψ::MPS) = dims(ψ, 1)

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
function productmps(
    N::Int,
    A::Q;
    T::Type = ComplexF64,
    normalize::Bool = false,
) where {Q<:AbstractArray}
    if ndims(A) == 1
        ψ = MPS(length(A), N; T = T)
        for i in Base.OneTo(N)
            ψ[i] = reshape(copy(A), (1, length(A), 1))
        end
    elseif ndims(A) == 3
        ψ = MPS(size(A, 2), N; T = T)
        ψ[1] = A[1:1, :, :]
        for i in Base.range(1+firstindex(ψ), lastindex(ψ))
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

"""
    productmps(lt::LatticeTypes, states::AbstractVector{String})

Create an MPS from a string of states.

# Example 

```julia-repl
julia> lt = Qubits();
julia> ψ = productmps(lt, ["up" for _ = 1:20]);
```
"""
function productmps(lt::LatticeTypes, states::AbstractVector{String})
    ψ = MPS(dim(lt), length(states); T = eltype(lt))
    for i in eachindex(states)
        ψ[i][1, :, 1] .= state(lt, states[i])
    end
    movecenter!(ψ, firstindex(ψ))
    return ψ
end

"""
    MPS(ψ::StateVector; kwargs...)

Write a StateVector as a MPS.

# Optional Keyword Arguments
    
    - `cutoff::Float64=0.0`: Truncation criteria to reduce the bond dimension.
    Good values range from 1e-8 to 1e-14.
    - `mindim::Int=1`: Minimum dimension for the truncation.
    - `maxdim::Int=0`: Maximum bond dimension for truncation. Set to 0 to have
    no limit.

# Examples
```julia-repl
julia> ψ = randomsv(2, 10);
julia> ψ = MPS(ψ; cutoff=1e-12);
```
"""
function MPS(ψ::GStateTensor{1}; kwargs...)
    return GMPS(ψ; kwargs...)
end
MPS(ψ::ConjGStateTensor{1}; kwargs...) = conj(GMPS(ψ.StateTensor; kwargs...))

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


### Conjugate of MPS 
entropy(ψ::ConjGMPS{1}, site::Int) = entropy(ψ.MPS, site)
sample(ψ::ConjGMPS{1}) = sample(ψ.MPS)
