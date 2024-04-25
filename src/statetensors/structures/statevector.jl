#=
    State vectors can be used to describes states of many-body systems with discrete
    degress of freedom.
=#

const StateVector = Union{GStateTensor{1}, ConjGStateTensor{1}}
export StateVector 

"""
    isstatevector(ψ)

Check to see if an object is a state vector.
"""
function isstatevector(ψ)
    return typeof(ψ) <: StateVector
end

### Initialising state vectors 
export randomsv, randomstatevector, productsv, productstatevector
"""
    StateVector(dim::Int, length::Int; kwargs...)

Create a state vector with physical dimension `dim` and `length` sites.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensor.
"""
function StateVector(dim::Int, length::Int; kwargs...)
    return GStateTensor(1, dim, length; kwargs...)
end

"""
    randomsv(dim::Int, length::Int)
    randomstatevector(dim::Int, length::Int)

Create a random state vector with dimensions `dim`, size `length`.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensors.
"""
function randomsv(dim::Int, length::Int; kwargs...)
    return randomgst(1, dim, length, kwargs...)
end
randomstatevector(d::Int, N::Int; kwargs...) = randomsv(d, N; kwargs...)


"""
    productsv(dim::Int, A::AbstractVector)
    productstatevector(dim::Int, A::AbstractVector)

Create a product state vector with length `length` and local state `A`.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensors.
"""
function productsv(N::Int, A::AbstractVector; kwargs...)
    return productgst(N, A; kwargs...)
end
productstatevector(N::Int, A::AbstractVector; kwargs...) = productsv(N, A; kwargs...)


### Inner products 
export inner, dot
"""
    inner(ψ::StateVector, ϕ::StateVector)
    dot(ψ::StateVector, ϕ::StateVector)
    *(ψ::StateVector, ϕ::StateVector)

Calculate the inner product of two StateVectors `ψ` and `ϕ`.
"""
function inner(ψ::StateVector, ϕ::StateVector)
    # Checks 
    if length(ψ) != length(ϕ) || !_sv_sv_product_checkdims(ψ, ϕ)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    return _sv_sv_product(ψ, ϕ)
end
dot(ψ::StateVector, ϕ::StateVector) = inner(ψ, ϕ)
import Base.*
*(ψ::StateVector, ϕ::StateVector) = inner(ψ, ϕ)

function _sv_sv_product(ψ::StateVector, ϕ::StateVector)
    return contract(reshape(tensor(ψ), length(tensor(ψ))),
             reshape(tensor(ϕ), length(tensor(ϕ))),
             1, 1, !isconj(ψ), isconj(ϕ))[]
end

function _sv_sv_product_checkdims(ψ::StateVector, ϕ::StateVector)
    for i = Base.OneTo(length(ψ))
        if size(tensor(ψ), i) != size(tensor(ϕ), i)
            return false
        end
    end
    return true
end

### Entanglement entropy 
export entropy
"""
    entropy(ψ::StateVector, site::Int)

Calcualte the bipartition entanglement entropy between sites `site` and `site+1`
for a StateVector `ψ`.
"""
function entropy(ψ::StateVector, site::Int)
    # Checks 
    if site < 0 || site >= length(ψ)
        throw(ArgumentError("site=$(site) is out of range."))
    end

    # Dimensions 
    inner_dims = Tuple(map(j->size(tensor(ψ), j), Base.OneTo(site))) # Allocations?
    outer_dims = Tuple(map(j->size(tensor(ψ), j), Base.range(site+1, length(ψ)))) # Allocations?

    # SVD
    _, S, _ = tsvd(reshape(tensor(ψ), prod(inner_dims), prod(outer_dims)), 2)
    S ./= norm(ψ)
    S2 = diag(S) .^ 2
    return -1 * sum(S2 .* log.(S2))
end

"""
    entropy(ψ::StateVector, sites)

Calcualte the bipartition entanglement entropy between two partitions, one with `sites`
and the other with the remaining sites for StateVector `ψ`.
"""
function entropy(ψ::StateVector, sites)
    # Checks 
    for site in sites 
        if site < 0 || site >= length(ψ)
            throw(ArgumentError("site=$(site) is out of range."))
        end
    end

    # SVD 
    _, S, _ = tsvd(tensor(ψ), sites)
    S ./= norm(ψ)
    S2 = diag(S) .^ 2
    return -1 * sum(S2 .* log.(S2))
end

### Sampling a state vector 
export sample
"""
    sample(ψ::MPS)

Sample the StateVector `ψ` with the interpretation that it is a wavefunction (or Born ansatz).
"""
function sample(ψ::StateVector)
    Ps = cache(Float64, size(tensor(ψ)), level=2, sublevel=1) .= abs.(tensor(ψ)) .^ 2 ./ norm(ψ)
    Ps = reshape(Ps, length(tensor(ψ)))
    cum = 0.0
    r = rand()
    idx = 0
    for i in eachindex(Ps)
        idx = i
        cum = cum + Ps[i]
        if cum > r
            break
        end
    end
    return Tuple(CartesianIndices(tensor(ψ))[idx]), Ps[idx]
end