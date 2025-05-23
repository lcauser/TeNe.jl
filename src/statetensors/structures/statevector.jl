#=
    State vectors can be used to describes states of many-body systems with discrete
    degress of freedom.
=#

const StateVector = Union{GStateTensor{1},ConjGStateTensor{1}}
export StateVector

"""
    isstatevector(ψ)

Check to see if an object is a state vector.
"""
function isstatevector(ψ)
    return typeof(ψ) <: StateVector
end

export dim
dim(ψ::StateVector, site::Int) = dim(ψ, 1, site)
dims(ψ::StateVector) = dims(ψ, 1)


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

"""
    productsv(lt::LatticeTypes, states::AbstractVector{String})
    productstatevector(lt::LatticeTypes, states::AbstractVector{String})

Create a product state from a string of state names.

# Example 

```julia-repl
julia> lt = Qubits();
julia> ψ = productsv(lt, ["up" for _ = 1:6]);
```
"""
function productsv(lt::LatticeTypes, states::AbstractVector{String})
    tensor = ones(eltype(lt))
    for i in eachindex(states)
        tensor = tensorproduct(tensor, state(lt, states[i]); tocache = i!=lastindex(states))
    end
    return GStateTensor(1, tensor)
end
productstatevector(lt::LatticeTypes, states::AbstractVector{String}) = productsv(lt, states)


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
    Ps =
        cache(Float64, size(tensor(ψ)), level = 2, sublevel = 1) .=
            abs.(tensor(ψ)) .^ 2 ./ norm(ψ)
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
