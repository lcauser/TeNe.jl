abstract type AbstractGate end

### Applying gates to StateVectors 
# Safe application of a gate 
export applygate!
"""
    applygate!(U::AbstractGate, ψ::StateVector, sites)

Apply a circuit gate `U` to the StateVector `ψ` at lattice sites `sites`.
"""
function applygate!(U::AbstractGate, ψ::StateVector, sites)
    # Check sites 
    ψlen = length(ψ)
    if any(map(j->(j>ψlen || j <= 0), sites))
        throw(ArgumentError("The list of sites $(sites) does not fall between 1 and $(length(ψ))."))
    end

    # Check dimensions of sites 
    if any(map(j->dim(ψ, j)!=dim(U), sites))
        throw(ArgumentError("The StateVector has the wrong physical dimensions."))
    end

    # Apply the gate
    _applygate!(U, ψ, sites)
end

# Unsafe gate application
function _applygate!(U::AbstractGate, ψ::StateVector, sites)
    rixs = Tuple(setdiff(Base.OneTo(length(ψ)), sites)) # Remove allocations?
    ψ′ = contract(tensor(U), tensor(ψ), Tuple(Base.range(2, ndims(tensor(U)), step=2)),
        sites, false, isconj(ψ))
    permutedims!(tensor(ψ), ψ′, reverseperms((sites..., rixs...)))
end


### A customizable gate 
mutable struct CircuitGate{d, n, T} <: AbstractGate where {d, n, T<:AbstractArray}
    gate::T
end
export creategate
"""
    creategate(gate::AbstractArray)

Create a gate using a tensor.
"""
function creategate(gate::AbstractArray)
    # Validation 
    num_qubits = fld(ndims(gate), 2)
    if ndims(gate) % 2 != 0
        throw(ArgumentError("The gates must have an even number of dimensions."))
    end
    sz = size(gate)
    if !all(j->j==sz[1], sz)
        throw(ArgumentError("The gate must have the same dimensions sizes."))
    end
    # 
    return CircuitGate{sz[1], num_qubits, typeof(gate)}(gate)
end

# Properties 
export tensor, dim, qubits
tensor(gate::CircuitGate) = gate.gate
dim(::CircuitGate{d}) where {d} = d 
Base.length(::CircuitGate{d, n}) where {d, n} = n