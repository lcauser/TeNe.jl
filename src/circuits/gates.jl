abstract type AbstractGate end

mutable struct CircuitGate{d, n, T} <: AbstractGate where {d, n, T<:AbstractArray}
    gate::T
end

export creategate
"""
    creategate(gate::AbstractArray)

Create a gate using a tensor,
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