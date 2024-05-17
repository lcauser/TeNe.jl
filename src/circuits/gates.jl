abstract type AbstractGate end

### A customizable gate 
mutable struct CircuitGate{d, n, T} <: AbstractGate where {d, n, T<:AbstractArray}
    gate::T
end
export creategate, makeunitary!, makeunitary
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

# Polar decomposition 
"""
    makeunitary!(gate::CircuitGate)

Make a CircuitGate unitary via the polar decomposition, which finds the unitary matrix 
which most closely resembles the gate.

# Examples

```julia-repl 
julia> U = creategate(randn(ComplexF64, 2, 2, 2, 2));
julia> makeunitary!(U)
```
"""
function makeunitary!(gate::CircuitGate)
    U, _, V = tsvd(tensor(gate), Base.range(2, 2*length(gate), step=2))
    ten = contract(U, V, ndims(U), 1)
    permutedims!(tensor(gate), ten, map(j->isodd(j) ? cld(j, 2) : length(gate)+fld(j, 2), Base.OneTo(2*length(gate))))
end

"""
    makeunitary(gate::CircuitGate)

Make a CircuitGate unitary via the polar decomposition, which finds the unitary matrix 
which most closely resembles the gate. Returns the unitary gate as a new allocation.

# Examples

```julia-repl 
julia> g = creategate(randn(ComplexF64, 2, 2, 2, 2));
julia> U = makeunitary(g)
```
"""
function makeunitary(gate::CircuitGate)
    newgate = creategate(deepcopy(tensor(gate)))
    makeunitary!(newgate)
    return newgate
end


### Making random unitaries 

function _unitary_close_to_id(d::Int, N::Int, ϵ::Number=1e-1)
    # identity
    id = LinearAlgebra.diagm(ones(ComplexF64, d))
    H = ones(ComplexF64, )
    for i = 1:N
        H = tensorproduct(H, id; tocache=!(i==N))
    end
    
    H .+= ϵ*randn(ComplexF64, [d for _ = 1:2*N]...)
    U = creategate(H)
    makeunitary!(U)
    return U
end