#=
    A circuit contains a full list of unitary operations and measurements that 
    can be applied to a tensor network state (or operator [later]).
=#

mutable struct Circuit{d}
    N::Int
    layers::Vector{CircuitLayer}
    con::CircuitConnectivity
end

export Circuit
Base.length(circuit::Circuit) = circuit.N 
dim(::Circuit{d}) where {d} = d
depth(circuit::Circuit) = length(circuit.layers)
Base.eltype(circuit::Circuit) = _promote_tensor_eltype(circuit.layers...)

### Adding to a circuit 
export add!, addlayer!
function addlayer!(circuit::Circuit)
    push!(circuit.layers, CircuitLayer(dim(circuit), length(circuit)))
end


function add!(circuit::Circuit, U::AbstractGate, sites)
    # We want to compress the circuit into as few layers as possible
    # We will check backwards to see how to add 
    idx = 0
    for i in reverse(eachindex(circuit.layers))
        if checkassigned(circuit.layers[i], sites, circuit.con)
            break
        end
        idx += 1
    end

    ### Add to the layer; if zero, create a new layer 
    if idx == 0
        addlayer!(circuit)
        add!(circuit.layers[end], U, sites, circuit.con)
    else
        add!(circuit.layers[end+1-idx], U, sites, circuit.con)
    end
end

### Initalising a circuit
export CircuitMPS, randombwcircuit

# Empty circuit 
function Circuit(d::Int, N::Int, connector::CircuitConnectivity=CircuitAll())
    return Circuit{d}(N, [], connector)
end


# Brickwall circuit 
"""
    randombwcircuit(d::Int, N::Int, depth::Int, width::Int=2; kwargs...)

Create a brickwall circuit composed of random unitary gates for a lattice
with physical dimension `d` and length `N`. The circuit has `depth` layers, and the gates
span `width` qubits.

# Optional Keyword Arguments
    - `ϵ::Number=0.01`: Random gates are generated as gates close to identity; the ϵ 
       is a paramater that controls how close to identity they are.
       - `connector::CircuitConnectivity=CircuitMPS()`: Add a circuit connector of choice.
       Default is CircuitMPS.
"""
function randombwcircuit(d::Int, N::Int, depth::Int, width::Int=2;
    ϵ::Number=0.01,
    connector::CircuitConnectivity=CircuitMPS()
)
    circuit = Circuit(d, N, connector)
    for m = Base.OneTo(depth)
        start = 1 + (m - 1) % 3
        for j = Base.range(start, N+1-width, step=width)
            add!(circuit, _unitary_close_to_id(d, width, ϵ), collect(Base.range(j, j+width-1)))
        end
    end

    return circuit
end

# Staircase circuit 
export randomstaircasecircuit

"""
    randomstaircasecircuit(d::Int, N::Int, depth::Int, width::Int=2; kwargs...)

Create a staircase circuit composed of random unitary gates for a lattice
with physical dimension `d` and length `N`. The circuit has `depth` layers, and the gates
span `width` qubits.

# Optional Keyword Arguments
    - `ϵ::Number=0.01`: Random gates are generated as gates close to identity; the ϵ 
       is a paramater that controls how close to identity they are.
       - `connector::CircuitConnectivity=CircuitMPS()`: Add a circuit connector of choice.
       Default is CircuitMPS.
"""
function randomstaircasecircuit(d::Int, N::Int, depth::Int, width::Int=2;
    ϵ::Number=0.01,
    connector::CircuitConnectivity=CircuitMPS()
)
    circuit = Circuit(d, N, connector)
    for m = Base.OneTo(depth)
        for j = Base.OneTo(N+1-width)
            add!(circuit, _unitary_close_to_id(d, width, ϵ), collect(Base.range(j, j+width-1)))
        end
    end

    return circuit
end