#=
    A circuit contains a full list of unitary operations and measurements that 
    can be applied to a tensor network state (or operator [later]).
=#

mutable struct Circuit{d}
    N::Int
    layers::Vector{CircuitLayer}
    conn::CircuitConnectivity
end

export Circuit 
Base.length(circuit::Circuit) = circuit.N 
dim(::Circuit{d}) where {d} = d
depth(circuit::Circuit) = length(circuit.layers)

##


### Initalising a circuit 
# Empty circuit 
function Circuit(d::Int, N::Int, connector::CircuitConnectivity=CircuitAll())
    return Circuit{d}(N, [], connector)
end

function CircuitMPS(d::Int, N::Int)
    return Circuit{d}(N, [], CircuitMPS())
end

# Brickwall circuit 
export randombwcircuit

"""
    randombwcircuit(d::Int, N::Int, depth::Int; ϵ::Number=0.01)

Create a brickwall circuit composed of random unitary two-body gates for a lattice
with physical dimension `d` and length `N`. The circuit has `depth` layers.

# Optional Keyword Arguments
    - `ϵ::Number=0.01`: Random gates are generated as gates close to identity; the ϵ 
       is a paramater that controls how close to identity they are.
"""
function randombwcircuit(d::Int, N::Int, depth::Int; ϵ::Number=0.01)
    circuit = Circuit{d}(N, CircuitLayer[])
    for m = Base.OneTo(depth)
        layer = CircuitLayer(d, N)
        for j = Base.range(isodd(m) ? 1 : 2, N-1, step=2)
            add!(layer, _unitary_close_to_id(d, 2, ϵ), (j, j+1))
        end
        push!(circuit.layers, layer)
    end
    return circuit
end