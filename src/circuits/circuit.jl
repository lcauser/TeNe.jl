#=
    A circuit contains a full list of unitary operations and measurements that 
    can be applied to a tensor network state (or operator [later]).
=#

mutable struct Circuit{d}
    N::Int
    layers::Vector{CircuitLayer}
end
Circuit(d::Int, N::Int) = Circuit{d}(N, CircuitLayer[])

export Circuit 
Base.length(circuit::Circuit) = circuit.N 
dim(::Circuit{d}) where {d} = d
depth(circuit::Circuit) = length(circuit.layers)


### Initalising a circuit 
export randombwcircuit
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