#=
    Each qubit in a circuit layer can only be acted on with a single unitary
    (or measurement, to come later). A circuit, and circuit layer must be equipped
    with a CircuitConnectivity trait that contains how qubits are connected.
=#

mutable struct CircuitLayer{d}
    gates::Vector{<:AbstractGate}
    qubits::Vector{Tuple{Vararg{Int}}}
    assigned::Vector{Bool} 
end

export CircuitLayer
"""
    CircuitLayer(dim::Int, N::Int)

Create a ciruit layer for physical dimension `dim` and size `N`.
"""
function CircuitLayer(dim::Int, N::Int)
    return CircuitLayer{dim}(AbstractGate[], Tuple{Vararg{Int}}[], [false for _ = 1:N])
end