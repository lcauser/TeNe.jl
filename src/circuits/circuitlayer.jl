#=
    Each qubit in a circuit layer can only be acted on with a single unitary
    (or measurement, to come later). A circuit, and circuit layer must be equipped
    with a CircuitConnectivity trait that contains how qubits are connected.
=#

mutable struct CircuitLayer{d}
    gates::Vector{<:AbstractGate}
    sites::Vector{Tuple{Vararg{Int}}}
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

dim(::CircuitLayer{d}) where {d} = d
Base.length(layer::CircuitLayer) = Base.length(layer.assigned)

### Connectors
abstract type CircuitConnectivity end
struct CircuitAll <: CircuitConnectivity end
export CircuitAll
export checkassigned, add!


"""
    checkassigned(layer::CircuitLayer, sites, [::CircuitConnectivity])

Check if the lattice sites `sites` are already used in a CircuitLayer.
"""
function checkassigned(layer::CircuitLayer, sites, ::CircuitConnectivity=CircuitAll())
    for site in sites
        if layer.assigned[site]
            return true
        end
    end
    return false
end

"""
    add!(layer::CircuitLayer, U::AbstractGate, sites, [con::CircuitConnectivity])

Add a gate `U` at lattice sites `sites` to a CircuitLayer. The connectivity of the lattice 
can be specified.
"""
function add!(layer::CircuitLayer, U::AbstractGate, sites, con::CircuitConnectivity)
    _validate_gate_layer(layer, U, sites)
    if checkassigned(layer, sites, con)
        throw(ArgumentError("One or more of the sites $(sites) in the circuit layer are already used."))
    end
    push!(layer.gates, U)
    push!(layer.sites, Tuple(sites))
    for site in sites 
        layer.assigned[site] = true
    end
end
add!(layer::CircuitLayer, U::AbstractGate, sites) = add!(layer, U, sites, CircuitAll())

function _validate_gate_layer(layer::CircuitLayer, gate::AbstractGate, sites)
    if length(sites) != length(gate)
        throw(ArgumentError("The length of the gate and sites do not match."))
    end
    if dim(layer) != dim(gate)
        throw(ArgumentError("The circuit layer and gate have contradicting physical dimensions."))
    end
    for site in sites 
        if site <= 0 || site > length(layer)
            throw(DomainError("Site $(site) is out of range."))
        end
    end
end

### Applying a circuit layer to an mps 
export applygates!
function applygates!(layer::CircuitLayer, ψ::Union{MPS, MPO}; kwargs...)
    # Decide on the sweeping direction 
    rev = center(ψ) > length(ψ) / 2

    # Loop through each gate
    lst = eachindex(layer.sites)
    lst = rev ? reverse(lst) : lst
    for i in lst
        applygate!(layer.gates[i], ψ, layer.sites[i][rev ? end : begin], rev; kwargs...)
    end
end

function applygates!(ψ::Union{MPS, MPO}, layer::CircuitLayer; kwargs...)
    # Decide on the sweeping direction 
    rev = center(ψ) > length(ψ) / 2

    # Loop through each gate
    lst = eachindex(layer.sites)
    lst = rev ? reverse(lst) : lst
    for i in lst
        applygate!(ψ, layer.gates[i], layer.sites[i][rev ? end : begin], rev; kwargs...)
    end
end

### Applying a circuit layer to a state vector
function applygates!(layer::CircuitLayer, ψ::Union{StateVector, StateOperator}; kwargs...)
    for i in eachindex(layer.sites)
        applygate!(layer.gates[i], ψ, layer.sites[i])
    end
end

function applygates!(ψ::Union{StateVector, StateOperator}, layer::CircuitLayer; kwargs...)
    for i in eachindex(layer.sites)
        applygate!(ψ, layer.gates[i], layer.sites[i])
    end
end