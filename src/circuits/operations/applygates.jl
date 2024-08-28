#=
    Applying layers, or a circuit, of gates to a tensor network
=#


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

### Multiplying by a circuit 
function applygates!(circuit::Circuit, ψ::Union{TensorNetworkVector, TensorNetworkOperator}; kwargs...)
    for m in eachindex(circuit.layers)
        applygates!(circuit.layers[m], ψ; kwargs...)
    end
end

### Applying gates and making a copy
export applygates
function applygates(
    circuit::Union{Circuit, CircuitLayer},
    ψ::Union{TensorNetworkVector, TensorNetworkOperator};
    kwargs...
    )
    # Create a copy
    ψ′ = deepcopy(ψ)
    applygates!(circuit, ψ′; kwargs...)
    return ψ′
end

function applygates(
    ψ::Union{TensorNetworkVector, TensorNetworkOperator},
    circuit::Union{Circuit, CircuitLayer};
    kwargs...
    )
    # Create a copy
    ψ′ = deepcopy(ψ)
    applygates!(ψ′, circuit; kwargs...)
    return ψ′
end