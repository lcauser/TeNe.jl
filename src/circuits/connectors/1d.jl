#=
    The circuit connector for a 1d lattice with neighbouring connectivity.
=#

### One-dimensional circuit 
struct Circuit1D <: CircuitConnectivity end

"""
    add!(layer::CircuitLayer, U::AbstractGate, sites, ::Circuit1D)

Add a gate `U` at lattice sites `sites` to a CircuitLayer. The connector only allows
for, at most, nearest neighbour gates.
"""
function add!(layer::CircuitLayer, U::AbstractGate, sites, con::Circuit1D)
    # Validation 
    _validate_gate_layer(layer, gate, sites)
    if length(sites) > 2
        throw(
            ArgumentError(
                "The Circuit1D connector only allows for nearest-neighbour gates.",
            ),
        )
    elseif length(sites) == 2
        if abs(sites[begin]-sites[end]) != 1
            throw(
                ArgumentError(
                    "The Circuit1D connector only allows for nearest-neighbour gates.",
                ),
            )
        end
    end

    # Check if assigned
    if checkassigned(layer, sites, con)
        throw(
            ArgumentError(
                "One or more of the sites $(sites) in the circuit layer are already used.",
            ),
        )
    end
    push!(layer.gates, U)
    push!(layer.sites, Tuple(sites))
    for site in sites
        layer.assigned[site] = true
    end
end
