#=
    The circuit connector for matrix product states.
=#

### One-dimensional circuit 
struct CircuitMPS <: CircuitConnectivity end

"""
    add!(layer::CircuitLayer, U::AbstractGate, sites, ::CircuitMPS)

Add a gate `U` at lattice sites `sites` to a CircuitLayer. The connector 
is used for MPS, which maps the lattice sites inbetween the given `sites`.
"""
function add!(layer::CircuitLayer, U::AbstractGate, sites, con::CircuitMPS)
    # Validation 
    _validate_gate_layer(layer, gate, sites)

    # Check if assigned
    sitesfull = Tuple(Base.range(min(sites...), max(sites...)))
    if checkassigned(layer, sitesfull, con)
        throw(ArgumentError("One or more of the sites $(sites) in the circuit layer are already used."))
    end

    # Assign the sites
    for site in sitesfull 
        layer.assigned[site] = true
    end

    # Add the tensor 
    U, sites = _order_gate(U, sites)
    push!(layer.gates, U)
    push!(layer.sites, Tuple(sites))
end

function _order_gate(U::AbstractGate, sites)
    sites = [sites...]
    perms = sortperm(sites)
    if all(j->perms[j]==j, Base.eachindex(perms))
        return U, sites
    end
    U2 = deepcopy(U)
    tperms = map(j->2*perms[cld(j, 2)] - isodd(j), Tuple(Base.OneTo(2*Base.length(sites))))
    tensor(U2) .= _permutedims(tensor(U), tperms)
    return U2, Tuple(sites[perms])
end