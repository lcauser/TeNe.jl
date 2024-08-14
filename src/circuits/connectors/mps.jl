#=
    The circuit connector for matrix product states.
=#

### One-dimensional circuit 
struct CircuitMPS <: CircuitConnectivity end
export CircuitMPS

"""
    add!(layer::CircuitLayer, U::AbstractGate, sites, ::CircuitMPS)

Add a gate `U` at lattice sites `sites` to a CircuitLayer. The connector 
is used for MPS, which maps the lattice sites inbetween the given `sites`.
"""
function add!(layer::CircuitLayer, U::AbstractGate, sites, con::CircuitMPS)
    # Validation 
    _validate_gate_layer(layer, U, sites)

    # Check if assigned
    sitesfull = Tuple(Base.range(min(sites...), max(sites...)))
    if checkassigned(layer, sitesfull, con)
        throw(ArgumentError("One or more of the sites $(sites) in the circuit layer are already used."))
    end

    # Assign the sites
    for site in sitesfull 
        layer.assigned[site] = true
    end

    # Order the sites in the gate
    U, sites = _order_gate(U, sites)

    # Find the position to insert 
    idx = 1
    for s in layer.sites 
        if s[begin] > sitesfull[end]
            break
        end
        idx += 1
    end
    
    # Insert the gate
    insert!(layer.gates, idx, U)
    insert!(layer.sites, idx, Tuple(sites))
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