#=
    Creates a circuit for time evolution using Trotterization
=#

"""
    trotterize(H::OpList, t::Number; kwargs...)

Trotterize the Hamiltonian `H` for time step `t`. That is, approximate exp(tH) by
a discrete circuit.

# Optional Keyword Arguments 

    - `order::Int=2`: The order of the approximation; either `1` or `2`.
    - `type::Symbol=:compressed`: Use `uncompressed` for a raw circuit. Use `compressed`
       to treat terms with a smaller site range as a term with a larger range.
       In theory, `compressed` should give some prefactor improvements on the error.

# Examples 

```julia-repl
julia> H = OpList(Qubits(), 20);
julia> for i = 1:20 add!(H, "x", i, 0.34); end
julia> for i = 1:19 add!(H, ["z", "z"], [i, i+1], 0.56); end
julia> circuit = trotterize(H, -0.01im);
```
"""
function trotterize(H::OpList, t::Number; order::Int=2, type::Symbol=:compressed)
    # Validation on the hyperparameters 
    if !(order > 0 && order <= 2)
        throw(ArgumentError("Only first and second order are currently supported."))
    end
    if type != :compressed && type != :uncompressed
        throw(ArgumentError("The type of circuit must be either compressed or uncompressed."))
    end

    # Choose the correct trotterization 
    if type == :compressed && order == 1
        return _trotter_compressed_first_order(H, t)
    elseif type == :compressed && order == 2
        return _trotter_compressed_second_order(H, t)
    elseif  type == :uncompressed && order == 1
        return _trotter_uncompressed_first_order(H, t)
    else
        return _trotter_uncompressed_second_order(H, t)
    end
end
export trotterize


### Creating the Trotter gate in a compressed way
# First-order decomposition with compressed gates
function _trotter_compressed_first_order(H::OpList, t::Number)
    # Create the circuit 
    circuit = Circuit(dim(H), length(H))

    # Add the layers 
    for i = 1:siterange(H)
        push!(circuit.layers, _create_trotter_layer_compressed(H, t, i))
    end
    return circuit
end

# Second-order decomposition with compressed gates 
function _trotter_compressed_second_order(H::OpList, t::Number)
    # Create the circuit 
    circuit = Circuit(dim(H), length(H))

    # Add the layers 
    for i = 1:siterange(H)-1
        push!(circuit.layers, _create_trotter_layer_compressed(H, t/2, i))
    end
    push!(circuit.layers, _create_trotter_layer_compressed(H, t, siterange(H)))
    for i = Base.range(siterange(H)-1, 1, step=-1)
        push!(circuit.layers, circuit.layers[i])
    end
    return circuit
end

# Create a compressed gate that starts at ``site``
function _create_trotter_gate_compressed(H::OpList, t::Number, site::Int)
    # Create a tensor to store the result
    rng = siterange(H) 
    ten = zeros(eltype(H.lt), [dim(H) for _ = 1:2*rng]...)

    # Loop through all sites and find the operators in the range 
    for i in Base.OneTo(rng)
        # Find operators which start at this index
        idxs = siteindexs(H, site + i - 1)
        for idx in idxs 
            # Make sure the operator has small enough range 
            rng_local = H.sites[idx][end] - H.sites[idx][begin] + 1
            padding = rng - (i - 1 + rng_local)
            if padding >= 0
                # Determine the coefficient to ensure the term is not over counted 
                first_site = max(site + i - 1 + rng_local - rng, 1)
                last_site = min(length(H) - rng + 1, site + i - 1)
                coeff = last_site - first_site + 1
                
                # Find the operator 
                op = totensor(H, idx, i-1, padding; tocache=true)
                op .*= (1 / coeff)
                ten .+= op
            end
        end
    end

    # Exponentiate the tensor 
    ten = TeNe.exp(ten, Tuple(Base.range(2, 2*rng, step=2)); prefactor=t)
    return ten    
end

# Create a layer of compressed gates; parity specifies the sequence of sites.
function _create_trotter_layer_compressed(H::OpList, t::Number, parity::Int)
    rng = siterange(H)
    layer = CircuitLayer(dim(H), length(H))
    for j = Base.range(parity, length(H)-rng+1, step=rng)
        add!(layer, creategate(_create_trotter_gate_compressed(H, t, j)),
            Tuple(Base.range(j, j+rng-1)))
    end
    return layer
end

### Creating the trotter gate without compression 
# First-order decomposition with uncompressed gates
function _trotter_uncompressed_first_order(H::OpList, t::Number)
    # Create the circuit 
    circuit = Circuit(dim(H), length(H))

    # Add the layers 
    for k = 1:siterange(H)
        for i = 1:k
            push!(circuit.layers, _create_trotter_layer_uncompressed(H, t, i, k))
        end
    end
    return circuit
end

# Second-order decomposition with uncompressed gates 
function _trotter_uncompressed_second_order(H::OpList, t::Number)
    # Create the circuit 
    circuit = Circuit(dim(H), length(H))

    # Add the layers 
    for k = 1:siterange(H)
        for i = 1:min(k, siterange(H)-1)
            push!(circuit.layers, _create_trotter_layer_uncompressed(H, t/2, i, k))
        end
    end
    push!(circuit.layers, _create_trotter_layer_uncompressed(H, t, siterange(H), siterange(H)))
    for k = Base.range(length(circuit.layers)-1, 1, step=-1)
        push!(circuit.layers, circuit.layers[k])
    end
    return circuit
end

# Creating an uncompressed gate 
function _create_trotter_gate_uncompressed(H::OpList, t::Number, site::Int, rng::Int)
    # Create a tensor to store the result 
    ten = zeros(eltype(H.lt), [dim(H) for _ = 1:2*rng]...)

    # Find operators which start at this index
    idxs = siteindexs(H, site)
    for idx in idxs 
        # Make sure the operator has small enough range 
        rng_local = H.sites[idx][end] - H.sites[idx][begin] + 1
        if rng_local == rng
            ten .+= totensor(H, idx, 0, 0; tocache=true)
        end
    end

    # Exponentiate the tensor 
    ten = TeNe.exp(ten, Tuple(Base.range(2, 2*rng, step=2)); prefactor=t)
    return ten  
end

# Create a layer of uncompressed gates; parity specifies the sequence of sites.
function _create_trotter_layer_uncompressed(H::OpList, t::Number, parity::Int, rng::Int)
    layer = CircuitLayer(dim(H), length(H))
    for j = Base.range(parity, length(H)-rng+1, step=rng)
        add!(layer, creategate(_create_trotter_gate_uncompressed(H, t, j, rng)),
            Tuple(Base.range(j, j+rng-1)))
    end
    return layer
end