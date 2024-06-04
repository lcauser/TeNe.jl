#=
    Creates a circuit for time evolution using Trotterization
=#

function trotterize(H::OpList, t::Number; order::Int=2, type::Symbol=:compressed)

end


### Creating the Trotter gate in a compressed way
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

### Creating the trotter gate without compression 
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