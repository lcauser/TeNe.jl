#=
    Creates a circuit for time evolution using Trotterization
=#

function trotterize(H::OpList, t::Number; order::Int=2, type::Symbol=:compressed)

end


### Creating the Trotter gate 

function _create_trotter_gate(H::OpList, t::Number, site::Int)
    # Create a tensor to store the result
    rng = siterange(H) 
    ten = zeros(eltype(H.lt), [dim(H) for _ = 1:2*rng]...)

    # Find all the gates that start at the site with the correct range 
    rng = siterange(H)
    idxs = siteindexs(H, site)
    


end