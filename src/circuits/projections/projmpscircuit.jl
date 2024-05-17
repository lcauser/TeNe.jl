#=
    A common strategy for updating quantum circuits requires us to partially evulate circuit expressions, such as
    <ϕ|U|ψ>. A ProjMPSCircuit calculates this projection onto a single gate within the expression (i.e. evaluates the
    inner product by contracting all the tensors in the network except for the gate in question), which can then be used
    in larger calcualtions, e.g., updating the gates to minimize some expectation.  

    This projection is calculated by applying circuit layers to the MPS approximately (depth-major).
=#

mutable struct ProjMPSCircuit{Q<:Number}
    # The terms in the inner product
    ϕ::MPS
    circuit::Circuit
    ψ::MPS
    λ::Q

    # Approximate layer applications 
    tops::Vector{MPS}
    bottoms::Vector{MPS}
    depth_center::Int
    cutoff::Float64
    maxdim::Int

    # Projection onto the gates using the MPS...
    lefts::Vector{<:AbstractArray}
    rights::Vector{<:AbstractArray}
    width_center::Int
end

function topblock(projU::ProjMPSCircuit, idx::Int)
    if idx < 0 || idx > depth(projU.circuit)
        throw(DomainError("The block must be between 0 and $(depth(projU.circuit))."))
    end
    if idx == 0 
        return projU.ϕ
    else
        return projU.tops[idx]
    end
end

function bottomblock(projU::ProjMPSCircuit, idx::Int)
    if idx < 1 || idx > depth(projU.circuit)+1
        throw(DomainError("The block must be between 1 and $(depth(projU.circuit)+1)."))
    end
    if idx == 0 
        return projU.ψ
    else
        return projU.bottoms[idx]
    end
end


function leftblock(projU::ProjMPSCircuit, idx::Int)
    if idx < 0 || idx > length(projU.circuit.layers[projU.depth_center])
        throw(DomainError("The block must be between 0 and $(length(projU.circuit.layers[projU.depth_center]))."))
    end
    return projU.lefts[idx+1]
end


function rightblock(projU::ProjMPSCircuit, idx::Int)
    if idx < 1 || idx > length(projU.circuit.layers[projU.depth_center]+1)
        throw(DomainError("The block must be between 0 and $(length(projU.circuit.layers[projU.depth_center]))."))
    end
    return projU.lefts[idx]
end


### Builing blocks
function buildtop!(projU::ProjMPSCircuit, idx::Int)
    projU.tops[idx+1] = deepcopy(topblock(projU, idx))
    applygates!(projU.tops[idx+1], projU.circuit.layers[end+1-idx]; cutoff=projU.cutoff, maxdim=projU.maxdim)
end

function buildbottom!(projU::ProjMPSCircuit, idx::Int)
    projU.bottoms[idx] = deepcopy(bottomblock(projU, idx+1))
    applygates!(projU.circuit.layers[end+1-idx], projU.bottoms[idx]; cutoff=projU.cutoff, maxdim=projU.maxdim)
end