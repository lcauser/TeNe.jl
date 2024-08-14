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
    tops::Vector{MPS} # ϕs
    bottoms::Vector{MPS} # ψs
    depth_center::Int
    cutoff::Float64
    maxdim::Int

    # Projection onto the gates using the MPS...
    lefts::Vector{<:AbstractArray}
    rights::Vector{<:AbstractArray}
    width_center::Int
end
export ProjMPSCircuit 

function ProjMPSCircuit(ϕ::MPS, circuit::Circuit, ψ::MPS; cutoff::Float64=1e-12,
    maxdim::Int=0, depth_center::Int=1, width_center::Int=1, λ::Number=1.0)

    tops = MPS[MPS(dim(ψ), length(ψ)) for _ = 1:depth(circuit)]
    bottoms = MPS[MPS(dim(ψ), length(ψ)) for _ = 1:depth(circuit)]
    projU = ProjMPSCircuit(ϕ, circuit, ψ, λ, tops, bottoms, 0, cutoff, maxdim, AbstractArray[], AbstractArray[], 0)
    movecenter!(projU, depth_center)
        
    return projU
end


### Fetch the blocks for the projection
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
    if idx == depth(projU.circuit)+1 
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
    projU.tops[idx] = _buildtop(deepcopy(topblock(projU, idx-1)), projU, idx)
end

function _buildtop(ϕ::MPS, projU::ProjMPSCircuit, idx::Int)
    applygates!(ϕ, projU.circuit.layers[end+1-idx]; cutoff=projU.cutoff, maxdim=projU.maxdim)
    return ϕ
end

function buildbottom!(projU::ProjMPSCircuit, idx::Int)
    projU.bottoms[idx] = _buildbottom(deepcopy(bottomblock(projU, idx+1)), projU, idx)
end

function _buildbottom(ψ::MPS, projU::ProjMPSCircuit, idx::Int)
    applygates!(projU.circuit.layers[end+1-idx], ψ; cutoff=projU.cutoff, maxdim=projU.maxdim)
    return ψ
end



### Moving the center 
function movecenter!(projU::ProjMPSCircuit, idx_depth::Int, idx_width::Int=0)
    _movecenter_depth!(projU, idx_depth)
end

function _movecenter_depth!(projU::ProjMPSCircuit, idx::Int)
    if projU.depth_center == 0
        for i in Base.OneTo(idx-1)
            buildtop!(projU, i)
        end
        for i in Base.range(depth(projU.circuit), idx+1, step=-1)
            buildbottom!(projU, i)
        end
    else
        if idx > projU.depth_center
            for i in Base.range(projU.depth_center, idx-1)
                buildtop!(projU, i)
            end
        elseif idx < projU.depth_center
            for i = Base.range(projU.depth_center, idx+1, step=-1)
                buildbottom!(projU, i)
            end
        end
    end
    projU.depth_center = idx
end
