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
    T::Type

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

    # Create the top and bottom MPSs
    tops = MPS[MPS(dim(ψ), length(ψ)) for _ = 1:depth(circuit)]
    bottoms = MPS[MPS(dim(ψ), length(ψ)) for _ = 1:depth(circuit)]

    projU = ProjMPSCircuit(
        conj(ϕ), circuit, ψ, λ, _promote_tensor_eltype(ϕ, circuit, ψ, λ),
        tops, bottoms, 0, cutoff, maxdim,
        AbstractArray[], AbstractArray[], 0
    )
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
    if idx < 0 || idx > length(getlayer(projU, projU.depth_center).gates)
        throw(DomainError("The block must be between 0 and $(length(getlayer(projU, projU.depth_center).gates))."))
    end
    if idx == 0 
        return ones(Float64, 1, 1)
    end
    return projU.lefts[idx]
end


function rightblock(projU::ProjMPSCircuit, idx::Int)
    if idx < 1 || idx > length(getlayer(projU, projU.depth_center).gates)+1
        throw(DomainError("The block must be between 1 and $(length(getlayer(projU, projU.depth_center).gates) + 1)."))
    end
    if idx == length(getlayer(projU, projU.depth_center).gates)+1
        return ones(Float64, 1, 1)
    end
    return projU.lefts[idx]
end


### Fetching the layer 
function getlayer(projU::ProjMPSCircuit, idx::Int)
    return projU.circuit.layers[end + 1 - idx]
end


### Builing top and bottom MPSs
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


### Building left and right blocks 
function _create_blocks!(projU::ProjMPSCircuit)
    lefts = []
    rights = []
    layer = getlayer(projU, projU.depth_center)
    for i = Base.OneTo(length(layer.gates))
        sitebefore = layer.sites[i][begin] - 1
        siteafter = layer.sites[i][end] + 1
        push!(lefts, ones(projU.T, bonddim(topblock(projU, projU.depth_center-1), sitebefore), 
            bonddim(bottomblock(projU, projU.depth_center+1), sitebefore)))
        push!(rights, ones(projU.T, bonddim(topblock(projU, projU.depth_center-1), siteafter), 
            bonddim(bottomblock(projU, projU.depth_center+1), siteafter)))
    end
    projU.lefts = lefts 
    projU.rights = rights
end


function buildleft!(projU::ProjMPSCircuit, idx::Int)
    projU.lefts[idx] .= _buildleft(projU, idx, leftblock(projU, idx-1))
end

function _buildleft(projU::ProjMPSCircuit, idx::Int, left)
    # Fetch information
    top_mps = topblock(projU, projU.depth_center-1)
    bottom_mps = topblock(projU, projU.depth_center+1)
    layer = getlayer(projU, projU.depth_center)
    finalidx = layer.sites[idx][end]

    # Contract in the previous gate 
    if idx != 1
        # Fetch the previous gate 
        sites_prev = layer.sites[idx-1]
        gate_prev = layer.gates[idx-1]

        # Contract with MPS tensors 
        for site = Base.range(sites_prev[begin], sites_prev[end])
            left = contract(left, top_mps[site], ndims(left)-1, 1, false, true)
            if site in sites_prev
                # Site is affected by gate
                left = contract(left,  bottom_mps[site], ndims(left)-2, 1)
                left = permutedim(left, ndims(left)-2, ndims(left)-1)
            else
                # Site is not affected by gate
                left = contract(left,  bottom_mps[site], (ndims(left)-2, ndims(left)-1), (1, 2))
            end
        end

        # Contract with gate 
        left = contract(left, gate_prev, Base.OneTo(ndims(gate_prev)), Base.OneTo(ndims(gate_prev)))

        # Set the starting
        startidx = sites_prev[end] + 1
    else
        startidx = 1
    end

    # Contract the next sites leading up to the final site 
    for site = Base.range(startidx, finalidx)
        left = contract(left, top_mps[site], ndims(left)-1, 1, false, true)
        left = contract(left, bottom_mps[site], (1, 2), (1, 2))
    end

end

function _buildright(projU::ProjMPSCircuit, idx::Int, right)
    # Fetch information
    top_mps = topblock(projU, projU.depth_center-1)
    bottom_mps = topblock(projU, projU.depth_center+1)
    layer = getlayer(projU, projU.depth_center)
    finalidx = laye.sites[idx][begin]

    # Contract previous gate
    if idx != length(layer.gates)
        # Fetch the previous gate 
        sites_prev = layer.sites[idx+1]
        gate_prev = layer.gates[idx+1]

        # Contract with MPS tensors 
        for site = Base.range(sites_prev[end], sites_prev[begin], step=-1)
            left = contract(left, top_mps[site], ndims(left)-1, 1, false, true)
            if site in sites_prev
                # Site is affected by gate
                left = contract(left,  bottom_mps[site], ndims(left)-2, 1)
                left = permutedim(left, ndims(left)-2, ndims(left)-1)
            else
                # Site is not affected by gate
                left = contract(left,  bottom_mps[site], (ndims(left)-2, ndims(left)-1), (1, 2))
            end
        end

        # Contract with gate 
        left = contract(left, gate_prev, Base.OneTo(ndims(gate_prev)), Base.OneTo(ndims(gate_prev)))

        # Set the starting
        startidx = sites_prev[end] + 1
    else
        startidx = 1
    end
end


### Moving the center 
function movecenter!(projU::ProjMPSCircuit, idx_depth::Int, idx_width::Int=0)
    _movecenter_depth!(projU, idx_depth)
end

function _movecenter_depth!(projU::ProjMPSCircuit, idx::Int)
    # Break if not chanign the depth 
    if idx == projU.depth_center
        return nothing
    end

    # Build the depth
    if projU.depth_center == 0
        for i in Base.OneTo(idx-1)
            buildtop!(projU, i)
        end
        for i in Base.range(depth(projU.circuit), idx+1, step=-1)
            buildbottom!(projU, i)
        end
    elseif idx > projU.depth_center
        for i in Base.range(projU.depth_center, idx-1)
            buildtop!(projU, i)
        end
    elseif idx < projU.depth_center
        for i = Base.range(projU.depth_center, idx+1, step=-1)
            buildbottom!(projU, i)
        end
    end
    projU.depth_center = idx

    # Reset the width 
    #_create_blocks!(projU)
    #projU.width_center = 0
end
