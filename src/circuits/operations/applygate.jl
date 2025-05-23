#=
    Contains functionality for applying a gate to tensor network objects/
=#

### Applying gates to StateVectors 
# Safe application of a gate 
export applygate!
"""
    applygate!(U::AbstractGate, ψ::StateVector, sites)
    applygate!(ψ::StateVector, U::AbstractGate, sites)

Apply a circuit gate `U` to the StateVector `ψ` at lattice sites `sites`.
"""
function applygate!(U::AbstractGate, ψ::StateVector, sites)
    _gate_vec_validation(U, ψ, sites)
    _applygate!(U, ψ, sites)
end

function applygate!(ψ::StateVector, U::AbstractGate, sites)
    _gate_vec_validation(U, ψ, sites)
    _applygate!(ψ, U, sites)
end

# Unsafe gate application
function _applygate!(U::AbstractGate, ψ::StateVector, sites)
    # Contraction indices
    ris = Tuple(setdiff(Base.OneTo(length(ψ)), sites)) # Remove allocations?
    uis = Base.range(2, ndims(tensor(U)), step = 2)

    # Do the contraction
    ψ′ = contract(tensor(U), tensor(ψ), uis, sites, isconj(ψ), false)
    permutedims!(tensor(ψ), ψ′, reverseperms((sites..., ris...)))
end

function _applygate!(ψ::StateVector, U::AbstractGate, sites)
    # Contraction indices
    ris = Tuple(setdiff(Base.OneTo(length(ψ)), sites)) # Remove allocations?
    uis = Base.range(1, ndims(tensor(U)), step = 2)

    # Do the contraction
    ψ′ = contract(tensor(U), tensor(ψ), uis, sites, isconj(ψ), false)
    permutedims!(tensor(ψ), ψ′, reverseperms((sites..., ris...)))
end


### Applying gates to StateOperators 

"""
    applygate!(U::AbstractGate, O::StateOperator, sites)
    applygate!(O::StateOperator, U::AbstractGate, sites)

Apply a circuit gate `U` to the StateOperator `O` at lattice sites `sites`.
"""
function applygate!(U::AbstractGate, O::StateOperator, sites)
    _gate_op_validation(U, O, sites)
    _applygate!(U, O, sites)
end

function applygate!(O::StateOperator, U::AbstractGate, sites)
    _gate_op_validation(U, O, sites)
    _applygate!(O, U, sites)
end

# Unsafe gate application 
function _applygate!(U::AbstractGate, O::StateOperator, sites)
    # Contraction indices
    Ois = Tuple(map(j->istranspose(O) ? outerind(O, j) : innerind(O, j), sites))
    ris = Tuple(setdiff(Base.OneTo(ndims(tensor(O))), Ois))
    uis = Base.range(istranspose(O) ? 1 : 2, ndims(tensor(U)), step = 2)

    # Do the contraction 
    O′ = contract(tensor(U), tensor(O), uis, Ois, isconj(O), false)
    permutedims!(tensor(O), O′, reverseperms((Ois..., ris...)))
end

function _applygate!(O::StateOperator, U::AbstractGate, sites)
    # Contraction indices
    Ois = Tuple(map(j->istranspose(O) ? innerind(O, j) : outerind(O, j), sites))
    ris = Tuple(setdiff(Base.OneTo(ndims(tensor(O))), Ois))
    uis = Base.range(istranspose(O) ? 2 : 1, ndims(tensor(U)), step = 2)

    # Do the contraction 
    O′ = contract(tensor(O), tensor(U), Ois, uis, false, isconj(O))
    permutedims!(tensor(O), O′, reverseperms((ris..., Ois...)))
end


### Applying gates to MPS 
"""
    applygate!(U::AbstractGate, ψ::MPS, site::Int, [rev::Bool=false]; kwargs...)
    applygate!(ψ::MPS, U::AbstractGate, site::Int, [rev::Bool=false]; kwargs...)

Apply a circuit gate `U` to an MPS `ψ`. The bool `rev` specifies the sweeping direction of 
the MPS: `rev=false` for sweeping right, and `rev=true` for sweeping left.
The first (or last for `rev=true`) site that the gate is applied to is `site`.
"""
function applygate!(U::AbstractGate, ψ::MPS, site::Int, rev::Bool = false; kwargs...)
    _gate_vec_validation(
        U,
        ψ,
        Tuple(
            rev ? Base.range(site+1-length(U), site) : Base.range(site, site+length(U)-1),
        ),
    )
    _applygate!(U, ψ, site, rev; kwargs...)
end

function applygate!(ψ::MPS, U::AbstractGate, site::Int, rev::Bool = false; kwargs...)
    _gate_vec_validation(
        U,
        ψ,
        Tuple(
            rev ? Base.range(site+1-length(U), site) : Base.range(site, site+length(U)-1),
        ),
    )
    _applygate!(U, ψ, site, rev, true; kwargs...)
end

# Unsafe application 
function _applygate!(
    U::AbstractGate,
    ψ::MPS,
    site::Int,
    rev::Bool = false,
    trans::Bool = false;
    kwargs...,
)
    # Properties 
    num_sites = length(U)
    firstsite = rev ? site - num_sites + 1 : site
    lastsite = rev ? site : site + num_sites - 1

    # Move the center 
    if !(center(ψ) >= firstsite && center(ψ) <= lastsite)
        if center(ψ) < firstsite
            movecenter!(ψ, firstsite)
        else
            movecenter!(ψ, lastsite)
        end
    end

    # Contract MPS tensors
    ten = ψ[firstsite]
    for i in Base.range(firstsite+1, lastsite)
        ten = contract(ten, ψ[i], ndims(ten), 1)
    end

    # Apply the gate 
    ten = contract(
        ten,
        tensor(U),
        Base.range(2, 1+num_sites),
        Base.range(trans ? 1 : 2, ndims(tensor(U)), step = 2),
        false,
        isconj(ψ),
    )
    ten = permutedim(ten, 2, ndims(ten))
    replacesites!(ψ, ten, site, rev; kwargs...)
end

### Apply a gate to an MPO 
"""
    applygate!(U::AbstractGate, O::MPO, site::Int, [rev::Bool=false]; kwargs...)
    applygate!(O::MPO, U::AbstractGate, site::Int, [rev::Bool=false]; kwargs...)

Apply a circuit gate `U` to an MPO `O`. The bool `rev` specifies the sweeping direction of 
the MPO: `rev=false` for sweeping right, and `rev=true` for sweeping left.
The first (or last for `rev=true`) site that the gate is applied to is `site`.
"""
function applygate!(U::AbstractGate, O::MPO, site::Int, rev::Bool = false; kwargs...)
    _gate_op_validation(
        U,
        O,
        Tuple(
            rev ? Base.range(site+1-length(U), site) : Base.range(site, site+length(U)-1),
        ),
    )
    _applygate!(U, O, site, rev; kwargs...)
end

function applygate!(O::MPO, U::AbstractGate, site::Int, rev::Bool = false; kwargs...)
    _gate_op_validation(
        U,
        O,
        Tuple(
            rev ? Base.range(site+1-length(U), site) : Base.range(site, site+length(U)-1),
        ),
    )
    _applygate!(U, O, site, rev, true; kwargs...)
end

# Unsafe application 
function _applygate!(
    U::AbstractGate,
    O::MPO,
    site::Int,
    rev::Bool = false,
    trans::Bool = false;
    kwargs...,
)
    # Properties 
    num_sites = length(U)
    firstsite = rev ? site - num_sites + 1 : site
    lastsite = rev ? site : site + num_sites - 1

    # Move the center 
    if !(center(O) >= firstsite && center(O) <= lastsite)
        if center(O) < firstsite
            movecenter!(O, firstsite)
        else
            movecenter!(O, lastsite)
        end
    end

    # Contract MPS tensors
    ten = O[firstsite]
    for i in Base.range(firstsite+1, lastsite)
        ten = contract(ten, O[i], ndims(ten), 1)
    end

    # Apply the gate 
    ten = contract(
        ten,
        tensor(U),
        Base.range(trans ⊻ istranspose(O) ? 3 : 2, ndims(ten)-1, step = 2),
        Base.range(trans ? 1 : 2, ndims(tensor(U)), step = 2),
        false,
        isconj(O),
    )

    # Permute the gates 
    perms = (
        1,
        Tuple(
            map(
                j -> isodd(j) ? 2 + num_sites + cld(j, 2) : 1 + cld(j, 2),
                Base.OneTo(2*num_sites),
            ),
        )...,
        2+num_sites,
    )
    ten = _permutedims(ten, perms)
    replacesites!(O, ten, site, rev; kwargs...)
end
