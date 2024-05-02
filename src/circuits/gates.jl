abstract type AbstractGate end

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
    uis = Base.range(2, ndims(tensor(U)), step=2)

    # Do the contraction
    ψ′ = contract(tensor(U), tensor(ψ), uis, sites, isconj(ψ), false)
    permutedims!(tensor(ψ), ψ′, reverseperms((sites..., ris...)))
end

function _applygate!(ψ::StateVector, U::AbstractGate, sites)
    # Contraction indices
    ris = Tuple(setdiff(Base.OneTo(length(ψ)), sites)) # Remove allocations?
    uis = Base.range(1, ndims(tensor(U)), step=2)

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
    uis = Base.range(istranspose(O) ? 1 : 2, ndims(tensor(U)), step=2)

    # Do the contraction 
    O′ = contract(tensor(U), tensor(O), uis, Ois, isconj(O), false)
    permutedims!(tensor(O), O′, reverseperms((Ois..., ris...)))
end

function _applygate!(O::StateOperator, U::AbstractGate, sites)
    # Contraction indices
    Ois = Tuple(map(j->istranspose(O) ? innerind(O, j) : outerind(O, j), sites))
    ris = Tuple(setdiff(Base.OneTo(ndims(tensor(O))), Ois))
    uis = Base.range(istranspose(O) ? 2 : 1, ndims(tensor(U)), step=2)

    # Do the contraction 
    O′ = contract(tensor(O), tensor(U), Ois, uis, false, isconj(O))
    permutedims!(tensor(O), O′, reverseperms((ris..., Ois...)))
end


### Applying gates to MPS 
function applygate!(U::AbstractGate, ψ::MPS, site::Int, rev::Bool=false; kwargs...)

end

# Unsafe application 
function _applygate!(U::AbstractGate, ψ::MPS, site::Int, rev::Bool=false; kwargs...)
    num_sites = length(U)
    firstsite = rev ? site - num_sites + 1 : site 
    lastsite = rev ? site : site + num_sites - 1

    # Contract MPS tensors
    ten = ψ[firstsite]
    for i = Base.range(firstsite+1, lastsite)
        ten = contract(ten, ψ[i], ndims(ten), 1)
    end

    # Apply the gate 
    ten = contract(ten, tensor(U), Base.range(2, 1+num_sites), Base.range(2, ndims(tensor(U)), step=2))
    ten = permutedim(ten, 2, ndims(ten))
    replacesites!(ψ, ten, site, rev; kwargs...)
end



### A customizable gate 
mutable struct CircuitGate{d, n, T} <: AbstractGate where {d, n, T<:AbstractArray}
    gate::T
end
export creategate
"""
    creategate(gate::AbstractArray)

Create a gate using a tensor.
"""
function creategate(gate::AbstractArray)
    # Validation 
    num_qubits = fld(ndims(gate), 2)
    if ndims(gate) % 2 != 0
        throw(ArgumentError("The gates must have an even number of dimensions."))
    end
    sz = size(gate)
    if !all(j->j==sz[1], sz)
        throw(ArgumentError("The gate must have the same dimensions sizes."))
    end
    # 
    return CircuitGate{sz[1], num_qubits, typeof(gate)}(gate)
end

# Properties 
export tensor, dim, qubits
tensor(gate::CircuitGate) = gate.gate
dim(::CircuitGate{d}) where {d} = d 
Base.length(::CircuitGate{d, n}) where {d, n} = n