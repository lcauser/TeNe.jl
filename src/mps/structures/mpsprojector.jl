#=
    An MPSProjector allows us to construct an MPO-like structure from the outer product 
    of two MPSs, O = |ψ><ϕ|.
=#

mutable struct MPSProjector{Q<:Number}
    ψ::MPS
    ϕ::MPS
    equal::Bool
    λ::Q
end
export MPSProjector

"""
    MPSProjector(ψ::MPS, ϕ::MPS)
    MPSProjector(ψ::MPS)

Create the projection O = |ψ><ϕ| from two MPSs `ψ` and `ϕ`.

The key word argument λ can be specified to give O = λ|ψ><ϕ|
"""
function MPSProjector(ψ::MPS, ϕ::MPS; λ::Number=1.0)
    if length(ψ) != length(ϕ)
        throw(ArgumentError("MPSs must have the same length."))
    end
    return MPSProjector(ψ, ϕ, false, λ)
end
MPSProjector(ψ::MPS; λ::Number=1.0) = MPSProjector(ψ, ψ, true, λ)


### Properties of an MPSProjector 
Base.firstindex(O::MPSProjector) = Base.firstindex(O.ψ)
Base.lastindex(O::MPSProjector) = Base.lastindex(O.ψ)
Base.eachindex(O::MPSProjector) = Base.eachindex(O.ψ)
TeNe.rank(O::MPSProjector) = 2
TeNe.norm(O::MPSProjector) = norm(O.ψ) * norm(O.ϕ)
Base.eltype(O::MPSProjector) = _promote_tensor_eltype(O.ψ, O.ϕ)
Base.length(O::MPSProjector) = length(O.ψ)
function dims(O::MPSProjector, which::Int)
    if which == 1
        return dims(O.ψ)
    elseif which == 2
        return dims(O.ϕ)
    else
        throw(ArgumentError("Invalid axis for rank-2 tensor."))
    end
end
function dim(O::MPSProjector, which::Int, i::Int)
    if which == 1
        return dim(O.ψ, i)
    elseif which == 2
        return dim(O.ϕ, i)
    else
        throw(ArgumentError("Invalid axis for rank-2 tensor."))
    end
end
innerdim(O::MPSProjector, i::Int) = dim(O.ψ, i)
outerdim(O::MPSProjector, i::Int) = dim(O.ϕ, i)
innerdims(O::MPSProjector) = dims(O.ψ)
outerdims(O::MPSProjector) = dims(O.ϕ)

### Manipulations of MPSProjector 
TeNe.conj(O::MPSProjector) = MPSProjector(conj(O.ψ), conj(O.ϕ), O.equal, conj(O.λ))
TeNe.transpose(O::MPSProjector) = MPSProjector(conj(O.ϕ), conj(O.ψ), O.equal, O.λ)
TeNe.adjoint(O::MPSProjector) = MPSProjector(O.ϕ, O.ψ, O.equal, conj(O.λ))

import Base: *, \
function *(a::Number, O::MPSProjector)
    return MPSProjector(O.ψ, O.ϕ, O.equal, a * O.λ)
end
*(O::MPSProjector, a::Number) = *(a, O)
/(O::MPSProjector, a::Number) = *(1/a, O)