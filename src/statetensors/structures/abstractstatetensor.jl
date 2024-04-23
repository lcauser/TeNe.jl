abstract type AbstractStateTensor end 

# Maybe for later: genalise to more general tensor network states?
"""
    issimilar(::AbstractStateTensor...)

Check to see if state tensors share the same properties.
"""
function issimilar(ψs::AbstractStateTensor...)
    for i = Base.range(2, length(ψs))
        length(ψs[i]) != length(ψs[1]) && return false 
        dim(ψs[1]) != dim(ψs[i]) && return false
    end
    return true
end
export issimilar

# Traits for State tensors 
abstract type GStateTensorTrait <: AbstractStateTensor end
ase.eltype(ψ::GStateTensorTrait) = Base.eltype(ψ.StateTensor)
Base.length(ψ::GStateTensorTrait) = length(ψ.StateTensor)
dim(ψ::GStateTensorTrait) = dim(ψ.StateTensor)
TeNe.rank(ψ::GStateTensorTrait) = rank(ψ.StateTensor)
TeNe.norm(ψ::GStateTensorTrait) = TeNe.norm(ψ.StateTensor)
