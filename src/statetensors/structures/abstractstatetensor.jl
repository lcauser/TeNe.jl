abstract type AbstractStateTensor <: TensorNetworkState end

# Maybe for later: genalise to more general tensor network states?
"""
    issimilar(::AbstractStateTensor...)

Check to see if state tensors share the same properties.
"""
function issimilar(ψs::AbstractStateTensor...)
    for i in Base.range(2, length(ψs))
        length(ψs[i]) != length(ψs[1]) && return false
        dim(ψs[1]) != dim(ψs[i]) && return false
    end
    return true
end
export issimilar

# Traits for State tensors 
abstract type GStateTensorTrait <: AbstractStateTensor end
Base.eltype(ψ::GStateTensorTrait) = Base.eltype(ψ.StateTensor)
Base.length(ψ::GStateTensorTrait) = length(ψ.StateTensor)
Base.getindex(ψ::GStateTensorTrait, is...) = ψ.StateTensor[i...]
Base.firstindex(ψ::GStateTensorTrait) = Base.firstindex(ψ.StateTensor)
Base.lastindex(ψ::GStateTensorTrait) = Base.lastindex(ψ.StateTensor)
Base.eachindex(ψ::GStateTensorTrait) = Base.eachindex(ψ.StateTensor)
Base.setindex!(ψ::GStateTensorTrait, x, is...) = Base.setindex!(ψ.StateTensor, x, is...)
dim(ψ::GStateTensorTrait) = dim(ψ.StateTensor)
dim(ψ::GStateTensorTrait, i...) = dim(ψ.StateTensor, i...)
dims(ψ::GStateTensorTrait, i::Int) = dims(ψ.StateTensor, i)
TeNe.rank(ψ::GStateTensorTrait) = rank(ψ.StateTensor)
TeNe.norm(ψ::GStateTensorTrait) = TeNe.norm(ψ.StateTensor)
tensor(ψ::GStateTensorTrait) = tensor(ψ.StateTensor)
