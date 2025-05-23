# Abstract MPS
abstract type AbstractMPS <: TensorNetworkState end
"""
    issimilar(::AbstractMPS...)

Check to see if MPSs share the same properties.
"""
function issimilar(ψs::AbstractMPS...)
    for i in Base.range(2, length(ψs))
        length(ψs[i]) != length(ψs[1]) && return false
        dim(ψs[1]) != dim(ψs[i]) && return false
    end
    return true
end
export issimilar

# Traits for GMPS
abstract type GMPSTrait <: AbstractMPS end
Base.eltype(ψ::GMPSTrait) = eltype(ψ.MPS)
Base.length(ψ::GMPSTrait) = length(ψ.MPS)
Base.getindex(ψ::GMPSTrait, i::Int) = ψ.MPS[i]
Base.firstindex(ψ::GMPSTrait) = firstindex(ψ.MPS)
Base.lastindex(ψ::GMPSTrait) = lastindex(ψ.MPS)
Base.eachindex(ψ::GMPSTrait) = eachindex(ψ.MPS)
Base.setindex!(ψ::GMPSTrait, x, i::Int) = setindex!(ψ.MPS, x, i)
dim(ψ::GMPSTrait) = dim(ψ.MPS)
dim(ψ::GMPSTrait, i...) = dim(ψ.MPS, i...)
dims(ψ::GMPSTrait, i::Int) = dims(ψ.MPS, i::Int)
TeNe.rank(ψ::GMPSTrait) = rank(ψ.MPS)
center(ψ::GMPSTrait) = center(ψ.MPS)
movecenter!(ψ::GMPSTrait, site::Int) = movecenter!(ψ.MPS, site)
bonddim(ψ::GMPSTrait, site::Int) = bonddim(ψ.MPS, site)
maxbonddim(ψ::GMPSTrait) = maxbonddim(ψ.MPS)
TeNe.norm(ψ::GMPSTrait) = TeNe.norm(ψ.MPS)

# Default properties to false, but override later for true...
isconj(ψ) = false
istranspose(ψ) = false
