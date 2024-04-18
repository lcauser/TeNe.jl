# Abstract MPS
abstract type AbstractMPS end
"""
    issimilar(::AbstractMPS...)

Check to see if MPSs share the same properties.
"""
function issimilar(ψs::AbstractMPS...)
    for i = Base.range(2, length(ψs))
        length(ψs[i]) != length(ψs[1]) && return false 
        if dim(ψs[1]) == 0 || dim(ψs[i]) == 0
            for j = 1:length(ψ)
                ψ_dims = map(k -> size(ψs[1][j], k), Base.range(2, 1+rank(ψs[1])))
                ϕ_dims = map(k -> size(ψs[i][j], k), Base.range(2, 1+rank(ψs[i])))
                ψ_dims != ϕ_dims && return false
            end
        elseif dim(ψs[1]) != dim(ψs[i])
            return false 
        end
    end
    return true
end
export issimilar

# Traits for GMPS
abstract type GMPSTrait <: AbstractMPS end
Base.eltype(ψ::GMPSTrait) = Base.eltype(ψ.MPS)
Base.length(ψ::GMPSTrait) = length(ψ.MPS)
Base.getindex(ψ::GMPSTrait, i::Int) = ψ.MPS[i]
Base.firstindex(ψ::GMPSTrait) = Base.firstindex(ψ.MPS)
Base.lastindex(ψ::GMPSTrait) = Base.lastindex(ψ.MPS)
Base.eachindex(ψ::GMPSTrait) = Base.eachindex(ψ.MPS)
Base.setindex!(ψ::GMPSTrait, x, i::Int) = Base.setindex!(ψ.MPS, x, i)
dim(ψ::GMPSTrait) = dim(ψ.MPS)
rank(ψ::GMPSTrait) = rank(ψ.MPS)
center(ψ::GMPSTrait) = center(ψ.MPS)
movecenter!(ψ::GMPSTrait, site::Int) = movecenter!(ψ.MPS, site)
bonddim(ψ::GMPSTrait, site::Int) = bonddim(ψ.MPS, site)
maxbonddim(ψ::GMPSTrait) = maxbonddim(ψ.MPS)
TeNe.norm(ψ::GMPSTrait) = TeNe.norm(ψ.MPS)
