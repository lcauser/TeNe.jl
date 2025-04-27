#=
    Dealing with MPS Projections with projectors
=#

mutable struct ProjMPSSquared{Q<:Number} <: MPSProjection
    ProjMPS::ProjMPS
    λ::Q
end

export ProjMPSSquared
"""
    ProjMPS(ψ::MPS, args::Union{MPS, MPO, MPSProjector}...; kwargs...)

Create a projection for the inner product of a string of GMPSs, including MPSSquared.

# Optional Keyword Arguments

    - `λ::Number=1.0`: Multiplies the projection by a constant.
    - `center::Int=0`: Initalise the projection at some center. By default, it defaults to 
    the center of the MPS `ψ` (`center=0`).
"""
function ProjMPSSquared(ψ::MPS, args::MPS...; λ::Number = 1.0, center::Int = 0)
    # Validation 
    _inner_validation(ψ, args...)
    _vec_op_vec_validation(ψ, args...)
    if center < 0 || center > length(ψ)
        throw(DomainError("Center $(center) is out of range 0 and $(length(ψ))."))
    end
    center = center == 0 ? args[end].center : center
    return ProjMPSSquared(ProjMPS(ψ, args...), λ)
end

### Length
Base.length(proj::ProjMPSSquared) = length(proj.ProjMPS)

### Moving the canonical center 
center(projψ::ProjMPSSquared) = center(projψ.ProjMPS)
movecenter!(projψ::ProjMPSSquared, idx::Int) = movecenter!(projψ.ProjMPS, idx)

### Taking the inner product
function inner(projψ::ProjMPSSquared)
    return projψ.λ * abs(inner(projψ.ProjMPS))^2
end

function inner(projψ::ProjMPSSquared, A::AbstractArray, dir::Bool = false)
    return projψ.λ * abs(inner(projψ.ProjMPS, A, dir))^2
end

### Gradients 
export gradient
function gradient(projψ::ProjMPSSquared, A::AbstractArray, dir::Bool = false)
    g = gradient(projψ.ProjMPS, A, dir)
    g .*= conj(inner(projψ.ProjMPS, A, dir)) * projψ.λ
    return g
end

### Products with the tensor 
function product(
    projψ::ProjMPSSquared,
    A::AbstractArray,
    dir::Bool = false;
    tocache::Bool = true,
)
    g = gradient(projψ, A, dir)
    g .= conj.(g)
    if !tocache
        g = copy(g)
    end
    return g
end
