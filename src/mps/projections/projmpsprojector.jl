#=
    Dealing with MPS Projections with projectors
=#

mutable struct ProjMPSProjectors{Q<:Number} <: MPSProjection
    projections::Vector{<:MPSProjection}
    center::Int
    sym::Bool
    λ::Q
end

export ProjMPS
function ProjMPS(ψ::MPS, args::Union{MPS, MPO, MPSProjector}...; λ::Number=1.0, center::Int=0, sym::Bool=false)
    # Validation 
    _inner_validation(ψ, args...)
    _vec_op_vec_validation(ψ, args[end], args[begin:end-1]...)
    if center < 0 || center > length(ψ)
        throw(DomainError("Center $(center) is out of range 0 and $(length(ψ))."))
    end

    # Create the projections


    # Move center 
    projψ.center = center == 0 ? ψ.center : center
    return projψ
end