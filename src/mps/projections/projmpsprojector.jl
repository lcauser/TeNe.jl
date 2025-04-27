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
"""
    ProjMPS(ψ::MPS, args::Union{MPS, MPO, MPSProjector}...; kwargs...)

Create a projection for the inner product of a string of GMPSs, including MPSProjectors.

# Optional Keyword Arguments

    - `λ::Number=1.0`: Multiplies the projection by a constant.
    - `center::Int=0`: Initalise the projection at some center. By default, it defaults to 
    the center of the MPS `ψ` (`center=0`).
    - `sym::Bool=false`: Is the projection symmetric in `ψ`? i.e. of the form `<ψ|O1…Ok|ψ>`?
"""
function ProjMPS(
    ψ::MPS,
    args::Union{MPS,MPO,MPSProjector}...;
    λ::Number = 1.0,
    center::Int = 0,
    sym::Bool = false,
)
    # Validation 
    _inner_validation(ψ, args...)
    _vec_op_vec_validation(ψ, args[end], args[begin:(end-1)]...)
    if center < 0 || center > length(ψ)
        throw(DomainError("Center $(center) is out of range 0 and $(length(ψ))."))
    end
    center = center == 0 ? args[end].center : center

    # Create the projections
    projections = MPSProjection[]
    projection = Union{MPS,MPO}[ψ]
    for arg in args
        if ismpsprojector(arg)
            push!(projection, arg.ψ)
            projection = ProjMPS(projection...; center = center)
            push!(projections, projection)
            projection = Union{MPS,MPO}[arg.ϕ]
            λ *= arg.λ
        else
            push!(projection, arg)
        end
    end
    projection = ProjMPS(projection...; center = center)
    push!(projections, projection)

    # Create the projector 
    projψ = ProjMPSProjectors(projections, center, sym, λ)
    return projψ
end

### Length
Base.length(proj::ProjMPSProjectors) = length(proj.projections[begin])

### Moving the canonical center 
export movecenter!
function movecenter!(projψ::ProjMPSProjectors, idx::Int)
    for proj in projψ.projections
        movecenter!(proj, idx)
    end
    projψ.center = idx
end

### Taking the inner product
export inner
function inner(projψ::ProjMPSProjectors)
    prod = projψ.λ
    for i in eachindex(projψ.projections)
        prod *= inner(projψ.projections[i])
    end
    return prod
end

function inner(projψ::ProjMPSProjectors, A::AbstractArray, dir::Bool = false)
    prod = projψ.λ
    # Contract with first
    if projψ.sym
        prod *= contract(
            A,
            gradientconj(projψ.projections[1], A, dir),
            Base.OneTo(ndims(A)),
            Base.OneTo(ndims(A)),
            !isconj(projψ.projections[1]),
        )[]
    else
        prod *= inner(projψ.projections[1])
    end

    # Contract with middles 
    for i in Base.range(2, length(projψ.projections)-1)
        prod *= inner(projψ.projections[i])
    end

    # Contract with the last 
    prod *= inner(projψ.projections[end], A, dir)
    return prod
end

### Gradients 
export gradient
function gradient(projψ::ProjMPSProjectors, A::AbstractArray, dir::Bool = false)
    g = gradient(projψ.projections[end], A, dir)
    for i in Base.range(firstindex(projψ.projections), lastindex(projψ.projections)-1)
        g .*= inner(projψ.projections[i])
    end
    return g
end

### Products with the tensor 
function product(
    projψ::ProjMPSProjectors,
    A::AbstractArray,
    dir::Bool = false;
    tocache::Bool = true,
)
    # If not symmetric, this is just the inner product!
    if !projψ.sym
        return inner(projψ, A, dir)
    end

    # Is symmetric 
    prod = gradientconj(projψ.projections[begin], A, dir; tocache = tocache)
    for i in Base.range(firstindex(projψ.projections)+1, lastindex(projψ.projections)-1)
        prod .*= inner(projψ.projections[i])
    end
    prod .*= inner(projψ.projections[end], A, dir)
    prod .*= projψ.λ
    return prod
end
