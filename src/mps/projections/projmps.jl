#=
    Projections for inner products of MPS
=#

mutable struct ProjMPS{Q<:Number} <: MPSProjection
    objects::Vector{Union{MPS,MPO}}
    lefts::Vector{<:AbstractArray}
    rights::Vector{<:AbstractArray}
    center::Int
    sym::Bool
    λ::Q
end

export ProjMPS
"""
    ProjMPS(ψ::MPS, args::Union{MPS, MPO}...; kwargs...)

Create a projection for the inner product of a string of GMPSs.

# Optional Keyword Arguments

    - `λ::Number=1.0`: Multiplies the projection by a constant.
    - `center::Int=0`: Initalise the projection at some center. By default, it defaults to 
    the center of the last MPS (`center=0`).
    - `sym::Bool=false`: Is the projection symmetric in `ψ`? i.e. of the form `<ψ|O1…Ok|ψ>`?
"""
function ProjMPS(
    ψ::MPS,
    args::Union{MPS,MPO}...;
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

    # Create blocks 
    T = _promote_tensor_eltype(ψ, args...)
    lefts = Array{T}[]
    rights = Array{T}[]
    for i in Base.range(firstindex(ψ), lastindex(ψ)-1)
        push!(lefts, zeros(T, size(ψ[i], 1), map(arg->size(arg[i+1], 1), args)...))
        push!(rights, zeros(T, size(ψ[end+1-i], 1), map(arg->size(arg[end-i], 1), args)...))
    end

    # Build blocks 
    center = center == 0 ? ψ.center : center
    projψ = ProjMPS(Union{MPS,MPO}[ψ, args...], lefts, rights, center, sym, λ)
    for i in Base.range(firstindex(ψ), lastindex(ψ)-1)
        buildleft!(projψ, i)
        buildright!(projψ, lastindex(ψ)+1-i)
    end

    # Move center 
    projψ.center = center == 0 ? args[end].center : center
    return projψ
end

### Building the blocks from the left 
function buildleft!(projψ::ProjMPS, site::Int)
    # Build on the previous left block
    ten = _buildleft(projψ, site)

    # Save to memory
    if size(ten) == size(leftblock(projψ, site))
        projψ.lefts[site] .= ten
    else
        projψ.lefts[site] = copy(ten)
    end
end

function _buildleft(projψ::ProjMPS, site::Int)
    ten = leftblock(projψ, site-1)
    ten = contract(
        ten,
        projψ.objects[begin][site],
        1,
        1,
        false,
        !isconj(projψ.objects[begin]),
    )
    for i in Base.OneTo(length(projψ.objects)-2)
        ten = contract(
            ten,
            projψ.objects[begin+i][site],
            (1, ndims(ten)-1),
            (1, innerind(projψ.objects[begin+i])),
            false,
            isconj(projψ.objects[begin+i]),
        )
    end
    ten = contract(
        ten,
        projψ.objects[end][site],
        (1, ndims(ten)-1),
        (1, 2),
        false,
        isconj(projψ.objects[end]),
    )
    return ten
end

### Build blocks from the right
function buildright!(projψ::ProjMPS, site::Int)
    # Build on the previous right block
    ten = _buildright(projψ, site)

    # Save to memory
    if size(ten) == size(rightblock(projψ, site))
        projψ.rights[site-1] .= ten
    else
        projψ.rights[site-1] = copy(ten)
    end
end

function _buildright(projψ::ProjMPS, site::Int)
    ten = rightblock(projψ, site+1)
    ten = contract(
        projψ.objects[end][site],
        ten,
        3,
        length(projψ.objects),
        isconj(projψ.objects[end]),
        false,
    )
    for i in Base.OneTo(length(projψ.objects)-2)
        ten = contract(
            projψ.objects[end-i][site],
            ten,
            (outerind(projψ.objects[end-i]), 4),
            (2, ndims(ten)),
            isconj(projψ.objects[end-i]),
            false,
        )
    end
    ten = contract(
        projψ.objects[begin][site],
        ten,
        (2, 3),
        (2, ndims(ten)),
        !isconj(projψ.objects[begin]),
        false,
    )
    return ten
end

### Inner products
export inner
"""
    inner(projψ::ProjMPS)

Calculate the inner product using the projection around the canonical center.
"""
function inner(projψ::ProjMPS)
    left = _buildleft(projψ, center(projψ))
    right = rightblock(projψ, center(projψ)+1)
    return projψ.λ *
           contract(left, right, Base.OneTo(ndims(left)), Base.OneTo(ndims(right)))[]
end


"""
    inner(projψ::ProjMPS, A::AbstractArray, dir::Bool=false)

Calculate the inner product using the projection around the canonical center,
with some replacement for the lattice sites `A`. Use `dir=true` if sweeping left.
"""
function inner(projψ::ProjMPS, A::AbstractArray, dir::Bool = false)
    grad = gradient(projψ, A, dir)
    return contract(
        grad,
        A,
        Base.OneTo(ndims(grad)),
        Base.OneTo(ndims(A)),
        false,
        isconj(projψ.objects[end]),
    )[]
end


### Calculating gradients with respect to a tensor 
export gradient
"""
    gradient(projψ::ProjMPS, A::AbstractArray, dir::Bool=false)

Calculate the gradient of an inner product around the center of a projection
with some replacement for the lattice sites `A`. Use `dir=true` if sweeping left.
"""
function gradient(projψ::ProjMPS, A::AbstractArray, dir::Bool = false)
    # Properties...
    nsites = ndims(A)-2
    firstsite = dir ? center(projψ) + 1 - nsites : center(projψ)
    lastsite = dir ? center(projψ) : center(projψ) + nsites - 1

    # Fetch the blocks 
    left = leftblock(projψ, firstsite-1)
    right = rightblock(projψ, lastsite+1)

    # Contract with the first MPS
    if projψ.sym
        left = contract(left, A, 1, 1, false, !isconj(projψ.objects[begin]))
        left = permutedim(left, ndims(left), 1)
    else
        for i in Base.range(firstsite, lastsite)
            left = contract(
                left,
                projψ.objects[begin][i],
                1,
                1,
                false,
                !isconj(projψ.objects[begin]),
            )
            left = permutedim(left, ndims(left), 1)
        end
    end

    # Contract with MPOs
    for i in Base.OneTo(length(projψ.objects)-2)
        for j in Base.range(firstsite, lastsite)
            left = contract(
                left,
                projψ.objects[begin+i][j],
                (1+i, length(projψ.objects)+1),
                (1, innerind(projψ.objects[begin+i])),
                false,
                isconj(projψ.objects[begin+i]),
            )
            left = permutedim(left, ndims(left), 1+i)
        end
    end

    # Contract with the right
    grad = contract(
        left,
        right,
        Base.OneTo(length(projψ.objects)-1),
        Base.OneTo(length(projψ.objects)-1),
    )
    grad .*= projψ.λ
    return grad
end

# Calculating the gradient on the conjugate MPS!
function gradientconj(
    projψ::ProjMPS,
    A::AbstractArray,
    dir::Bool = false;
    tocache::Bool = false,
)
    # Properties...
    nsites = ndims(A)-2
    firstsite = dir ? center(projψ) + 1 - nsites : center(projψ)
    lastsite = dir ? center(projψ) : center(projψ) + nsites - 1

    # Fetch the blocks 
    left = leftblock(projψ, firstsite-1)
    right = rightblock(projψ, lastsite+1)

    # Contract with the first MPS
    if projψ.sym
        right = contract(A, right, ndims(A), ndims(right), isconj(projψ.objects[end]))
        right = permutedim(right, 1, ndims(right))
    else
        for i in Base.range(lastsite, firstsite, step = -1)
            right = contract(
                projψ.objects[end][i],
                right,
                3,
                ndims(right),
                isconj(projψ.objects[end]),
            )
            right = permutedim(right, 1, ndims(right))
        end
    end

    # Contract with MPOs 
    for i in Base.OneTo(length(projψ.objects)-2)
        for j in Base.range(lastsite, firstsite, step = -1)
            right = contract(
                projψ.objects[end-i][j],
                right,
                (outerind(projψ.objects[end-i]), 4),
                (nsites, ndims(right)-i),
                isconj(projψ.objects[end-i]),
            )
            right = permutedim(right, 1, ndims(right)-i)
        end
    end

    # Contract the left with the right
    grad = contract(
        left,
        right,
        Base.range(2, ndims(left)),
        Base.range(2+nsites, ndims(right));
        tocache = tocache,
    )
    grad .*= projψ.λ
    return grad
end


export product

"""
    product(projψ::ProjMPS, A::AbstractArray, dir::Bool=false; kwargs...)

Calculate the product of an MPS projection with a tensor `A`.

# Optional Keyword Arguments 

    - `tocache::Bool=false`: Save the product to the cache?
"""
function product(projψ::ProjMPS, A::AbstractArray, dir::Bool = false; tocache::Bool = true)
    # If not symmetric, this is just the inner product!
    if !projψ.sym
        return inner(projψ, A, dir)
    end

    # Properties...
    nsites = ndims(A)-2
    nobjects = length(projψ.objects)
    firstsite = dir ? center(projψ) + 1 - nsites : center(projψ)
    lastsite = dir ? center(projψ) : center(projψ) + nsites - 1

    # Fetch the blocks 
    left = leftblock(projψ, firstsite-1)
    right = rightblock(projψ, lastsite+1)

    # Contract with A 
    right = contract(A, right, 2+nsites, nobjects, isconj(projψ.objects[end]))
    right = permutedim(right, 1, nsites+nobjects)

    # Contract with the middle blocks 
    for i in Base.OneTo(nobjects-2)
        for j in Base.range(lastsite, firstsite, step = -1)
            right = contract(
                projψ.objects[end-i][j],
                right,
                (outerind(projψ.objects[end-i]), 4),
                (nsites, nsites+nobjects-i),
                isconj(projψ.objects[end-i]),
            )
            right = permutedim(right, 1, nsites+nobjects-i)
        end
    end

    # Contract with the left block 
    prod = contract(
        left,
        right,
        Base.range(2, ndims(left)),
        Base.range(nsites+2, ndims(right));
        tocache = tocache,
    )
    prod .*= projψ.λ
    return prod
end
