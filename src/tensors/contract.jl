#=
    Functions for contracting tensors 
=#

export contract, contract!

### General contractions 

"""
    contract(x, y, cix, ciy, [conjx=false, conjy=false]; kwargs...)

Contract tensors `x` and `y` across dimensions `cix` and `ciy`, and returns it as `z`.

# Arguments

    - `x`: first tensor to contract.
    - `y': second tensor to contract.
    - `cix`: the dimensions of the first tensor to contract.
    - `ciy`: the dimensions of the second tensor to contract.
    - `conjx::Bool=false`: Take the complex conjugate of argument x?
    - `conjy::Bool=false`: Take the complex conjugate of argument y?

# Optional Keyword Arguments
    
    - `tocache::Bool=true`: store the result in the second level of the cache?
    - `sublevel=:auto`: if stored in cache, at which sublevel? :auto finds non-aliased memory

# Examples 

```jldoctest
julia> x = randn(ComplexF64, 2, 3, 4);
julia> y = randn(ComplexF64, 3, 5, 6);
julia> z = contract(x, y, 2, 1);
julia> size(z)
(2, 4, 5, 6)
```

```jldoctest
julia> x = randn(ComplexF64, 2, 3, 4, 5);
julia> y = randn(ComplexF64, 6, 5, 2, 7);
julia> z = contract(x, y, (1, 4), (3, 2));
julia> size(z)
(3, 4, 6, 7)
```
"""
function contract(
    x,
    y,
    cix,
    ciy,
    conjx::Bool = false,
    conjy::Bool = false;
    tocache::Bool = true,
    sublevel = :auto,
)
    # Dimensions of the problem
    _contract_checkdims(x, y, cix, ciy)
    sx, sy, rix, riy, pix, piy = _contract_permuted_dimensions(x, y, cix, ciy)

    # Create new matrix to store result 
    dims = (_contract_dims(sx, rix)..., _contract_dims(sy, riy)...)
    if tocache
        z = cache(dims, x, y; level = 2, sublevel = sublevel)
    else
        z = promote_tensor(dims, x, y)
    end
    _contract!(z, x, y, sx, sy, cix, ciy, rix, riy, pix, piy, conjx, conjy)
    return z
end

function contract(
    x,
    y,
    cix::Int,
    ciy::Int,
    conjx::Bool = false,
    conjy::Bool = false;
    kwargs...,
)
    return contract(x, y, (cix,), (ciy,), conjx, conjy; kwargs...)
end


"""
    contract!(z, x, y, cix, ciy, [conjx=false, conjy=false])

Contract tensors `x` and `y` across dimensions `cix` and `ciy`, and store the result in `z`.
In-place version of contract.

# Arguments

    - `z`: tensor to store the result.
    - `x`: first tensor to contract.
    - `y': second tensor to contract.
    - `cix`: the dimensions of the first tensor to contract.
    - `ciy`: the dimensions of the second tensor to contract.
    - `conjx::Bool=false`: Take the complex conjugate of argument x?
    - `conjy::Bool=false`: Take the complex conjugate of argument y?

# Examples 

```jldoctest
julia> x = randn(ComplexF64, 2, 3, 4);
julia> y = randn(ComplexF64, 3, 5, 6);
julia> z = similar(x, (2, 4, 5, 6));
julia> contract!(z, x, y, 2, 1)
```
"""
function contract!(z, x, y, cix, ciy, conjx::Bool = false, conjy::Bool = false)
    _contract_checkdims(x, y, cix, ciy)
    sx, sy, rix, riy, pix, piy = _contract_permuted_dimensions(x, y, cix, ciy)
    _contract_checkreturn(z, sx, sy, rix, riy)
    _contract!(z, x, y, sx, sy, cix, ciy, rix, riy, pix, piy, conjx, conjy)
end

function contract!(z, x, y, cix::Int, ciy::Int, conjx::Bool = false, conjy::Bool = false)
    contract!(z, x, y, (cix,), (ciy,), conjx, conjy)
end


### Matrix-vector contractions to reduce overhead 
function contract!(
    z,
    x::S,
    y::T,
    cix::Int = 2,
    ciy::Int = 1,
    conjx::Bool = false,
    conjy::Bool = false,
) where {S<:AbstractArray{<:Number,2},T<:AbstractArray{<:Number,1}}
    _contract_checkdims(x, y, cix, ciy)
    if (ndims(z) != 1) || length(z) != size(x, cix == 1 ? 2 : 1)
        throw(ArgumentError("Destination tensor has the wrong dimensions."))
    end
    if conjx
        x = (cache(eltype(x), size(x), 1) .= conj.(x))
    end
    if conjy
        y =
            cache(eltype(y), size(y), 1, (!conjx || (length(x) != length(y))) ? 1 : 2) .=
                conj.(y)
    end
    mul!(z, cix == 2 ? x : transpose(x), y)
end

function contract!(
    z,
    x::S,
    y::T,
    cix::Int = 1,
    ciy::Int = 1,
    conjx::Bool = false,
    conjy::Bool = false,
) where {S<:AbstractArray{<:Number,1},T<:AbstractArray{<:Number,2}}
    contract!(z, y, x, ciy, cix, conjx, conjy)
end

function contract(
    x::S,
    y::T,
    cix::Int = 2,
    ciy::Int = 1,
    conjx::Bool = false,
    conjy::Bool = false,
) where {S<:AbstractArray{<:Number,2},T<:AbstractArray{<:Number,1}}
    z = zeros(Base.promote_op(*, eltype(x), eltype(y)), size(x, cix == 2 ? 1 : 2))
    contract!(z, x, y, cix, ciy, conjx, conjy)
    return z
end

function contract(
    x::S,
    y::T,
    cix::Int = 1,
    ciy::Int = 1,
    conjx::Bool = false,
    conjy::Bool = false,
) where {S<:AbstractArray{<:Number,1},T<:AbstractArray{<:Number,2}}
    z = zeros(Base.promote_op(*, eltype(x), eltype(y)), size(y, ciy == 2 ? 1 : 2))
    contract!(z, x, y, cix, ciy, conjx, conjy)
    return z
end

### Matrix-matrix multiplications to reduce overhead
function contract!(
    z,
    x::S,
    y::T,
    cix::Int = 2,
    ciy::Int = 1,
    conjx::Bool = false,
    conjy::Bool = false,
) where {S<:AbstractArray{<:Number,2},T<:AbstractArray{<:Number,2}}
    _contract_checkdims(x, y, cix, ciy)
    if (ndims(z) != 2) ||
       (size(z, 1) != size(x, cix == 1 ? 2 : 1)) ||
       (size(z, 2) != size(y, ciy == 1 ? 2 : 1))
        throw(ArgumentError("Destination tensor has the wrong dimensions."))
    end
    if conjx
        x = (cache(x, 1, 1) .= conj.(x))
    end
    if conjy
        y = cache(y, 1, (!conjx || (length(x) != length(y))) ? 1 : 2) .= conj.(y)
    end
    mul!(z, cix == 2 ? x : transpose(x), ciy == 1 ? y : transpose(y))
end

function contract(
    x::S,
    y::T,
    cix::Int = 2,
    ciy::Int = 1,
    conjx::Bool = false,
    conjy::Bool = false,
) where {S<:AbstractArray{<:Number,2},T<:AbstractArray{<:Number,2}}
    z = zeros(
        Base.promote_op(*, eltype(x), eltype(y)),
        size(x, cix == 2 ? 1 : 2),
        size(y, ciy == 2 ? 1 : 2),
    )
    contract!(z, x, y, cix, ciy, conjx, conjy)
    return z
end


### Unsafe contractions
# Contraction with permutation 
function _contract!(z, x, y, sx, sy, cix, ciy, rix, riy, pix, piy, conjx, conjy)
    # Permute tensors
    px = permutedims!(cache(x, _contract_dims(sx, pix)), x, pix)
    py = permutedims!(
        cache(y, _contract_dims(sy, piy), 1, length(x) == length(y) ? 2 : 1),
        y,
        piy,
    )

    # Conjugations 
    if conjx
        px .= conj.(px)
    end
    if conjy
        py .= conj.(py)
    end

    # Contract
    mul!(
        reshape(z, (prod(_contract_dims(sx, rix)), prod(_contract_dims(sy, riy)))),
        reshape(px, (prod(_contract_dims(sx, rix)), prod(_contract_dims(sx, cix)))),
        reshape(py, (prod(_contract_dims(sy, ciy)), prod(_contract_dims(sy, riy)))),
    )
end


### Calculates all the needed dimensions and permutations for contractions 
function _contract_permuted_dimensions(x, y, cix, ciy)
    #=
        the setdiff function returns as a vector and thus makes allocations...
        can we do this in a different way which returns a tuple?
    =#
    # Find the dimensions of the tensors 
    sx = size(x)
    sy = size(y)

    # Find the indices which aren't contracted & permutation indices 
    rix = tuple(setdiff(1:length(sx), cix)...)
    riy = tuple(setdiff(1:length(sy), ciy)...)

    # Find the permutation for the tensors
    pix = (rix..., cix...)
    piy = (ciy..., riy...)

    return sx, sy, rix, riy, pix, piy
end

function _contract_dims(s, i)
    return map(i->s[i], i)
end

### Checks 
# Check to see if contraction tensors are the right size
function _contract_checkdims(x, y, cix, ciy)
    #if typeof(get_backend(x)) != typeof(get_backend(y))
    #    throw(ArgumentError("Tensors are stored on different backends."))
    #end
    if length(cix) != length(ciy)
        throw(ArgumentError("Contraction indices have different lengths."))
    end
    for i in eachindex(cix)
        if cix[i] < 1 || cix[i] > ndims(x) || ciy[i] < 1 || ciy[i] > ndims(y)
            throw(ArgumentError("Contraction indices are out of range."))
        end
        if size(x, cix[i]) != size(y, ciy[i])
            throw(ArgumentError("Contraction indices do not match."))
        end
    end
end

# Check to see if the tensor the contraction is stored in is correct
function _contract_checkreturn(z, sx, sy, rix, riy)
    if (size(z) != (_contract_dims(sx, rix)..., _contract_dims(sy, riy)...))
        throw(ArgumentError("Destination tensor has the wrong dimensions."))
    end
end
