#=
    Functions for contracting tensors 
=#

### General contractions 
function contract(x, y, cix, ciy, conjx=false, conjy=false)
    _contract_checkdims!(x, y, cix, ciy)
    sx, sy, rix, riy, pix, piy = _contract_permuted_dimensions(x, y, cix, ciy)
    z = zeros(Base.promote_op(*, eltype(x), eltype(y)), _contract_dims(sx, rix)..., _contract_dims(sy, riy)...)
    _contract!(z, x, y, sx, sy, cix, ciy, rix, riy, pix, piy, conjx, conjy)  
    return z
end
contract(x, y, cix::Int, ciy::Int, conjx::Bool=false, conjy::Bool=false) = contract(x, y, (cix,), (ciy,), conjx, conjy)

function contract!(z, x, y, cix, ciy, conjx::Bool=false, conjy::Bool=false)
    _contract_checkdims!(x, y, cix, ciy)
    sx, sy, rix, riy, pix, piy= _contract_permuted_dimensions(x, y, cix, ciy)
    _contract_checkreturn!(z, sx, sy, rix, riy)
    _contract!(z, x, y, sx, sy, cix, ciy, rix, riy, pix, piy, conjx, conjy)  
end
contract!(z, x, y, cix::Int, ciy::Int, conjx::Bool=false, conjy::Bool=false) = contract!(z, x, y, (cix,), (ciy,), conjx, conjy)


### Matrix-vector contractions to reduce overhead 
function contract!(z, x::S, y::T, cix::Int=2, ciy::Int=1, conjx::Bool=false, conjy::Bool=false) where {S<:AbstractArray{<:Number, 2}, T<:AbstractArray{<:Number, 1}}
    _contract_checkdims!(x, y, cix, ciy)
    if (ndims(z) != 1) || length(z) != size(x, cix == 1 ? 2 : 1)
        throw(ArgumentError("Destination tensor has the wrong dimensions."))
    end
    if conjx x = (cache(eltype(x), size(x), 1) .= conj.(x)) end
    if conjy y = cache(eltype(y), size(y), (!conjx || (length(x) != length(y))) ? 1 : 2) .= conj.(y) end
    mul!(z, cix == 2 ? x : transpose(x), y)
end
function contract!(z, x::S, y::T, cix::Int=1, ciy::Int=1, conjx::Bool=false, conjy::Bool=false) where {S<:AbstractArray{<:Number, 1}, T<:AbstractArray{<:Number, 2}}
    contract!(z, y, x, ciy, cix, conjx, conjy)
end
function contract(x::S, y::T, cix::Int=2, ciy::Int=1, conjx::Bool=false, conjy::Bool=false) where {S<:AbstractArray{<:Number, 2}, T<:AbstractArray{<:Number, 1}}
    z = zeros(Base.promote_op(*, eltype(x), eltype(y)), size(x, cix == 2 ? 1 : 2))
    contract!(z, x, y, cix, ciy, conjx, conjy)
    return z
end
function contract(x::S, y::T, cix::Int=1, ciy::Int=1, conjx::Bool=false, conjy::Bool=false) where {S<:AbstractArray{<:Number, 1}, T<:AbstractArray{<:Number, 2}}
    z = zeros(Base.promote_op(*, eltype(x), eltype(y)), size(y, ciy == 2 ? 1 : 2))
    contract!(z, x, y, cix, ciy, conjx, conjy)
    return z
end

### Matrix-matrix multiplications to reduce overhead
function contract!(z, x::S, y::T, cix::Int=2, ciy::Int=1, conjx::Bool=false, conjy::Bool=false) where {S<:AbstractArray{<:Number, 2}, T<:AbstractArray{<:Number, 2}}
    _contract_checkdims!(x, y, cix, ciy)
    if (ndims(z) != 2) || (size(z, 1) != size(x, cix == 1 ? 2 : 1)) || (size(z, 2) != size(y, ciy == 1 ? 2 : 1))
        throw(ArgumentError("Destination tensor has the wrong dimensions."))
    end
    if conjx x = (cache(eltype(x), size(x), 1) .= conj.(x)) end
    if conjy y = cache(eltype(y), size(y), (!conjx || (length(x) != length(y))) ? 1 : 2) .= conj.(y) end
    mul!(z, cix == 2 ? x : transpose(x), ciy == 1 ? y : transpose(y))
end
function contract(x::S, y::T, cix::Int=2, ciy::Int=1, conjx::Bool=false, conjy::Bool=false) where {S<:AbstractArray{<:Number, 2}, T<:AbstractArray{<:Number, 2}}
    z = zeros(Base.promote_op(*, eltype(x), eltype(y)), size(x, cix == 2 ? 1 : 2), size(y, ciy == 2 ? 1 : 2))
    contract!(z, x, y, cix, ciy, conjx, conjy)
    return z
end


### Unsafe contractions
# Contraction with permutation 
function _contract!(z, x, y, sx, sy, cix, ciy, rix, riy, pix, piy, conjx, conjy)
    # Permute tensors
    px = permutedims!(cache(eltype(x), _contract_dims(sx, pix)), x, pix)
    py = permutedims!(cache(eltype(y), _contract_dims(sy, piy), length(x) == length(y) ? 2 : 1), y, piy)

    # Conjugations 
    if conjx px .= conj.(px) end
    if conjy py .= conj.(py) end

    # Contract
    mul!(reshape(z, (prod(_contract_dims(sx, rix)), prod(_contract_dims(sy, riy)))),
         reshape(px, (prod(_contract_dims(sx, rix)), prod(_contract_dims(sx, cix)))),
         reshape(py, (prod(_contract_dims(sy, ciy)), prod(_contract_dims(sy, riy)))))
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
function _contract_checkdims!(x, y, cix, ciy)
    if length(cix) != length(ciy)
        throw(ArgumentError("Contraction indices have different lengths"))
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
function _contract_checkreturn!(z, sx, sy, rix, riy)
    if (size(z) != (_contract_dims(sx, rix)..., _contract_dims(sy, riy)...))
        throw(ArgumentError("Destination tensor has the wrong dimensions."))
    end
end