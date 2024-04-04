#=
    Functions for combining dimensions
=#

export combinedims, combinedims!
export uncombinedims, uncombinedims!


### Functions for combining dimensions
"""
    combinedims(y, x, cixs; kwargs...)

Combine the dimensions cixs in tensor x.

# Key arguments

    - `return_copy=false`: Return the result in newly allocated memory? Only
       necessary if the combined dimensions are the last dimensions of x.
"""
function combinedims(x, cixs; return_copy=false)
    # Checks and permutations 
    _combinedims_check_dims(x, cixs)
    pixs, sr, sc = _combinedims_permuted_dims(x, cixs)
    shape = (sr..., prod(sc))

    # Combine the dimensions
    if !_combinedims_check_perm(pixs)
        y = reshape(x, shape)
        return (return_copy ? copy(y) : y, (cixs, sc))
    else
        y = reshape(similar(x), shape)
        _combinedims!(y, x, pixs, sr, sc)
        return (y, (cixs, sc))
    end
end

"""
    combinedims!(y, x, cixs)

Combine the dimensions cixs in tensor x, and store the result in y.
"""
function combinedims!(y, x, cixs)
    # Checks and permutations 
    _combinedims_check_dims(x, cixs)
    pixs, sr, sc = _combinedims_permuted_dims(x, cixs)
    _combinedims_check_result(y, sr, sc)

    # Combine the indices
    _combinedims!(y, x, pixs, sr, sc)
    return (cixs, sc)
end


### Unsafe combination
function _combinedims!(y, x, pixs, sr, sc)
    # Permute the dimensions
    y = reshape(y, (sr..., sc...))
    permutedims!(y, x, pixs)
end


### Find the permuted dimensions 
function _combinedims_permuted_dims(x, cixs)
    rixs = tuple(setdiff(Base.OneTo(ndims(x)), cixs)...)
    pixs = (rixs..., cixs...)
    sr = map(i -> size(x, i), rixs)
    sc = map(i -> size(x, i), cixs)
    return pixs, sr, sc
end


### Checks on the indices
function _combinedims_check_dims(x, cixs)
    if !all(y -> y <= ndims(x), cixs)
        throw(ArgumentError("The specified dimensions $(cixs) are out of range for the rank-$(ndims(x)) tensor x."))
    end
end

function _combinedims_check_result(y, sr, sc)
    if size(y) != (sr..., prod(sc))
        throw(ArgumentError("Destination tensor has the wrong dimensions."))
    end
end

function _combinedims_check_perm(pixs)
    # Check to see if permutation is necessary 
    if all(i -> i == pixs[i], Base.OneTo(length(pixs)))
        return false
    else
        return true 
    end
end

### Functions for uncombining dimensions
function uncombinedims!(y, x, cixs, scs)
    _uncombinedims_check_dims(x, cixs, scs)
    sxs, perms = _uncombinedims_perms(x, cixs, scs)
    _uncombinedims_check_result(y, sxs, perms)
    _uncombinedims(y, x, sxs, perms)
end


### Unsafe uncombination of dimensions 
function _uncombinedims!(y, x, sxs, perms)
    x = reshape(x, sxs)
    permutedims!(y, x, perms)
end



### Finding the permutations to uncombine 
mutable struct _uncombinedims_counter
    j::Int64
end
function _uncombinedims_perms(x, cixs, scs)
    # Reshape sizes 
    sxs = (map(i -> size(x, i), Base.OneTo(ndims(x)-1))..., scs...)

    # Dimensions to permute 
    j = _uncombinedims_counter(1)
    function f(i, cixs, j, k)
        if i in cixs
            return k + findfirst(i .== cixs)
        else
            j.j += 1
            return j.j - 1
        end
    end
    perms = map(i -> f(i, cixs, j, ndims(x)-length(scs)+1), Base.OneTo(ndims(x)-1+length(cixs)))

    return sxs, perms 
end

### Checks 
function _uncombinedims_check_dims(x, cixs, scs)
    if prod(scs) != size(x, ndims(x))
        throw(ArgumentError("Incorrect sizes for combined dimensions."))
    end
    max_len = ndims(x) - 1 + len(cixs)
    if !all(i -> (i <= max_len && i >= 1), cixs)
        throw(ArgumentError("Incorrect permutation of dimensions."))
    end
end

function _uncombineidxs_check_result(y, sxs, perms)
    if ndims(y) != length(sxs)
        throw(ArgumentError("Desination tensor has incorrect number of dimensions."))
    end
    if !all(i -> sxs[perms[i]] == size(y, i), Base.OneTo(length(sxs)))
        throw(ArgumentError("Desination tensor has wrong dimensions: $(size(y)) != $(sxs[perms])"))
    end
end