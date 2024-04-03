#=
    Functions for combining dimensions
=#

export combinedims, combinedims!


### Functions for combining dimensions


function combinedims(x, cixs; return_copy=false)
    # Checks and permutations 
    _combinedims_check_dims(x, cixs)
    pixs, sr, sc = _combinedims_permuted_dims(x, cixs)
    shape = (map(sr)..., prod(sc))

    # Combine the dimensions
    if _combinedims_check_perm(pixs)
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
    rixs = tuple(setdiff(1:ndims(x), cixs)...)
    pixs = (rixs..., cixs...)
    sr = map(i -> size(x, i), rixs)
    sc = map(i -> size(x, i), cixs)
    return pixs, sr, sc
end


### Checks on the indices
function _combinedims_check_dims(x, cixs)
    if !all(y -> y > ndims(x), cixs)
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
    if all(i -> i == pixs[i] for i = 1:length(pixs))
        return false
    else
        return true 
    end
end