#=
    QR decomposition of tensors.
=#


### Unsafe QR 
function _qr(x, dims; enforce_positive = false)
    # Permute and reshape 
    y, sinners, souters = _svd_reshape(x, dims)

    # Do the QR using LinearAlgebra 
    Q, R = LinearAlgebra.qr(y)

    # Make the diagonal elements of R positive 
    if enforce_positive
        M = diagm([real(R[i, i]) < 0.0 ? -1 : 1 for i in eachindex(R[1, :])])
        Q = contract(Q, M, 2, 1; tocache = false)
        R = contract(M, R, 2, 1; tocache = false)
    end

    # Ungroup the dimensions
    U = reshape(U, (sinners..., vals))
    V = reshape(V, (vals, souters...))
end


### Permuting dimensions 
function _qr_permute_dims(x, dims)
    inner_dims = tuple(setdiff(1:ndims(x), dims)...)
    sinners = map(i->size(x, i), inner_dims)
    souters = map(i->size(x, i), dims)
    return sinners, souters, tuple(inner_dims..., dims...)
end

### Checks 
function _qr_checkdims(x, dims)
    if ndims(x) < length(dims)
        throw(ArgumentError("The number of dimensions given exceed that of the tensor."))
    end
    if !all(i -> isinteger(i), dims)
        throw(ArgumentError("Dimensions $(dims) are not integers."))
    end
    if !all(i -> (i > 0 && i <= ndims(x)), dims)
        throw(
            ArgumentError("Dimensions $(dims) exceed the size of the tensor: $(ndims(x))."),
        )
    end
end
