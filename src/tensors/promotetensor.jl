#=
    Functions for promoting tensors.
=#

function promote_tensor(dims, args...)
    T = _promote_tensor_eltype(args...)
    return zeros(T, dims)
end