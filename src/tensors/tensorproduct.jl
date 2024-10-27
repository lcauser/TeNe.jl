#=
    Taking the tensor product of two tensors 
=#

export tensorproduct, tensorproduct!


### Tensor products
"""
    tensorproduct!(z, x, y, [conjx=false, conjy=false])

Compute the tensor product of the two tensors `x` and `y`, and store the 
result in `z`. Optionally, do the tensor product using the conjugate of
the tensors.

# Arguments

    - `z`: tensor to store the result.
    - `x`: first tensor.
    - `y': second tensor.
    - `conjx::Bool=false`: Take the complex conjugate of argument x?
    - `conjy::Bool=false`: Take the complex conjugate of argument y?

# Examples 

```jldoctest
julia> x = randn(ComplexF64, 2, 3);
julia> y = randn(ComplexF64, 4, 5);
julia> z = similar(x, (2, 3, 4, 5));
julia> tensorproduct!(z, x, y);
```
"""
function tensorproduct!(z, x, y, conjx::Bool=false, conjy::Bool=false)
    # Checks on args 
    _tensorproduct_check_result(z, x, y)
    _tensorproduct!(z, x, y, conjx, conjy)
end


"""
    tensorproduct(x, y, [conjx=false, conjy=false])

Compute the tensor product of the two tensors `x` and `y`, and store the 
result in `z`. Optionally, do the tensor product using the conjugate of
the tensors.

# Arguments

    - `x`: first tensor.
    - `y': second tensor.
    - `conjx::Bool=false`: Take the complex conjugate of argument x?
    - `conjy::Bool=false`: Take the complex conjugate of argument y?

# Optional Keyword Arguments
    
    - `tocache::Bool=true`: store the result in the second level of the cache?
    - `sublevel=:auto`: if stored in cache, at which sublevel? :auto finds non-aliased memory

# Examples 

```jldoctest
julia> x = randn(ComplexF64, 2, 3);
julia> y = randn(ComplexF64, 4, 5);
julia> z = tensorproduct(x, y);
julia> size(z)
(2, 3, 4, 5)
```
"""
function tensorproduct(x, y, conjx::Bool=false, conjy::Bool=false;
                       tocache::Bool=true, sublevel=:auto)
    # Checks on arguments 
    _tensorproduct_check_args(x, y)

    # Create the tensor to store the result 
    dims = (size(x)...,  size(y)...)
    if tocache
        z = cache(dims, x, y; level=2, sublevel=sublevel)
    else
        z = promote_tensor(dims, x, y)
    end
    _tensorproduct!(z, x, y, conjx, conjy)
    return z
end


### Unsafe tensor product
function _tensorproduct!(z, x, y, conjx::Bool, conjy::Bool)
    # Take the conjugates?
    if conjx x = cache(x, 1) .= conj.(x) end
    if conjy y = cache(y, 1, length(x) == length(y) ? 2 : 1) .= conj(y) end 

    # Tensor product 
    mul!(reshape(z, (prod(size(x)), prod(size(y)))),
         reshape(x, (prod(size(x)), 1)),
         reshape(y, (1, prod(size(y)))))
end

### Checks 
function _tensorproduct_check_args(x, y)
    #if typeof(get_backend(x)) != typeof(get_backend(y))
    #    throw(ArgumentError("Tensors are stored on different backends."))
    #end
end

function _tensorproduct_check_result(z, x, y)
    #if typeof(get_backend(x)) != typeof(get_backend(y) != typeof(get_backend(z)))
    #    throw(ArgumentError("Tensors are stored on different backends."))
    #end
    if size(z) != (size(x)..., size(y))
        throw(ArgumentError("Destination tensor has the wrong dimensions."))
    end
end