# Tensors
```@contents
Pages = ["tensors.md"]
Depth = 3
```

## Tensor operations

### Tensor contraction
Contract two tensors `x` and `y` over dimensions `cix` and `ciy` with a shared size.
This can be done in place with a pre-allocated tensor `z` using `contract!(z, x, y, cix, ciy)`, or otherwise `contract(x, y, cix, ciy)`.
The dimensions `cix` and `ciy` can be specified as integers or tuples of integers for contractions over multiple dimensions.
Optionally, the tensors used in the contraction can be conjugated using the arguments `conjx` and `conjy`.
Note that, by default, the result will be stored in the memory cache (using keyword argument `tocache`).
If the contraction is not some intermediate step, and you would like to save the resulting tensor for future use, then use `tocache=false`.
```@docs
contract!
contract
```

### Tensor product 
Take the tensor product over two tensors `x` and `y` to give a single tensor.
This can be done in place with a pre-allocated tensor `z` using `tensorproduct!(z, x, y)`, or otherwise `tensorproduct(x, y)`.
Optionally, the tensors used in the product can be conjugated using the arguments `conjx` and `conjy`.
Note that, by default, the result will be stored in the memory cache (using keyword argument `tocache`).
If the result is not some intermediate step, and you would like to save the resulting tensor for future use, then use `tocache=false`.
```@docs
tensorproduct!
tensorproduct
```

### Tensor trace 
Compute the trace over multiple dimensions `cix` in a tensor `x`.
This can be done in place with a pre-allocated tensor `z` using `trace!(z, x, cix...)`, or otherwise `trace(x, cix...)`.
Optionally, the tensor used in the trace can be conjugated using the key word argument `conj`.
Note that, by default, the result will be stored in the memory cache (using keyword argument `tocache`).
If the result is not some intermediate step, and you would like to save the resulting tensor for future use, then use `tocache=false`.
```@docs
trace!
trace(x, cix::Int...)
```

### Permuting a single dimension
Permute the dimension at position `i` in tensor `x` to position `j`.
This can be done in place with a pre-allocated tensor `z` using `permutedim!(z, x, i, j)`, or otherwise `permutedim(x, i, j)`.
Note that, by default, the result will be stored in the memory cache (using keyword argument `tocache`).
If the result is not some intermediate step, and you would like to save the resulting tensor for future use, then use `tocache=false`.
```@docs
permutedim!
permutedim
```

### Combining & restoring dimensions 
Dimensions in a tensor can be combined into a single dimension, and restored using a key.
This allows us to make efficient use of BLAS and LAPACK routines involving matrix operations. 

Combine the dimensions `cixs` of tensor `x`.
This can be done in place with a pre-allocated tensor `z` using `key = combinedims!(z, x, cixs)`, or otherwise `z, key = combinedims(x, cixs)`.
Note that, by default, the result will be stored in the memory cache (using keyword argument `tocache`).
If the result is not some intermediate step, and you would like to save the resulting tensor for future use, then use `tocache=false`.
```@docs
combinedims!
combinedims
```

After combining the dimensions, which returns a key, the dimensions of the tensor can be restored.
This can be done in place with a pre-allocated tensor `y` using `uncombinedims!(y, x, key)`, or otherwise `y = uncombinedims(y, x, key)`.
Note that, by default, the result will be stored in the memory cache (using keyword argument `tocache`).
If the result is not some intermediate step, and you would like to save the resulting tensor for future use, then use `tocache=false`.
```@docs
uncombinedims!
uncombinedims
```

## Tensor factorisations

### Singular value decomposition
The singular value decomposition ``M = USV`` is typically applied to a matrix, but can equally be applied to a tensor to split the dimensions into separate tensors.
This is done by permuting and reshaping the tensor into a matrix representation and applying the SVD.
The dimensions to be contained in ``V`` are specified by `dims`, with `U, S, V = tsvd(x, dims)`.

```@docs 
tsvd
```

At a later date, we would like to improve the SVD to pre-allocate memory in the cache for calculating the returns (and optional parameters to pre-allocated the memory to restore the results.)