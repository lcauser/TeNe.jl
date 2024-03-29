#=
    Taking the trace of tensors 
=#

export trace, trace!


### Tracing 
"""
    trace!(z, x, cix::Int...; kwargs)

Compute the trace of x over dimensions cix, and store the result in z.
In place version of trace.

# Optional Keyword Arguments
    
    - 'conj::Bool=false': take the conjugate?
"""
function trace!(z, x, cix::Int...; conj::Bool=false)
    # Fetch the dimensions & do checks 
    sx, rix, pos = _trace_dimensions(x, cix)
    _trace_check_args(sx, cix)
    _trace_check_result(z, x, sx, rix)

    # Do the trace 
    _trace!(z, x, sx, cix, rix, pos, conj)
end

"""
    trace(x, cix::Int...; kwargs)

Compute the trace of x over dimensions cix.

# Optional Keyword Arguments
    
    - 'conj::Bool=false': take the conjugate?
    - 'tocache::Bool=false': store the result in the second level of the cache?
    - 'sublevel::Int=1': if stored in cache, at which sublevel?
"""
function trace(x, cix::Int...; conj::Bool=false, tocache::Bool=false, sublevel::Int=1)
    # Fetch the dimensions & do checks 
    sx, rix, pos = _trace_dimensions(x, cix)
    _trace_check_args(sx, cix)
    
    # Create the tensor to store the result 
    if tocache 
        z = cache(eltype(x), map(x->sx[x], rix), 2, sublevel; backend=typeof(get_backend(x)))
    else
        z = zeros(eltype(x), map(x->sx[x], rix)...)
    end

    # Do the trace 
    _trace!(z, x, sx, cix, rix, pos, conj)
    return z
end


### Unsafe trace
function _trace!(z, x, sx, cix, rix, pos, conjx)
    # Can probably be improved by only looping over contraction indices once...
    z .= zero(eltype(z))
    for idxs in Base.Iterators.product((Base.OneTo(sx[rix[i]]) for i in eachindex(rix))...)
        for k = 1:sx[cix[1]]
            full_idxs = (pos[i] == 0 ? k : idxs[pos[i]] for i in Base.OneTo(ndims(x)))
            z[idxs...] += conjx ? conj(x[full_idxs...]) : x[full_idxs...]
        end
    end
end


### Finding the traced and non-traced dimensions 
function _trace_dimensions(x, cix)
    sx = size(x)
    rix = tuple(setdiff(1:ndims(x), cix)...)
    pos = zeros(Int, ndims(x))
    j = 1
    for i in Base.OneTo(ndims(x))
        if i in cix 
            pos[i] = 0
        else
            pos[i] = j 
            j += 1
        end
    end
    return sx, rix, tuple(pos...)
end

### Checks 
function _trace_check_args(sx, cix)
    if !all(x->sx[x]==sx[cix[begin]], cix)
        throw(ArgumentError("Trace dimensions have different lengths."))
    end
end

function _trace_check_result(z, x, sx, rix)
    if typeof(get_backend(x)) != typeof(get_backend(z))
        throw(ArgumentError("Tensors are stored on different backends."))
    end
    if size(z) != sx[rix]
        throw(ArgumentError("Destination tensor has the wrong dimensions."))
    end
end