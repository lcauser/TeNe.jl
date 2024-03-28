#=
    TeNe.jl uses an LRU cache to store temporary allocations to reduce
    allocations of memory. The cache stores data in a one-dimensional format,
    and can be accessed by specifying the eltype and the length of the array.
    Additionally, there is support for different levels, such that different
    arrays do not point to the same memory. Also allows for different entries
    for different threads, so that it is multi-threading safe.
=#

# Returning & allocation memory from / to the cache given the data type and size(s)
function cache(T::DataType, length::Int, level::Int, sublevel::Int, threadid::Int; kwargs...)
    backend::DataType = get(kwargs, :backend, CPU)
    if !haskey(_CACHE, (T, length, level, sublevel, threadid))
        _CACHE[(T, length, level, sublevel, threadid)] = zeros(T, length)
    end
    if backend != CPU
        return _CACHE[(T, length, level, sublevel, threadid)] |> backend
    else
        _CACHE[(T, length, level, sublevel, threadid)]
    end
end

function cache(T::DataType, dims, level::Int, sublevel::Int, threadid::Int; kwargs...)
    return reshape(cache(T, prod(dims), level, sublevel, threadid; kwargs...), dims)
end

function cache(T::DataType, length, level::Int=1, sublevel::Int=1; kwargs...) 
    return cache(T, length, level, sublevel, Threads.threadid())
end

# Accessing from an array 
function cache(x::Q, level::Int, sublevel::Int, threadid::Int) where {Q<:AbstractArray}
    backend = typeof(get_backend(x))
    return cache(eltype(x), size(x), level, sublevel, threadid; backend=backend)
end

function cache(x::Q, level::Int=1, sublevel::Int=1) where {Q<:AbstractArray}
    return cache(x, level, sublevel, Threads.threadid())
end

function cache(x::Q, dims, level::Int, sublevel::Int, threadid::Int) where {Q<:AbstractArray}
    return reshape(cache(x, level, sublevel, threadid), dims)
end

function cache(x::Q, dims, level::Int=1, sublevel::Int=1) where {Q<:AbstractArray}
    return reshape(cache(x, level, sublevel, Threads.threadid()), dims)
end


# Delete from the cache to prompt GC to remove allocations
function deletecache!(T::DataType, length::Int, level::Int, sublevel::Int, threadid::Int)
    _CACHE[(T, length, level, sublevel, threadid)] = nothing
end

function deletecache!(T::DataType, length::Int, level::Int=1, sublevel::Int=1) 
    deletecache!(T, length, level, sublevel, Threads.threadid())
end