#=
    SimpleTensors uses an LRU cache to store temporary allocations to reduce
    allocations of memory. The cache stores data in a one-dimensional format,
    and can be accessed by specifying the eltype and the length of the array.
    Additionally, there is support for different levels, such that different
    arrays do not point to the same memory. Also allows for different entries
    for different threads, so that it is multi-threading safe.
=#

# Returning & allocation memory from / to the cache!
function cache(T::DataType, length::Int, level::Int, threadid::Int)
    if !haskey(_cache, (T, length, level, threadid))
        _cache[(T, length, level, threadid)] = zeros(T, length)
    end
    return _cache[(T, length, level, threadid)]
end

function cache(T::DataType, dims, level::Int, threadid::Int)
    return reshape(cache(T, prod(dims), level, threadid), dims)
end

cache(T::DataType, length, level) = cache(T, length, level, Threads.threadid())
cache(T::DataType, length) = cache(T, length, 1, Threads.threadid())

# Delete from the cache to prompt GC to remove allocations
function deletecache!(T::DataType, length::Int, level::Int, threadid::Int)
    _cache[(T, length, level, threadid)] = nothing
end
deletecache!(T::DataType, length::Int, level::Int) = deletecache!(T, length, level, Threads.threadid())
deletecache!(T::DataType, length::Int) = deletecache!(T, length, 1, Threads.threadid())