module SimpleTensors

# Dependancies
using LinearAlgebra 
using LRUCache

# Create the cache 
const _cache = LRU{Tuple{DataType, Int64, Int64, Int64}, Any}(maxsize=10^5)

# Files 
include("cache.jl")
include("contract.jl")

end
