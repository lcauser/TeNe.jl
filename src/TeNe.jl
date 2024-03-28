module TeNe

# Dependancies
using LinearAlgebra 
using LRUCache

# Create the cache 
const _cache = LRU{Tuple{DataType, Int64, Int64, Int64}, Any}(maxsize=10^5)

# exports 
export contract, contract!

# Files 
include("cache.jl")
include("contract.jl")

end
