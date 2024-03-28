module TeNe

# Dependancies
using LinearAlgebra 
using LRUCache
using KernelAbstractions

# Create the cache 
const _CACHE_MEM_LIM = 4294967296
const _CACHE = LRU{Tuple{DataType, Int64, Int64, Int64, Int64}, Any}(maxsize=_CACHE_MEM_LIM, by=Base.summarysize)

# exports 
export contract, contract!

# Files 
include("cache.jl")
include("contract.jl")

end
