module TeNe

# Dependancies
using LinearAlgebra 
using LRUCache
using KernelAbstractions

# Caching; intermediate contractions memory can be pre-allocated and reused.
const _CACHE_MEM_LIM = 4294967296
const _CACHE = LRU{Tuple{DataType, Int64, Int64, Int64, Int64}, Any}(maxsize=_CACHE_MEM_LIM, by=Base.summarysize)
include("cache.jl")

# Tensors
include("tensors/contract.jl")
include("tensors/tensorproduct.jl")
include("tensors/trace.jl")
include("tensors/permutedim.jl")
include("tensors/combinedims.jl")

end
