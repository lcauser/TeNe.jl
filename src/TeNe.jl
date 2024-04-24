module TeNe

# Dependancies
using LinearAlgebra 
using LRUCache
using StaticArrays
using KernelAbstractions
using HDF5

# Caching; intermediate contractions memory can be pre-allocated and reused.
const _CACHE_MEM_LIM = 4294967296
const _CACHE = LRU{Tuple{DataType, Int64, Int64, Int64, Int64}, Any}(maxsize=_CACHE_MEM_LIM, by=Base.summarysize)
include("cache.jl")

# Default settings 
const _TeNe_cutoff = 1e-16

# Tensors
include("tensors/contract.jl")
include("tensors/tensorproduct.jl")
include("tensors/trace.jl")
include("tensors/permutedim.jl")
include("tensors/combinedims.jl")
include("tensors/exp.jl")
include("tensors/svd.jl")

# State vectors 
include("statetensors/structures/abstractstatetensor.jl")
include("statetensors/structures/gstatetensor.jl")
include("statetensors/structures/statevector.jl")
include("statetensors/structures/stateoperator.jl")

# MPS 
include("mps/structures/abstractmps.jl")
include("mps/structures/gmps.jl")
include("mps/structures/mps.jl")
include("mps/structures/mpo.jl")
end
