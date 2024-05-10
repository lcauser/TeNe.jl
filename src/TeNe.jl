module TeNe

# Dependancies
using LinearAlgebra 
using LRUCache
using StaticArrays
using KernelAbstractions
using HDF5
using Printf
using KrylovKit

# Imports 
import Base: +, -, *, /

# Caching; intermediate contractions memory can be pre-allocated and reused.
const _CACHE_MEM_LIM = 4294967296
const _CACHE = LRU{Tuple{DataType, Int64, Int64, Int64, Int64}, Any}(maxsize=_CACHE_MEM_LIM, by=Base.summarysize)
include("cache.jl")
#KernelAbstractions.get_backend(::SVector) = CPU
#KernelAbstractions.get_backend(::SMatrix) = CPU

# Default settings 
const _TeNe_cutoff = 1e-12

abstract type TensorNetworkState end
function issimilar(ψs::TensorNetworkState...)
    for i = Base.range(2, length(ψs))
        length(ψs[i]) != length(ψs[1]) && return false 
        dim(ψs[1]) != dim(ψs[i]) && return false
    end
    return true
end
export issimilar

### Tensors
include("tensors/promotetensor.jl")
include("tensors/contract.jl")
include("tensors/tensorproduct.jl")
include("tensors/trace.jl")
include("tensors/permutedim.jl")
include("tensors/combinedims.jl")
include("tensors/exp.jl")
include("tensors/svd.jl")

### LatticeTypes 
include("latticetypes/latticetype.jl")
include("latticetypes/qubits.jl")
include("latticetypes/oplist.jl")

### Tensor network states
# State Tensors
include("statetensors/structures/abstractstatetensor.jl")
include("statetensors/structures/gstatetensor.jl")
include("statetensors/structures/statevector.jl")
include("statetensors/structures/stateoperator.jl")

# MPS 
include("mps/structures/abstractmps.jl")
include("mps/structures/gmps.jl")
include("mps/structures/mps.jl")
include("mps/structures/mpo.jl")
include("mps/structures/mpsprojector.jl")

### Circuits 
include("circuits/gates.jl")
include("circuits/qubitgates.jl")
include("circuits/circuitlayer.jl")

### Type validation 
# Type abstraction 
const TensorNetworkVector = Union{StateVector, MPS}
const TensorNetworkOperator = Union{StateOperator, MPO, MPSProjector, OpList}

# Tensor networks 
include("validation/vec_vec.jl")
include("validation/op_vec.jl")
include("validation/vec_op_vec.jl")
include("validation/op_op.jl")
include("validation/op_trace.jl")

# Circuits
include("validation/gate_vec.jl")
include("validation/gate_op.jl")

### Tensor Network Projections 
# MPS 
include("mps/projections/abstractmpsprojection.jl")
include("mps/projections/projmps.jl")
include("mps/projections/projmpssquared.jl")
#include("mps/projections/projmpsprojector.jl") # Disabled for now due to unknown bug

### Tensor Network Operations 
# State Tensors 
include("statetensors/operations/applyso.jl")
include("statetensors/operations/inner.jl")
include("statetensors/operations/trace.jl")

# MPS
include("mps/operations/applympo.jl")
include("mps/operations/applyprojector.jl")
include("mps/operations/inner.jl")
include("mps/operations/trace.jl")
include("mps/operations/creatempo.jl")

### Optmisation methods 
# MPS 
include("mps/optimisers/mpsoptimiser.jl")
include("mps/optimisers/dmrg.jl")
end
