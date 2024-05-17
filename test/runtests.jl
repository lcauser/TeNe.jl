using TeNe
using LinearAlgebra
using Test

# Tensors 
@testset "Tensors" begin
    include("tensors/contract.jl")
    include("tensors/tensorproduct.jl")
    include("tensors/trace.jl")
    include("tensors/permutedim.jl")
    include("tensors/combinedims.jl")
    include("tensors/exp.jl")
    include("tensors/svd.jl")
end

# Lattice types
@testset "Lattice Types" begin 
    include("latticetypes/qubits.jl")
    include("latticetypes/oplist.jl")
end

# State vectors 
@testset "State vectors" begin 
    include("statetensors/stateoperator.jl")
    include("statetensors/applyso.jl")
    include("statetensors/inner.jl")
    include("statetensors/trace.jl")
end

# MPS 
@testset "MPS" begin 
    include("mps/gmps.jl")
    include("mps/applympo.jl")
    include("mps/inner.jl")
    include("mps/trace.jl")
    include("mps/creatempo.jl")
    include("mps/projmps.jl")
    include("mps/dmrg.jl")
end

# Circuits 
@testset "Circuits" begin 
    include("circuits/qubitgates.jl")
    include("circuits/gates.jl")
    include("circuits/circuitlayers.jl")
end