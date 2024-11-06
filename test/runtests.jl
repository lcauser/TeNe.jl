using TeNe
using LinearAlgebra
using Test

### Unit tests: 
@testset "Unit Tests" begin
    # Tensors 
    @testset "Tensors" begin
        include("unit/tensors/contract.jl")
        include("unit/tensors/tensorproduct.jl")
        include("unit/tensors/trace.jl")
        include("unit/tensors/permutedim.jl")
        include("unit/tensors/combinedims.jl")
        include("unit/tensors/exp.jl")
        include("unit/tensors/svd.jl")
    end

    # Lattice types
    @testset "Lattice Types" begin 
        include("unit/latticetypes/oplist.jl")
        include("unit/latticetypes/qubits.jl")
        include("unit/latticetypes/bosons.jl")
    end

    # State vectors 
    @testset "State vectors" begin 
        include("unit/statetensors/structures/stateoperator.jl")
        include("unit/statetensors/operations/applyso.jl")
        include("unit/statetensors/operations/inner.jl")
        include("unit/statetensors/operations/trace.jl")
    end

    # MPS 
    @testset "MPS" begin 
        include("unit/mps/structures/gmps.jl")
        include("unit/mps/operations/applympo.jl")
        include("unit/mps/operations/inner.jl")
        include("unit/mps/operations/trace.jl")
        include("unit/mps/operations/creatempo.jl")
        include("unit/mps/projections/projmps.jl")
    end

    # Circuits 
    @testset "Circuits" begin 
        include("unit/circuits/qubitgates.jl")
        include("unit/circuits/gates.jl")
        include("unit/circuits/circuitlayers.jl")
        include("unit/circuits/circuits.jl")
        include("unit/circuits/projections/projmpscircuit.jl")
        include("unit/circuits/operations/trotterization.jl")
    end
end

### End-to-end testing
@testset "End-to-end" begin 
    include("end-to-end/mps/dmrg.jl")
    include("end-to-end/mps/stateoptimiser.jl")
end