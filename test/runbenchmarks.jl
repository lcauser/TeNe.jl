using BenchmarkTools
using TeNe 
include("benchmarks/dmrg.jl")
include("benchmarks/stateoptimiser.jl")

suite = BenchmarkGroup()
suite["dmrg"] = @benchmarkable benchmark_dmrg()
suite["state_optimiser"] = @benchmarkable benchmark_dmrg()

tune!(suite)
results = run(suite, verbose = true)

BenchmarkTools.save("output.json", median(results))