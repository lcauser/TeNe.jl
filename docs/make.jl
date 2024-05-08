push!(LOAD_PATH,"../src/")
using Documenter, TeNe

makedocs(
    sitename="TeNe.jl",
    pages = [
        "Introduction" => "index.md",
        "Examples" => ["examples/dmrg.md"],
        "Manual" => ["manual/tensors.md", "manual/statetensors.md", "manual/mps.md"]
    ]
)
deploydocs(
    repo = "github.com/lcauser/TeNe.jl.git",
)