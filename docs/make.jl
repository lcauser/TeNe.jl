push!(LOAD_PATH,"../src/")
using Documenter, TeNe

makedocs(
    sitename="TeNe.jl",
    pages = [
        "Introduction" => "index.md",
        "Examples" => "examples.md",
        "Tensors" => "tensors.md",
        "State tensors" => "statetensors.md",
        "Matrix product states" => "mps.md"
    ]
)
deploydocs(
    repo = "github.com/lcauser/TeNe.jl.git",
)