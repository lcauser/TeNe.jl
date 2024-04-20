push!(LOAD_PATH,"../src/")
using Documenter, TeNe

makedocs(sitename="TeNe.jl")
deploydocs(
    repo = "github.com/lcauser/TeNe.jl.git",
)