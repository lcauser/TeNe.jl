push!(LOAD_PATH,"../src/")
using TeNe
using Documenter
makedocs(
         sitename = "TeNe.jl",
         modules  = [TeNe],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/lcauser/TeNe.jl",
)