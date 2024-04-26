#=
    Operator lists are a high-level interface for strings of operators,
    such as Hamiltonians.
=#

mutable struct OpList{Q<:Number}
    lt::LatticeTypes
    length::Int
    ops::Vector{Vector{String}}
    sites::Vector{Vector{Int}}
    coeffs::Vector{Q}
end


### Initiate
export OpList
function OpList(lt::LatticeTypes, length::Int)
    return OpList{eltype(lt)}(lt, length, Vector{String}[], Vector{Int}[], eltype(lt)[])
end