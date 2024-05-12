#=
    A circuit contains a full list of unitary operations and measurements that 
    can be applied to a tensor network state (or operator [later]).
=#

mutable struct Circuit{d}
    N::Int
    layers::Vector{CircuitLayer}
    gates::Vector{<:AbstractGate}
    sites::Vector{Tuple{Vararg{Int}}}
    assigned::Vector{Bool} 
end