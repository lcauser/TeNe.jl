#=
    State operators are an exact representation of matrices for many-body systems with discrete
    degrees of freedom.
=#

### Traits for State Operators 

# Transpose
struct TransposeStateOperator <: GStateTensorTrait
    StateTensor::GStateTensor{2}
end
TeNe.transpose(O::GStateTensor{2}) = TransposeStateOperator(O)
TeNe.transpose(O::TransposeStateOperator) = O.StateTensor
istranspose(ψ::TransposeStateOperator) = true

# Adjoint
struct AdjointStateOperator <: GStateTensorTrait
    StateTensor::GStateTensor{2}
end
TeNe.adjoint(O::GStateTensor{2}) = AdjointStateOperator(O)
TeNe.adjoint(O::AdjointStateOperator) = O.StateTensor
isconj(ψ::AdjointStateOperator) = true
istranspose(ψ::AdjointStateOperator) = true

# Trait rules
TeNe.conj(O::TransposeStateOperator) = AdjointStateOperator(O.StateTensor)
TeNe.conj(O::AdjointStateOperator) = TransposeStateOperator(O.StateTensor)
TeNe.transpose(O::ConjGStateTensor{2}) = AdjointStateOperator(O.StateTensor)
TeNe.transpose(O::AdjointStateOperator) = ConjGStateTensor(O.StateTensor)
TeNe.adjoint(O::ConjGStateTensor{2}) = TransposeStateOperator(O.StateTensor)
TeNe.adjoint(O::TransposeStateOperator) = ConjGStateTensor(O.StateTensor)

const StateOperator = Union{GStateTensor{2}, ConjGStateTensor{2}, TransposeStateOperator, AdjointStateOperator}
export StateOpeartor

export isstateoperator 
"""
    isstateoperator(O)

Check to see if an object is a state operator.
"""
function isstateoperator(O)
    return typoeof(O) <: StateOperator
end

### Initialising StateOperators 
export randomso, randomstateoperator, productso, productstateoperator
"""
    StateOperator(dim::Int, length::Int)

Create a StateOperator with physical dimensions `dim` and `length` sites.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensor.
"""
function StateOperator(dim::Int, length::Int; kwargs...)
    return GStateTensor(2, dim, length; kwargs...)
end


"""
    randomso(dim::Int, length::Int; kwargs...)
    randomstateoperator(dim::Int, length::Int; kwargs...)

Create a state operator with physical dimensions `dim` and `length` sites.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensor.
"""
function randomso(dim::Int, length::Int; kwargs...)
    return randomgst(2, dim, length; kwargs...)
end
randomstateoperator(dim::Int, length::Int; kwargs...) = randomso(dim, length; kwargs...)

"""
    productso(N::Int, A::AbstractMatrix; kwargs...)
    productstateoperator(N::Int, A::AbstractMatrix; kwargs...)

Create a state operator as a product state with size `N`, and composed of array `A`.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensor.
"""
function productso(N::Int, A::AbstractMatrix; kwargs...)
    return productgst(N, A; kwargs...)
end
productstateoperator(N::Int, A::AbstractMatrix; kwargs...) = productso(N, A; kwargs...)


### Applying a StateOperator to a StateVector 
