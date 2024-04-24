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
export StateOperator

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

# Examples 

```julia-repl
julia> O = randomso(2, 10);
```
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

# Examples 

```julia-repl
julia> O = productso(10, [1 0; 0 1]);
```
"""
function productso(N::Int, A::AbstractMatrix; kwargs...)
    return productgst(N, A; kwargs...)
end
productstateoperator(N::Int, A::AbstractMatrix; kwargs...) = productso(N, A; kwargs...)


### Applying a StateOperator to a StateVector 
export applyso, applyso!
"""
    applyso(O::StateOperator, ψ::StateVector)
    applyso(ψ::StateVector, O::StateOperator)
    *(O::StateOperator, ψ::StateVector)
    *(ψ::StateVector, O::StateOperator)

Multiply the StateOperator `O` to the StateVector `ψ`.

# Examples 

```julia-repl
julia> O = productso(10, [0 1; 1 0])
julia> ψ = productsv(10, [1, 0])
julia> ϕ = O * ψ;
```
"""
function applyso(O::StateOperator, ψ::StateVector)
    ϕ = StateVector(dim(ψ), length(ψ); T=_promote_tensor_eltype(O, ψ))
    applyso!(ϕ, O, ψ)
    return ϕ
end
applyso(ψ::StateVector, O::StateOperator) = applyso(transpose(O), ψ)
*(O::StateOperator, ψ::StateVector) = applyso(O, ψ)
*(ψ::StateVector, O::StateOperator) = applyso(transpose(O), ψ)

"""
    applyso!(ϕ::StateVector, O::StateOperator, ψ::StateVector)
    applyso!(ϕ::StateVector, ψ::StateVector, O::StateOperator)

Multiply the StateOperator `O` to the StateVector `ψ`, and store the result in `ϕ`.

# Examples 

```julia-repl
julia> O = productso(10, [0 1; 1 0])
julia> ψ = productsv(10, [1, 0])
julia> ϕ = StateVector(2, 10)
julia> applyso!(ϕ, O, ψ);
```
"""
function applyso!(ϕ::StateVector, O::StateOperator, ψ::StateVector)
    if !issimilar(O, ψ, ϕ)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    iO = !istranspose(O) ? collect(2:2:ndims(tensor(O))) : collect(1:2:ndims(tensor(O))) # Remove allocations?...
    iψ = collect(1:ndims(tensor(ψ))) # Remove allocations?...
    conjO = isconj(ϕ) ? !isconj(O) : isconj(O)
    conjψ = isconj(ϕ) ? !isconj(ψ) : isconj(ψ) 
    contract!(tensor(ψ), tensor(O), tensor(ψ), iO, iψ, conjO, conjψ)
end
applyso!(ϕ::StateVector, ψ::StateVector, O::StateOperator) = applyso!(ϕ, transpose(O), ψ)

"""
    applyso!(O::StateOperator, ψ::StateVector)
    applyso!(ψ::StateVector, O::StateOperator)

Multiply the StateOperator `O` to the StateVector `ψ`. In-place version of `applyso(O, ψ)`.

# Examples 

```julia-repl
julia> O = productso(10, [0 1; 1 0])
julia> ψ = productsv(10, [1, 0])
julia> applyso!(O, ψ);
```
"""
function applyso!(O::StateOperator, ψ::StateVector)
    if !issimilar(O, ψ)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    iO = !istranspose(O) ? collect(2:2:ndims(tensor(O))) : collect(1:2:ndims(tensor(O))) # Remove allocations?...
    iψ = collect(1:ndims(tensor(ψ))) # Remove allocations?...
    conjO = isconj(ψ) ? !isconj(O) : isconj(O)
    ten = contract(tensor(O), tensor(ψ), iO, iψ, conjO, false)
    tensor(ψ) .= ten
end

### Applying a StateOperator to a StateOperator

"""
    applyso(O1::StateOperator, O2::StateOperator)

Calculate the product of two StateOperators, `O1` and `O2`.

# Examples 

```julia-repl
julia> O1 = productso(10, [0 1; 1 0])
julia> O2 = productso(10, [0 1im; -1im 0])
julia> O = O1 * O2;
```
"""
function applyso(O1::StateOperator, O2::StateOperator)
    O = StateOperator(dim(O1), length(O1); T=_promote_tensor_eltype(O1, O2))
    applyso!(O, O1, O2)
    return O
end
*(O1::StateOperator, O2::StateOperator) = applyso(O1, O2)

"""
    applyso!(O::StateOperator, O1::StateOperator, O2::StateOperator)

Calculate the product of two StateOperators, `O1` and `O2`. Store the result 
in StateOperator `O`.

# Examples 

```julia-repl
julia> O1 = productso(10, [0 1; 1 0])
julia> O2 = productso(10, [0 1im; -1im 0])
julia> O = StateOperator(2, 10)
julia> applyso!(O, O1, O2);
```
"""
function applyso!(O::StateOperator, O1::StateOperator, O2::StateOperator)
    # Checks 
    if !issimilar(O, O1, O2)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    conjO1 = isconj(O) ? !isconj(O1) : isconj(O1)
    conjO2 = isconj(O) ? !isconj(O2) : isconj(O2)
    if istranspose(O)
        ten = contract(tensor(O2), tensor(O1),
            collect((istranspose(O2) ? 1 : 2):2:ndims(tensor(O2))),
            collect((istranspose(O1) ? 2 : 1):2:ndims(tensor(O1))), conjO2, conjO1)
    else
        ten = contract(tensor(O1), tensor(O2),
            collect((istranspose(O1) ? 1 : 2):2:ndims(tensor(O1))),
            collect((istranspose(O2) ? 2 : 1):2:ndims(tensor(O2))), conjO1, conjO2)
    end
    permutedims!(tensor(O), ten, _applyso_perm_dims(length(O)))
end


function _applyso_perm_dims(N::Int)
    ### Remove allocations??
    return map(j->isodd(j) ? cld(j, 2) : N + cld(j, 2), Base.OneTo(2*N))
end