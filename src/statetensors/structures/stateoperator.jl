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

const StateOperator = Union{GStateTensor{2}, ConjGStateTensor{2},
    TransposeStateOperator, AdjointStateOperator}
export StateOperator

export isstateoperator 
"""
    isstateoperator(O)

Check to see if an object is a state operator.
"""
function isstateoperator(O)
    return typeof(O) <: StateOperator
end

# Indices
innerind(::Union{GStateTensor{2}, ConjGStateTensor{2}}, site::Int) = 2*site-1
innerind(::Union{AdjointStateOperator, TransposeStateOperator}, site::Int) = 2*site
outerind(::Union{GStateTensor{2}, ConjGStateTensor{2}}, site::Int) = 2*site
outerind(::Union{AdjointStateOperator, TransposeStateOperator}, site::Int) = 2*site-1
innerinds(O::Union{GStateTensor{2}, ConjGStateTensor{2}}) = Tuple(1:2:2*length(O))
innerinds(O::Union{AdjointStateOperator, TransposeStateOperator}) = Tuple(2:2:2*length(O))
outerinds(O::Union{GStateTensor{2}, ConjGStateTensor{2}}) = Tuple(2:2:2*length(O))
outerinds(O::Union{AdjointStateOperator, TransposeStateOperator}) = Tuple(1:2:2*length(O))

# Physical dimensions 
export innerdim, outerdim, innerdims, outerdims
innerdim(O::Union{GStateTensor{2}, ConjGStateTensor{2}}, site::Int) = size(tensor(O), 2*site-1)
innerdim(O::Union{AdjointStateOperator, TransposeStateOperator}, site::Int) = size(tensor(O), 2*site)
outerdim(O::Union{GStateTensor{2}, ConjGStateTensor{2}}, site::Int) = size(tensor(O), 2*site)
outerdim(O::Union{AdjointStateOperator, TransposeStateOperator}, site::Int) = size(tensor(O), 2*site-1)
innerdims(O::Union{GStateTensor{2}, ConjGStateTensor{2}}) = Tuple(map(j->size(tensor(O), 2*j-1), Base.OneTo(length(O))))
innerdims(O::Union{AdjointStateOperator, TransposeStateOperator}) = Tuple(map(j->size(tensor(O), 2*j), Base.OneTo(length(O))))
outerdims(O::Union{GStateTensor{2}, ConjGStateTensor{2}}) = Tuple(map(j->size(tensor(O), 2*j), Base.OneTo(length(O))))
outerdims(O::Union{AdjointStateOperator, TransposeStateOperator}) = Tuple(map(j->size(tensor(O), 2*j-1), Base.OneTo(length(O))))

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


"""
    productso(lt::LatticeTypes, states::AbstractVector{String})
    productstateoperator(lt::LatticeTypes, states::AbstractVector{String})

Create a product operator from a string of operator names.

# Example 

```julia-repl
julia> lt = Qubits();
julia> ψ = productso(lt, ["x" for _ = 1:6]);
```
"""
function productso(lt::LatticeTypes, ops::AbstractVector{String})
    tensor = ones(eltype(lt), )
    for i in eachindex(ops)
        tensor = tensorproduct(tensor, op(lt, ops[i]); tocache = i!=lastindex(ops))
    end
    return GStateTensor(2, tensor)
end
productstateoperator(lt::LatticeTypes, states::AbstractVector{String}) = productso(lt, states)

### Initalising StateOperator from OpList 
"""
    StateOperator(ops::OpList)

Create a StateOperator from an OpList.

# Examples 

```julia-repl
julia> lt = Qubits();
julia> H = OpList(lt, 10);
julia> for i = 1:10 add!(H, "x", i) end;
julia> for i = 1:9 add!(H, ["z", "z"], [i, i+1]) end;
julia> H = StateOperator(H);
```
"""
function StateOperator(ops::OpList)
    O = StateOperator(dim(ops.lt), ops.length; T=eltype(ops.lt))
    for i in eachindex(ops.ops)
        ten = ones(eltype(ops.lt), )
        k = 1
        for j = Base.OneTo(length(O))
            if k <= length(ops.sites[i]) && ops.sites[i][k] == j
                oper = ops.ops[i][k]
                k += 1
            else
                oper = "id"
            end
            ten = tensorproduct(ten, op(ops.lt, oper))
        end
        tensor(O) .+= ten
    end
    return O
end


### Exponential of a state tensor 
export exp
"""
    exp(O::StateOperator; kwargs...)

Exponentiate a StateOperator.

# Optional Keyword Arguments

    - `prefactor=1.0`: Multiply the StateOperator by some value before the exponetiation.
"""
function TeNe.exp(O::StateOperator; prefactor=1.0)
    if istranspose(O)
        ten = cache(size(tensor(O)), tensor(O))
        dims = Tuple(map(j->isodd(j) ? j+1 : j-1, Base.OneTo(ndims(ten))))
        permutedims!(ten, tensor(O), dims)
    else
        ten = tensor(O)
    end
    if isconj(O)
        ten = cache(size(ten), ten) .= conj.(ten)
    end
    ten *= prefactor
    ten = TeNe.exp(ten, Base.range(2, 2*length(O), step=2))
    return GStateTensor(2, dim(O), ten)
end