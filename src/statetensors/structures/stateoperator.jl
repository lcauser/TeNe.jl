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
    return typeof(O) <: StateOperator
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
julia> O = productso(10, [0 1; 1 0]);
julia> ψ = productsv(10, [1, 0]);
julia> ϕ = O * ψ;
```
"""
function applyso(O::StateOperator, ψ::StateVector)
    ϕ = GStateTensor(1, dim(ψ), promote_tensor(_so_sv_product_dims(O), tensor(O), tensor(ψ)))
    if !(length(O) == length(ψ)) || !_so_sv_product_checkdims(ϕ, O, ψ)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    _so_sv_product!(ϕ, O ,ψ)
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
julia> O = productso(10, [0 1; 1 0]);
julia> ψ = productsv(10, [1, 0]);
julia> ϕ = StateVector(2, 10);
julia> applyso!(ϕ, O, ψ);
```
"""
function applyso!(ϕ::StateVector, O::StateOperator, ψ::StateVector)
    if !(length(ϕ) == length(O) == length(ψ)) || !_so_sv_product_checkdims(ϕ, O, ψ)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    _so_sv_product!(ϕ, O, ψ)
end
applyso!(ϕ::StateVector, ψ::StateVector, O::StateOperator) = applyso!(ϕ, transpose(O), ψ)


function _so_sv_product!(ϕ::StateVector, O::StateOperator, ψ::StateVector)
    # Properties
    iO = !istranspose(O) ? collect(2:2:ndims(tensor(O))) : collect(1:2:ndims(tensor(O))) # Remove allocations?...
    iψ = collect(1:ndims(tensor(ψ))) # Remove allocations?...
    conjO = isconj(ϕ) ? !isconj(O) : isconj(O)
    conjψ = isconj(ϕ) ? !isconj(ψ) : isconj(ψ) 

    # Contraction
    contract!(tensor(ϕ), tensor(O), tensor(ψ), iO, iψ, conjO, conjψ)
end

function _so_sv_product_dims(O::StateOperator)
    if !istranspose(O)
        return Tuple(map(j->size(tensor(O), 2*j-1), Base.OneTo(length(O))))
    else
        return Tuple(map(j->size(tensor(O), 2*j), Base.OneTo(length(O))))
    end
end

function _so_sv_product_checkdims(ϕ::StateVector, O::StateOperator, ψ::StateVector)
    for i in Base.OneTo(length(ψ))
        if size(tensor(ψ), i) != size(tensor(O), istranspose(O) ? 2*i : 2*i-1)
            return false
        end
        if size(tensor(ϕ), i) != size(tensor(O), istranspose(O) ? 2*i-1 : 2*i)
            return false
        end
    end
    return true
end

### Applying a StateOperator to a StateOperator
"""
    applyso(O1::StateOperator, O2::StateOperator)

Calculate the product of two StateOperators, `O1` and `O2`.

# Examples 

```julia-repl
julia> O1 = productso(10, [0 1; 1 0]);
julia> O2 = productso(10, [0 1im; -1im 0]);
julia> O = O1 * O2;
```
"""
function applyso(O1::StateOperator, O2::StateOperator)
    O = GStateTensor(2, dim(O1), promote_tensor(_so_so_product_dims(O1, O2), O1, O2))
    if !(length(O) == length(O1) == length(O2)) && _so_so_product_checkdims(O, O1, O2)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    _so_so_product!(O, O1, O2)
    return O
end
*(O1::StateOperator, O2::StateOperator) = applyso(O1, O2)

"""
    applyso!(O::StateOperator, O1::StateOperator, O2::StateOperator)

Calculate the product of two StateOperators, `O1` and `O2`. Store the result 
in StateOperator `O`.

# Examples 

```julia-repl
julia> O1 = productso(10, [0 1; 1 0]);
julia> O2 = productso(10, [0 1im; -1im 0]);
julia> O = StateOperator(2, 10);
julia> applyso!(O, O1, O2);
```
"""
function applyso!(O::StateOperator, O1::StateOperator, O2::StateOperator)
    if !(length(O) == length(O1) == length(O2)) && _so_so_product_checkdims(O, O1, O2)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    _so_so_product!(O, O1, O2)
end

function _so_so_product_checkdims(O::StateOperator, O1::StateOperator, O2::StateOperator)
    for i = Base.OneTo(length(O))
        if size(tensor(O), istranspose(O) ? 2*i : 2*i-1) != size(tensor(O1), istranspose(O1) ? 2*i : 2*i-1)
            return false
        end
        if size(tensor(O), istranspose(O) ? 2*i-1 : 2*i) != size(tensor(O2), istranspose(O2) ? 2*i-1 : 2*i)
            return false
        end
        if size(tensor(O1), istranspose(O1) ? 2*i-1 : 2*i) != size(tensor(O2), istranspose(O2) ? 2*i : 2*i-1)
            return false
        end
    end
    return true
end

function _so_so_product_dims(O1::StateOperator, O2::StateOperator)
    dims = map(j->isodd(j) ? size(tensor(O1), istranspose(O1) ? j+1 : j) :
        size(tensor(O2), istranspose(O2) ? j-1 : j), Base.OneTo(2*length(O1)))
    return Tuple(dims)
end

function _so_so_product_perms(N::Int)
    ### Remove allocations??
    return Tuple(map(j->isodd(j) ? cld(j, 2) : N + cld(j, 2), Base.OneTo(2*N)))
end

function _so_so_product!(O::StateOperator, O1::StateOperator, O2::StateOperator)
    conjO1 = isconj(O) ? !isconj(O1) : isconj(O1)
    conjO2 = isconj(O) ? !isconj(O2) : isconj(O2)
    if istranspose(O)
        ten = contract(tensor(O2), tensor(O1),
            Tuple((istranspose(O2) ? 1 : 2):2:ndims(tensor(O2))),
            Tuple((istranspose(O1) ? 2 : 1):2:ndims(tensor(O1))), conjO2, conjO1)
    else
        ten = contract(tensor(O1), tensor(O2),
            Tuple((istranspose(O1) ? 1 : 2):2:ndims(tensor(O1))),
            Tuple((istranspose(O2) ? 2 : 1):2:ndims(tensor(O2))), conjO1, conjO2)
    end
    permutedims!(tensor(O), ten, _so_so_product_perms(length(O)))
end


### Inner products 
export inner
"""
    inner(ψ::StateVector, O::StateOperator, ϕ::StateVector)
    inner(ψ::StateVector, O1::StateOperator, O2::MPO, ϕ::StateVector)

Calculate the expectation of a string of StateOperators `Os` with respect to StateVectors `ψ` and `ϕ`.

# Examples 

```jldoctest
julia> ψ = randomsv(2, 10);
julia> O = productso(10, [0 1; 1 0]);
julia> inner(ψ, O, O, ψ)
1.0 + 0.0im
```
"""
function inner(ψ::StateVector, ϕs::Union{StateVector, StateOperator}...)
    # Checks 
    for i = Base.OneTo(length(ϕs)-1)
        if !isstateoperator(ϕs[i])
            throw(ArgumentError("The inner terms in the braket must be MPOs."))
        end
    end
    if !isstatevector(ϕs[end])
        throw(ArgumentError("The last term in the MPS must be an MPO."))
    end

    # Create a StateVector from the cache 
    ϕ = GStateTensor(1, dim(ψ), cache(size(tensor(ϕs[end])), tensor(ϕs[end])))

    # Contraction
    applyso!(ϕ, ϕs[end-1], ϕs[end])
    for i = Base.range(length(ϕs)-2, 1, step=-1)
        ϕnew = GStateTensor(1, dim(ψ), cache(size(tensor(ϕs[end])), tensor(ϕs[end])))
        applyso!(ϕnew, ϕs[i], ϕ)
        ϕ = ϕnew
    end
    return inner(ψ, ϕ)
end


### Trace of StateOperators 
export trace
"""
    trace(Os::StateOperator...)

Compute the trace of a string of StateOperators.

# Examples 

```jldoctest
julia> O = productso(10, [0 1; 1 0]);
julia> trace(O, O)
1024.0 + 0.0im
```
"""
function trace(Os::StateOperator...)
    # Checks 
    if !issimilar(Os...)
        throw(ArgumentError("Arguments have properties that do not match."))
    end

    if length(Os) == 1
        ten = tensor(Os[1])
        for _ = 1:length(Os[1])
            ten = trace(ten, 1, 2)
        end
        return ten[]
    end

    # Create a tensor from the cache
    T = _promote_tensor_eltype(Os...) 
    perms = tuple(map(j->isodd(j) ? j + 1 : j - 1, Base.OneTo(ndims(tensor(Os[1]))))...)
    if istranspose(Os[1])
        dims = map(j->size(tensor(Os[1]), j), perms)
        ten = cache(T, dims, 2, 1)
        permutedims!(ten, tensor(Os[1]), perms)
    else
        dims = size(tensor(Os[1]))
        ten = cache(T, dims, 2, 1) .= tensor(Os[1])
    end
    if isconj(Os[1])
        ten .= conj.(ten)
    end
    
    # Contract with center operators 
    perms2 = _so_so_product_perms(length(Os[1]))
    for i in range(2, length(Os)-1)
        ten = contract(ten, tensor(Os[i]), collect(2:2:ndims(ten)),
            collect((istranspose(Os[i]) ? 2 : 1):2:ndims(tensor(Os[i]))),
            false, isconj(Os[i]))
        ten = permutedims(ten, perms2)
    end

    # Contract with the last 
    dims = Tuple(Base.OneTo(ndims(ten)))
    if istranspose(Os[end])
        return contract(ten, tensor(Os[end]), dims, dims, false, isconj(Os[end]))[]
    else
        return contract(ten, tensor(Os[end]), dims, perms, false, isconj(Os[end]))[]
    end
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