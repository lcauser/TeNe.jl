#=
    Contains functionality for applying StateOperators to StateVectors.
=#


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
    if !(length(O) == length(O1) == length(O2)) || !_so_so_product_checkdims(O, O1, O2)
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
    if !(length(O) == length(O1) == length(O2)) || !_so_so_product_checkdims(O, O1, O2)
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
