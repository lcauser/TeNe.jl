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
    _op_vec_validation(O, ψ)
    ϕ = GStateTensor(1, promote_tensor(innerdims(O), tensor(O), tensor(ψ)))
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
    _op_vec_validation(ϕ, O, ψ)
    _so_sv_product!(ϕ, O, ψ)
end
applyso!(ϕ::StateVector, ψ::StateVector, O::StateOperator) = applyso!(ϕ, transpose(O), ψ)


function _so_sv_product!(ϕ::StateVector, O::StateOperator, ψ::StateVector)
    contract!(tensor(ϕ), tensor(O), tensor(ψ),
        outerinds(O), Tuple(1:ndims(tensor(ψ))),
        isconj(ϕ) ? !isconj(O) : isconj(O), isconj(ϕ) ? !isconj(ψ) : isconj(ψ) )
end

function _so_sv_product_dims(O::StateOperator)
    return Tuple(map(j->innerdim(O, j), Base.OneTo(length(O))))
end

function _so_sv_product_check(ϕ::StateVector, O::StateOperator, ψ::StateVector)
    if length(ψ) != length(O) != length(ϕ)
        return false
    end
    for i in Base.OneTo(length(ψ))
        if dim(ψ, i) != outerdim(O, i) || dim(ϕ, i) != innerdim(O, i)
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
    _op_op_validation(O1, O2)
    O = GStateTensor(2, promote_tensor(_so_so_product_dims(O1, O2), O1, O2))
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
    _op_op_validation(O, O1, O2)
    _so_so_product!(O, O1, O2)
end

function _so_so_product_dims(O1::StateOperator, O2::StateOperator)
    ### Remove allocations??
    return Tuple(map(j->isodd(j) ? innerdim(O1, cld(j, 2)) : outerdim(O2, fld(j, 2)), Base.OneTo(2*length(O1))))
end

function _so_so_product_perms(N::Int)
    ### Remove allocations??
    return Tuple(map(j->isodd(j) ? cld(j, 2) : N + cld(j, 2), Base.OneTo(2*N)))
end

function _so_so_product!(O::StateOperator, O1::StateOperator, O2::StateOperator)
    conjO1 = isconj(O) ? !isconj(O1) : isconj(O1)
    conjO2 = isconj(O) ? !isconj(O2) : isconj(O2)
    if istranspose(O)
        ten = contract(tensor(O2), tensor(O1), outerinds(O2), innerinds(O1), conjO2, conjO1)
    else
        ten = contract(tensor(O1), tensor(O2), outerinds(O1), innerinds(O2), conjO1, conjO2)
    end
    permutedims!(tensor(O), ten, _so_so_product_perms(length(O)))
end
