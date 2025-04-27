#=
    Calculating inner products of StateVectors and StateOperators.
=#

export inner, dot

### Inner products of StateVectors
"""
    inner(ψ::StateVector, ϕ::StateVector)
    dot(ψ::StateVector, ϕ::StateVector)
    *(ψ::StateVector, ϕ::StateVector)

Calculate the inner product of two StateVectors `ψ` and `ϕ`.
"""
function inner(ψ::StateVector, ϕ::StateVector)
    _vec_vec_validation(ψ, ϕ)
    return _sv_sv_product(ψ, ϕ)
end
dot(ψ::StateVector, ϕ::StateVector) = inner(ψ, ϕ)
*(ψ::StateVector, ϕ::StateVector) = inner(ψ, ϕ)

function _sv_sv_product(ψ::StateVector, ϕ::StateVector)
    return contract(
        reshape(tensor(ψ), length(tensor(ψ))),
        reshape(tensor(ϕ), length(tensor(ϕ))),
        1,
        1,
        !isconj(ψ),
        isconj(ϕ),
    )[]
end

### Inner products with StateOperators
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
function inner(ψ::StateVector, ϕs::Union{StateVector,StateOperator}...)
    # Checks 
    _inner_validation(ψ, ϕs...)
    _vec_op_vec_validation(ψ, ϕs[end], ϕs[begin:(end-1)]...)

    # Create a StateVector from the cache 
    ϕ = GStateTensor(1, cache(outerdims(ϕs[end-1]), tensor(ϕs[end-1]), tensor(ϕs[end])))

    # Contraction
    applyso!(ϕ, ϕs[end-1], ϕs[end])
    for i in Base.range(length(ϕs)-2, 1, step = -1)
        ϕnew = GStateTensor(1, cache(outerdims(ϕs[i]), tensor(ϕs[i])))
        _so_sv_product!(ϕnew, ϕs[i], ϕ)
        ϕ = ϕnew
    end
    return _sv_sv_product(ψ, ϕ)
end
