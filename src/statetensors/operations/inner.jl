#=
    Calculating inner products of StateVectors and StateOperators.
=#

import Base.*
export inner, dot

### Inner products of StateVectors
"""
    inner(ψ::StateVector, ϕ::StateVector)
    dot(ψ::StateVector, ϕ::StateVector)
    *(ψ::StateVector, ϕ::StateVector)

Calculate the inner product of two StateVectors `ψ` and `ϕ`.
"""
function inner(ψ::StateVector, ϕ::StateVector)
    # Checks 
    if length(ψ) != length(ϕ) || !_sv_sv_product_checkdims(ψ, ϕ)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    return _sv_sv_product(ψ, ϕ)
end
dot(ψ::StateVector, ϕ::StateVector) = inner(ψ, ϕ)
*(ψ::StateVector, ϕ::StateVector) = inner(ψ, ϕ)

function _sv_sv_product(ψ::StateVector, ϕ::StateVector)
    return contract(reshape(tensor(ψ), length(tensor(ψ))),
             reshape(tensor(ϕ), length(tensor(ϕ))),
             1, 1, !isconj(ψ), isconj(ϕ))[]
end

function _sv_sv_product_checkdims(ψ::StateVector, ϕ::StateVector)
    for i = Base.OneTo(length(ψ))
        if size(tensor(ψ), i) != size(tensor(ϕ), i)
            return false
        end
    end
    return true
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
        ϕnew = GStateTensor(1, dim(ψ), cache(size(tensor(ϕs[end])), tensor(ϕs[end]), tensor(ϕ)))
        applyso!(ϕnew, ϕs[i], ϕ)
        ϕ = ϕnew
    end
    return inner(ψ, ϕ)
end