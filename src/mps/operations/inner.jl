#=
    Inner products with MPS and MPOs 
=#

import Base.* 
export dot, inner 

### Inner product of MPS 
"""
    inner(ψ::MPS, ϕ::MPS)
    dot(ψ::MPS, ϕ::MPS)
    *(ψ::MPS, ϕ::MPS)

Calculate the inner product of two MPSs `ψ` and `ϕ`.

# Examples

```jldoctest
julia> ψ = productmps(10, [1, 1]);
julia> ϕ = productmps(10, [1, 0]);
julia> ϕ * ψ
1
```
"""
function inner(ψ::MPS, ϕ::MPS)
    # Checks 
    if !_mps_mps_product_check(ψ::MPS, ϕ::MPS)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    return _mps_mps_product(ψ, ϕ)
end
dot(ψ::MPS, ϕ::MPS) = inner(ψ, ϕ)
import Base.*
*(ψ::MPS, ϕ::MPS) = inner(ψ, ϕ)


function _mps_mps_product(ψ::MPS, ϕ::MPS)
    # Type info...
    T= _promote_tensor_eltype(ψ, ϕ)

    # Contract the network...
    block = cache(T, (size(ψ[begin], 1), size(ϕ[begin], 1)), 2, 1) .= 1
    for i in eachindex(ψ)
        # Contract the new block 
        block_new = contract(block, ψ[i], 1, 1, false, !isconj(ψ))
        block = contract(block_new, ϕ[i], (1, 2), (1, 2), false, isconj(ϕ))
    end
    return block[]
end

function _mps_mps_product_check(ψ::MPS, ϕ::MPS)
    if length(ψ) != length(ϕ)
        return false
    end
    for i in eachindex(ψ)
        if dim(ψ, i) != dim(ϕ, i)
            return false
        end
    end
    return true
end

### Inner product of MPS with MPOs
"""
    inner(ψ::MPS, O::MPO, ϕ::MPS)
    inner(ψ::MPS, O1::MPO, O2::MPO, ϕ::MPS)

Calculate the expectation of a string of operators `Os` with respect to MPSs `ψ` and `ϕ`.
"""
function inner(ψ::MPS, ϕs::Union{MPS, MPO}...)
    # Checks 
    for i = 1:length(ϕs)-1
        if !ismpo(ϕs[i])
            throw(ArgumentError("The inner terms in the braket must be MPOs."))
        end
    end
    if !ismps(ϕs[end])
        throw(ArgumentError("The last term in the MPS must be an MPO."))
    end
    if !_mps_mpo_mps_product_check(ψ, ϕs[end], ϕs[begin:end-1]...)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    return _mps_mpo_mps_product(ψ, ϕs[end], ϕs[begin:end-1]...)
end

function _mps_mpo_mps_product(ψ::MPS, ϕ::MPS, Os::MPO...)
    # Type info...
    T = Base.promote_op(*, eltype(ψ), eltype(ϕ), eltype.(Os)...)

    # Contraction 
    dims = (size(ψ[begin], 1), map(O->size(O[begin], 1), Os)..., size(ϕ[begin], 1))
    block = cache(T, dims, 2, 1) .= 1
    for i in eachindex(ψ)
        block = contract(block, ψ[i], 1, 1, false, !isconj(ψ))
        for j in eachindex(Os)
            block = contract(block, Os[j][i], (1, ndims(block)-1), (1, innerind(Os[j])), false, isconj(Os[j]))
        end
        block = contract(block, ϕ[i], (1, ndims(block)-1), (1, 2), false, isconj(ϕ))
    end

    return block[]
end

function _mps_mpo_mps_product_check(ψ::MPS, ϕ::MPS, Os::MPO...)
    if length(ψ) != length(ϕ)
        return false 
    end
    for i in eachindex(Os)
        if length(Os[i]) != length(ψ)
            return false 
        end
    end
    for i in eachindex(ψ)
        if dim(ψ, i) != innerdim(Os[1], i)
            return false 
        end
        if dim(ϕ, i) != outerdim(Os[end], i)
            return false
        end
        for j in Base.OneTo(length(Os)-1)
            if outerdim(Os[j], i) != innerdim(Os[j+1], i)
                return false
            end
        end
    end
    return true
end

### Inner products MPOs with StateVectors
"""
    inner(ψ::StateVector, O::MPO, ϕ::StateVector)
    inner(ψ::StateVector, O1::MPO, O2::MPO, ϕ::StateVector)

Calculate the expectation of a string of operators `Os` with respect to StateVectors `ψ` and `ϕ`.
"""
function inner(ψ::StateVector, ϕs::Union{StateVector, MPO}...)
    # Checks 
    for i = 1:length(ϕs)-1
        if !ismpo(ϕs[i])
            throw(ArgumentError("The inner terms in the braket must be MPOs."))
        end
    end
    if !isstatevector(ϕs[end])
        throw(ArgumentError("The last term in the MPS must be an StateVector."))
    end
    if !_sv_mpo_sv_product_check(ψ, ϕs[end], ϕs[begin:end-1]...)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    return _sv_mpo_sv_product(ψ, ϕs[end], ϕs[begin:end-1]...)
end

function _sv_mpo_sv_product(ψ::StateVector, ϕ::StateVector, Os::MPO...)
    for i in reverse(eachindex(Os))
        ϕ′ = GStateTensor(1, dim(ψ), cache(innerdims(Os[i]), tensor(ψ), tensor(ϕ), Os[i][begin]))
        _mpo_sv_product!(ϕ′, Os[i], ϕ)
        ϕ = ϕ′
    end
    return _sv_sv_product(ψ, ϕ)
end

function _sv_mpo_sv_product_check(ψ::StateVector, ϕ::StateVector, Os::MPO...)
    if length(ψ) != length(ϕ)
        return false 
    end
    for i in eachindex(Os)
        if length(Os[i]) != length(ψ)
            return false 
        end
    end
    for i in Base.OneTo(length(ψ))
        if dim(ψ, i) != innerdim(Os[1], i)
            return false 
        end
        if dim(ϕ, i) != outerdim(Os[end], i)
            return false
        end
        for j in Base.OneTo(length(Os)-1)
            if outerdim(Os[j], i) != innerdim(Os[j+1], i)
                return false
            end
        end
    end
    return true
end