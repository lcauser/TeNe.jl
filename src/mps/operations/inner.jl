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
    if !issimilar(ψ, ϕ)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    return _mps_mps_product(ψ, ϕ)
end
dot(ψ::MPS, ϕ::MPS) = inner(ψ, ϕ)
import Base.*
*(ψ::MPS, ϕ::MPS) = inner(ψ, ϕ)


function _mps_mps_product(ψ::MPS, ϕ::MPS)
    # Type info...
    T = Base.promote_op(*, eltype(ψ), eltype(ϕ))
    conjψ = !isconj(ψ) # Not because inner product has conj on bra by default...
    conjϕ = isconj(ϕ)

    # Contract the network...
    block = cache(T, (size(ψ[begin], 1), size(ϕ[begin], 1)), 2, 1) .= 1
    for i in eachindex(ψ)
        # Contract the new block 
        block_new = contract(block, ψ[i], 1, 1, false, conjψ)
        block = contract(block_new, ϕ[i], (1, 2), (1, 2), false, conjϕ)
    end
    return block[]
end

### Inner product of MPS with MPOs
"""
    inner(ψ::MPS, O::MPO, ϕ::MPS)
    inner(ψ::MPS, O1::MPO, O2::MPO, ϕ::MPS)

Calculate the expectation of a string of operators `Os` with respect to MPSs `ψ` and `ϕ`.
"""
function inner(ψ::MPS, ϕs::Union{MPS, MPO}...)
    # Checks 
    if !issimilar(ψ, ϕs...)
        throw(ArgumentError("Arguments have properties that do not match."))
    end
    for i = 1:length(ϕs)-1
        if !ismpo(ϕs[i])
            throw(ArgumentError("The inner terms in the braket must be MPOs."))
        end
    end
    if !ismps(ϕs[end])
        throw(ArgumentError("The last term in the MPS must be an MPO."))
    end
    return _mps_mpo_mps_product(ψ, ϕs[end], ϕs[1:end-1]...)
end

function _mps_mpo_mps_product(ψ::MPS, ϕ::MPS, Os::MPO...)
    # Type info...
    T = Base.promote_op(*, eltype(ψ), eltype(ϕ), eltype.(Os)...)
    conjψ = !isconj(ψ) # Not because inner product has conj on bra by default...
    conjϕ = isconj(ϕ)
    conjOs = map(O->isconj(O), Os)
    transOs = map(O->istranspose(O), Os)

    # Contraction 
    dims = (size(ψ[begin], 1), map(O->size(O[begin], 1), Os)..., size(ϕ[begin], 1))
    block = cache(T, dims, 2, 1) .= 1
    for i in eachindex(ψ)
        block = contract(block, ψ[i], 1, 1, false, conjψ)
        for j in eachindex(Os)
            block = contract(block, Os[j][i], (1, ndims(block)-1), (1, transOs[j] ? 3 : 2), false, conjOs[j])
        end
        block = contract(block, ϕ[i], (1, ndims(block)-1), (1, 2), false, conjϕ)
    end

    return block[]
end

### Inner products MPOs with StateVectors