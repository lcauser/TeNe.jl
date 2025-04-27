#=
    Inner products with MPS and MPOs 
=#

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
    _vec_vec_validation(ψ, ϕ)
    return _mps_mps_product(ψ, ϕ)
end
dot(ψ::MPS, ϕ::MPS) = inner(ψ, ϕ)
*(ψ::MPS, ϕ::MPS) = inner(ψ, ϕ)


function _mps_mps_product(ψ::MPS, ϕ::MPS)
    # Type info...
    T = _promote_tensor_eltype(ψ, ϕ)

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
function inner(ψ::MPS, ϕs::Union{MPS,MPO}...)
    # Checks 
    _inner_validation(ψ, ϕs...)
    _vec_op_vec_validation(ψ, ϕs[end], ϕs[begin:(end-1)]...)
    return _mps_mpo_mps_product(ψ, ϕs[end], ϕs[begin:(end-1)]...)
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
            block = contract(
                block,
                Os[j][i],
                (1, ndims(block)-1),
                (1, innerind(Os[j])),
                false,
                isconj(Os[j]),
            )
        end
        block = contract(block, ϕ[i], (1, ndims(block)-1), (1, 2), false, isconj(ϕ))
    end

    return block[]
end

### Inner products of MPSs with both/either MPSProjectors and MPOs 
function inner(ψ::MPS, ϕs::Union{MPS,MPO,MPSProjector}...)
    # Checks 
    _inner_validation(ψ, ϕs...)
    _vec_op_vec_validation(ψ, ϕs[end], ϕs[begin:(end-1)]...)

    prod = 1.0
    prod_string = Union{MPS,MPO}[ψ]
    for ϕ in ϕs
        if ismpsprojector(ϕ)
            push!(prod_string, ϕ.ψ)
            prod *= ϕ.λ * inner(prod_string...)
            prod_string = Union{MPS,MPO}[ϕ.ϕ]
        else
            push!(prod_string, ϕ)
        end
    end
    prod *= inner(prod_string...)
    return prod
end



### Inner products MPOs with StateVectors
"""
    inner(ψ::StateVector, O::MPO, ϕ::StateVector)
    inner(ψ::StateVector, O1::MPO, O2::MPO, ϕ::StateVector)

Calculate the expectation of a string of operators `Os` with respect to StateVectors `ψ` and `ϕ`.
"""
function inner(ψ::StateVector, ϕs::Union{StateVector,MPO}...)
    # Checks 
    _inner_validation(ψ, ϕs...)
    _vec_op_vec_validation(ψ, ϕs[end], ϕs[begin:(end-1)]...)
    return _sv_mpo_sv_product(ψ, ϕs[end], ϕs[begin:(end-1)]...)
end

function _sv_mpo_sv_product(ψ::StateVector, ϕ::StateVector, Os::MPO...)
    for i in reverse(eachindex(Os))
        ϕ′ = GStateTensor(
            1,
            dim(ψ),
            cache(innerdims(Os[i]), tensor(ψ), tensor(ϕ), Os[i][begin]),
        )
        _mpo_sv_product!(ϕ′, Os[i], ϕ)
        ϕ = ϕ′
    end
    return _sv_sv_product(ψ, ϕ)
end


### Inner products of MPS with operator lists 
"""
    inner(ψ::MPS, O::OpList, ϕ::MPS)

Measure the expectation value for each operator in an OpList with respect to MPSs
`ψ` and `ϕ`.
"""
function inner(ψ::MPS, O::OpList, ϕ::MPS)
    # Validation 
    _vec_op_vec_validation(ψ, ϕ, O)

    # Build the environment 
    proj = ProjMPS(ψ, ϕ)

    # Calculate the expectation for each operator in the list 
    expectations = zeros(_promote_tensor_eltype(O, ψ, ϕ), length(O.ops))
    for i in eachindex(O.ops)
        expectations[i] = _inner(proj, O.lt, O.ops[i], O.sites[i]) * O.coeffs[i]
    end
    return expectations
end

function _inner(proj::ProjMPS, lt::LatticeTypes, ops::Vector{String}, sites::Vector{Int})
    # Fetch the blocks
    left = leftblock(proj, sites[begin]-1)
    right = rightblock(proj, sites[end]+1)

    # Contract with the center sites and the relavent operator 
    ctr = 1
    for i in Base.range(sites[begin], sites[end])
        left = contract(
            left,
            proj.objects[begin][i],
            1,
            1,
            false,
            !isconj(proj.objects[begin]),
        )
        if ctr <= length(sites) && i == sites[ctr]
            left = contract(left, op(lt, ops[ctr]), 2, 1)
            left = permutedim(left, 3, 2)
            ctr += 1
        end
        left = contract(
            left,
            proj.objects[end][i],
            (1, 2),
            (1, 2),
            false,
            isconj(proj.objects[end]),
        )
    end
    return contract(left, right, (1, 2), (1, 2))[]
end
