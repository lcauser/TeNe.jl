#=
    Applying MPOs to MPS and StateTensors...
=#
export applympo, applympo!

### Applying an MPO to an MPS
# TODO: add densitymatrix method and variational method...
"""
    applympo(O::MPO, ψ::MPS; kwargs...)
    applympo(ψ::MPS, O::MPO; kwargs...)
    *(O::MPO, ψ::MPS)
    *(ψ::MPS, O::MPO)
    

Apply MPO `O` to MPS `ψ`.

# Optional Keyword Arguments

    - `alg=:naive`: The algorithm to carry out the multiplication. Use :naive 
      for an exact application, followed by SVD truncation. Use :zipup for quicker
      but less precise applications of the MPO. Use :densitymatrix for a one-shot method 
      with high accuracy. Use :variational for the slowest, but must optimal application.
    - `cutoff::Float64=0.0`: Truncation criteria to reduce the bond dimension.
      Good values range from 1e-8 to 1e-14.
    - `mindim::Int=1`: Minimum dimension for the truncation.
    - `maxdim::Int=0`: Maximum bond dimension for truncation. Set to 0 to have
      no limit.
"""
function applympo(O::MPO, ψ::MPS; alg = :naive, kwargs...)
    # Checks on the arguments
    _op_vec_validation(O, ψ)

    # Create a new MPS to store the result 
    ϕ = MPS(dim(ψ), length(ψ); T = _promote_tensor_eltype(O, ψ))

    # Run the algorithm
    if alg==:naive
        _mpo_mps_naive!(ϕ, O, ψ; kwargs...)
    elseif alg==:zipup
        _mpo_mps_zipup!(ϕ, O, ψ; kwargs...)
    else
        throw(ArgumentError("The algorithm $(alg) is unknown."))
    end
    return ϕ
end
applympo(ψ::MPS, O::MPO; kwargs...) = applympo(transpose(O), ψ; kwargs...)
*(O::MPO, ψ::MPS) = applympo(O, ψ; cutoff = _TeNe_cutoff)
*(ψ::MPS, O::MPO) = applympo(ψ, O; cutoff = _TeNe_cutoff)

# Naive method; do the contraction exactly and then truncate
function _mpo_mps_naive!(ϕ::MPS, O::MPO, ψ::MPS; kwargs...)
    # Move canonical centre of both to the first site 
    movecenter!(O, firstindex(O))
    movecenter!(ψ, firstindex(ψ))

    # Do the contraction 
    ϕ.center = firstindex(ϕ)
    for i in eachindex(ϕ)
        tensor = contract(O[i], ψ[i], outerind(O), 2, isconj(O), isconj(ψ))
        tensor = _permutedims(tensor, (1, 4, 2, 3, 5))
        tensor = reshape(
            tensor,
            (
                size(tensor, 1)*size(tensor, 2),
                size(tensor, 3),
                size(tensor, 4)*size(tensor, 5),
            ),
        )
        ϕ[i] = copy(tensor)
        movecenter!(ϕ, i; kwargs...)
    end
    movecenter!(ϕ, firstindex(ϕ); kwargs...)
end

# Zip-up method: see ``E.M. Stoudenmire, Steven R. White, New J. Phys. 12, 055026 (2010)``
function _mpo_mps_zipup!(ϕ::MPS, O::MPO, ψ::MPS; kwargs...)
    # Move canonical centre of both to the first site 
    movecenter!(O, firstindex(O))
    movecenter!(ψ, firstindex(ψ))

    # Do the contraction
    block = cache(eltype(ϕ), (1, size(O[begin], 1), size(ψ[begin], 1)), 2, 1) .= 1
    for i in eachindex(O)
        block = contract(block, O[i], 2, 1, false, isconj(O))
        block = contract(block, ψ[i], (2, 1 + outerind(O)), (1, 2), false, isconj(ψ))
        if i < length(O)
            U, S, block = tsvd(block, (3, 4); kwargs...)
            block = contract(S, block, 2, 1)
            ϕ[i] = U
        else
            block = reshape(block, (size(block, 1), size(block, 2), 1))
            ϕ[i] = copy(block)
        end
    end
    ϕ.center = length(ϕ)
    movecenter!(ϕ, firstindex(ϕ); kwargs...)
end

### Applying an MPO to an MPO

"""
    applympo(O1::MPO, O2::MPO; kwargs...)

Apply MPO `O1` to MPO `O2`.

# Optional Keyword Arguments

    - `alg=:naive`: The algorithm to carry out the multiplication. Use :naive 
      for an exact application, followed by SVD truncation. Use :zipup for quicker
      but less precise applications of the MPO. Use :densitymatrix for a one-shot method 
      with high accuracy. Use :variational for the slowest, but must optimal application.
    - `cutoff::Float64=0.0`: Truncation criteria to reduce the bond dimension.
      Good values range from 1e-8 to 1e-14.
    - `mindim::Int=1`: Minimum dimension for the truncation.
    - `maxdim::Int=0`: Maximum bond dimension for truncation. Set to 0 to have
      no limit.
"""
function applympo(O1::MPO, O2::MPO; alg = :naive, kwargs...)
    # Checks on the arguments
    _op_op_validation(O1, O2)

    # Create MPO to store the result 
    O = MPO(dim(O1), length(O2); T = _promote_tensor_eltype(O1, O2))

    # Run the algorithm
    if alg==:naive
        _mpo_mpo_naive!(O, O1, O2; kwargs...)
    elseif alg==:zipup
        _mpo_mpo_zipup!(O, O1, O2; kwargs...)
    else
        throw(ArgumentError("The algorithm $(alg) is unknown."))
    end
    return O
end
*(O1::MPO, O2::MPO) = applympo(O1, O2; cutoff = _TeNe_cutoff)

function _mpo_mpo_naive!(O::MPO, O1::MPO, O2::MPO; kwargs...)
    # Move canonical centre of both to the first site 
    movecenter!(O1, 1)
    movecenter!(O2, 1)

    # Do the contraction 
    O.center = firstindex(O)
    for i in eachindex(O)
        tensor = contract(O1[i], O2[i], outerind(O1), innerind(O2), isconj(O1), isconj(O2))
        tensor = _permutedims(tensor, (1, 4, 2, 5, 3, 6))
        tensor = reshape(
            tensor,
            (
                size(tensor, 1)*size(tensor, 2),
                size(tensor, 3),
                size(tensor, 4),
                size(tensor, 5)*size(tensor, 6),
            ),
        )
        O[i] = copy(tensor)
        movecenter!(O, i; kwargs...)
    end
    movecenter!(O, firstindex(O); kwargs...)
end

function _mpo_mpo_zipup!(O::MPO, O1::MPO, O2::MPO; kwargs...)
    # Move canonical centre of both to the first site 
    movecenter!(O1, firstindex(O1))
    movecenter!(O2, firstindex(O2))

    # Do the contraction
    block = cache(eltype(O), (1, size(O1[begin], 1), size(O2[begin], 1)), 2, 1) .= 1
    for i in eachindex(O)
        block = contract(block, O1[i], 2, 1, false, isconj(O1))
        block = contract(
            block,
            O2[i],
            (2, 1 + outerind(O1)),
            (1, innerind(O2)),
            false,
            isconj(O2),
        )
        if i < length(O)
            U, S, block = tsvd(block, (3, 5); kwargs...)
            block = contract(S, block, 2, 1)
            O[i] = U
        else
            block = reshape(block, (size(block, 1), size(block, 2), size(block, 4), 1))
            O[i] = copy(block)
        end
    end
    O.center = length(O)
    movecenter!(O, firstindex(O); kwargs...)
end

### Applying an MPO to a state vector
"""
    applympo(O::MPO, ψ::StateVector)
    applympo(ψ::StateVector, O::MPO)
    *(O::MPO, ψ::StateVector)
    *(ψ::StateVector, O::MPO)

Apply MPO `O` to StateVector `ψ`.

# Examples 

```julia-repl
julia> ψ = randomsv(2, 10);
julia> O = productmpo(10, [0 1; 1 0]);
julia> ϕ = O * ψ;
"""
function applympo(O::MPO, ψ::StateVector)
    _op_vec_validation(O, ψ)
    ϕ = GStateTensor(1, promote_tensor(innerdims(O), tensor(ψ), O[begin]))
    _mpo_sv_product!(ϕ, O, ψ)
    return ϕ
end
applympo(ψ::StateVector, O::MPO) = applympo(transpose(O), ψ)
*(O::MPO, ψ::StateVector) = applympo(O, ψ)
*(ψ::StateVector, O::MPO) = applympo(ψ, O)


"""
    applympo!(ϕ::StateVector, O::MPO, ψ::StateVector)
    applympo!(ϕ::StateVector, ψ::StateVector, O::MPO)

Apply MPO `O` to StateVector `ψ`, and store the result in `ϕ`.

# Examples 

```julia-repl
julia> ψ = randomsv(2, 10);
julia> O = productmpo(10, [0 1; 1 0]);
julia> ϕ = StateVector(2, 10);
julia> applympo!(ϕ, O, ψ);
"""
function applympo!(ϕ::StateVector, O::MPO, ψ::StateVector)
    _op_vec_validation(ϕ, O, ψ)
    _mpo_sv_product!(ϕ, O, ψ)
end
applympo!(ϕ::StateVector, ψ::StateVector, O::MPO) = applyMPO(ϕ, transpose(O), ψ)

function _mpo_sv_product!(ϕ::StateVector, O::MPO, ψ::StateVector)
    # Apply the first tensor 
    ten = reshape(tensor(ψ), (size(tensor(ψ))..., 1))
    ten = contract(ten, O[1], (ndims(ten), 1), (1, outerind(O)), isconj(ψ), isconj(O))

    # Loop through the remaining tensors 
    for i in range(firstindex(O)+1, lastindex(O))
        ten = contract(ten, O[i], (ndims(ten), 1), (1, outerind(O)), false, isconj(O))
    end
    ten = reshape(ten, size(ten)[begin:(end-1)])
    tensor(ϕ) .= isconj(ϕ) ? conj.(ten) : ten
end


### Applying an MPO to a StateOperator
"""
    applympo(O1::MPO, O2::StateOperator)
    applympo(O1::StateOperator, O2::MPO)

Multiply a StateOperator by an MPO.
"""
function applympo(O1::MPO, O2::StateOperator)
    _op_op_validation(O1, O2)
    O = GStateTensor(2, promote_tensor(_mpo_so_dims(O1, O2), O1, O2))
    _mpo_so_product!(O, O1, O2)
    return O
end
function applympo(O1::StateOperator, O2::MPO)
    _op_op_validation(O1, O2)
    O = GStateTensor(2, dim(O2), promote_tensor(_mpo_so_dims(O1, O2, true), O1, O2))
    _mpo_so_product!(O, O2, O1; tp = true)
    return O
end
*(O1::MPO, O2::StateOperator) = applympo(O1, O2)
*(O1::StateOperator, O2::MPO) = applympo(O1, O2)

"""
    applympo!(O::StateOperator, O1::MPO, O2::StateOperator)
    applympo!(O::StateOperator, O1::StateOperator, O2::MPO)

Multiply a StateOperator by an MPO, and store the result in `O`.
"""
function applympo!(O::StateOperator, O1::MPO, O2::StateOperator)
    _op_op_validation(O, O1, O2)
    _mpo_so_product!(O, O1, O2)
end

function applympo!(O::StateOperator, O1::StateOperator, O2::MPO)
    _op_op_validation(O, O1, O2)
    _mpo_so_product!(O, transpose(O2), transpose(O1); tp = true)
end

function _mpo_so_product!(O::StateOperator, O1::MPO, O2::StateOperator; tp::Bool = false)
    # Apply the first tensor 
    ten = reshape(tensor(O2), (size(tensor(O2))..., 1))
    ten = contract(
        ten,
        O1[1],
        (ndims(ten), istranspose(O2) ? 2 : 1),
        (1, istranspose(O1) ? 2 : 3),
        isconj(O2),
        isconj(O1),
    )

    # Apply intermediate tensors
    # Loop through the remaining tensors 
    for i in range(firstindex(O1)+1, lastindex(O1))
        ten = contract(
            ten,
            O1[i],
            (ndims(ten), istranspose(O2) ? i+1 : i),
            (1, istranspose(O1) ? 2 : 3),
            false,
            isconj(O1),
        )
    end

    # Manipulate to fit O
    ten = reshape(ten, size(ten)[begin:(end-1)])
    perms = _applyso_perm_dims(length(O))
    if istranspose(O) != tp
        perms = Tuple(map(j->isodd(j) ? perms[j+1] : perms[j-1], Base.eachindex(perms)))
    end
    permutedims!(tensor(O), ten, perms)
    if isconj(O)
        tensor(O) .= conj.(tensor(O))
    end
end

function _mpo_so_dims(O1::MPO, O2::StateOperator, reverse::Bool = false)
    if reverse
        dims = Tuple(
            map(
                j->isodd(j) ? size(tensor(O2), istranspose(O2) ? j+1 : j) :
                   size(O1[cld(j, 2)], outerind(O1)),
                Base.OneTo(2*length(O1)),
            ),
        )
    else
        dims = Tuple(
            map(
                j->isodd(j) ? size(O1[cld(j, 2)], innerind(O1)) :
                   size(tensor(O2), istranspose(O2) ? j-1 : j),
                Base.OneTo(2*length(O1)),
            ),
        )
    end
    return dims
end
