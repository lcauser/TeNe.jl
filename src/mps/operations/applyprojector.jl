#= 
    Functionality for multiplying by MPSProjectors 
=#

### Applying an MPSProjector to an MPS 
"""
    applympo(O::MPSProjector, ψ::MPS; kwargs...)
    applympo(ψ::MPS, O::MPSProjector; kwargs...)
    *(O::MPSProjector, ψ::MPS)
    *(ψ::MPS, O::MPSProjector)

# Examples 

```julia-repl
julia> ψ = randommps(2, 20, 3);
julia> O = MPSProjector(randommps(2, 20, 4), randommps(2, 20, 5));
julia> ψ = O * ψ;
"""
function applympo(O::MPSProjector, ψ::MPS)
    _op_vec_validation(O, ψ)
    return (O.λ * inner(O.ϕ, ψ)) * O.ψ 
end
function applympo(ψ::MPS, O::MPSProjector)
    _op_vec_validation(adjoint(O), ψ)
    return (O.λ * inner(ψ, O.ψ)) * O.ϕ 
end
*(O::MPSProjector, ψ::MPS) = applympo(O, ψ)
*(ψ::MPS, O::MPSProjector) = applympo(ψ, O)


### Applying an MPSProjector to an MPSProjector
"""
    applympo(O1::MPSProjector, O2::MPSProjector; kwargs...)
    *(O1::MPSProjector, O2::MPSProjector)

# Examples 

```julia-repl
julia> O1 = MPSProjector(randommps(2, 20, 4), randommps(2, 20, 5));
julia> O2 = MPSProjector(randommps(2, 20, 4), randommps(2, 20, 5));
julia> O3 = O1 * O2;
"""
function applympo(O1::MPSProjector, O2::MPSProjector)
    _op_op_validation(O1, O2)
    return MPSProjector(O1.ψ, O2.ϕ, false, inner(O1.ϕ, O2.ψ) * O1.λ * O2.λ)
end
*(O1::MPSProjector, O2::MPSProjector) = applympo(O1, O2)


### Applying an MPSProjector to an MPO 
"""
    applympo(O1::MPSProjector, O2::MPO; kwargs...)
    applympo(O1::MPO, O2::MPSProjector; kwargs...)
    *(O1::MPSProjector, O2::MPO; kwargs...)
    *(O1::MPO, O2::MPSProjector; kwargs...)

Apply an MPO to an MPS Projector. 

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
function applympo(O1::MPSProjector, O2::MPO; kwargs...)
    ϕ = applympo(O1.ϕ, O2; kwargs...)
    return MPSProjector(O1.ψ, ϕ, false, O1.λ)
end

function applympo(O1::MPO, O2::MPSProjector; kwargs...)
    ψ = applympo(O1, O2.ψ; kwargs...)
    return MPSProjector(ψ, O2.ϕ, false, O2.λ)
end
*(O1::MPSProjector, O2::MPO) = applympo(O1, O2; cutoff=_TeNe_cutoff)
*(O1::MPO, O2::MPSProjector) = applympo(O1, O2; cutoff=_TeNe_cutoff)