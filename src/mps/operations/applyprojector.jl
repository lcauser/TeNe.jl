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