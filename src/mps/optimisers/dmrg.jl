#=
    The density-matrix renormalization group on an MPS ansatz.
=#

export dmrg
"""
    dmrg(ψ::MPS, Hs::Union{MPO, MPS}...; kwargs...)

Run the DMRG algorithm to find extremal eigenvalues of the sum of operators `Hs`.

# Optional Keyword Arguments

    - `nsite::Int=2`: The number of sequential tensors to update in an optisation.
    Set to `1` for variational MPS, and `2` for DMRG. Beware that making this bigger
    will result in a slower optimisation.
    - `cutoff::Float64=0.0`: Truncation criteria to reduce the bond dimension.
    Good values range from 1e-8 to 1e-14. Not needed if `nsites=1`.
    - `mindim::Int=1`: Mininum dimension for truncated. Not needed if `nsites=1`.
    - `maxdim::Int=0`: Maximum bond dimension for truncation. Set to 0 to have
    no limit. Not needed if `nsites=1`.
    - `tol::Float64=1e-8`: Convergence tolerance per sweep.
    - `which::Symbol=:SR`: Which eigenvalue to target? Use `:SR` for the smallest real
    component, `:LR` for the largest real component.
    - `ishermitian::Bool=true`: Is the operator Hermitian?
    - `kryloviters::Int=2`: The number of Krylov iterations per update.
    - `krylovdim::Int=3`: The size of the Krylov subspace per update.
    - `krylovtol::Float=1e-14`: The tolerance for the Krylov iterator.
    - `verbose::Boo=true`: Output the information during optimisation?
"""
function dmrg(ψ::MPS, Hs::Union{MPO,MPS,MPSProjector}...; kwargs...)
    # Construct the projections 
    projψs = MPSProjection[]
    λproj::Number = get(kwargs, :λproj, 100.0)
    for H in Hs
        if ismpo(H) || ismpsprojector(H)
            push!(projψs, ProjMPS(ψ, H, ψ; sym = true))
        else
            push!(projψs, ProjMPSSquared(H, ψ; λ = λproj))
        end
    end

    # Create the updater 
    update = MPSUpdateDMRG(;
        which = get(kwargs, :which, :SR),
        ishermitian = get(kwargs, :ishermitian, true),
        krylovtol = Float64(get(kwargs, :krylovtol, 1e-14)),
        kryloviters = get(kwargs, :kryloviters, 2),
        krylovdim = get(kwargs, :krylovdim, 3),
    )

    # Create the optimiser 
    optim = MPSOptimiser(
        ψ,
        projψs,
        update,
        MPSObjectiveDMRG();
        nsites = get(kwargs, :nsites, 2),
        verbose = get(kwargs, :verbose, true),
        cutoff = Float64(get(kwargs, :cutoff, _TeNe_cutoff)),
        mindim = get(kwargs, :mindim, 0),
        maxdim = get(kwargs, :maxdim, 0),
        tol = Float64(get(kwargs, :tol, 1e-8)),
    )

    # Run the optimiser 
    sweep!(optim, get(kwargs, :maxsweeps, 0))

    return optim.costs[end], optim
end


### The updater
struct MPSUpdateDMRG <: MPSUpdate
    kryloviters::Int
    krylovdim::Int
    ishermitian::Bool
    tol::Float64
    which::Symbol
end

function MPSUpdateDMRG(;
    which::Symbol = :SR,
    kryloviters::Int = 2,
    krylovdim::Int = 3,
    ishermitian::Bool = false,
    krylovtol::Float64 = 1e-14,
)
    return MPSUpdateDMRG(kryloviters, krylovdim, ishermitian, krylovtol, which)
end

function update(updater::MPSUpdateDMRG, optim::MPSOptimiser, A)
    f(x) = _dmrg_applyH(optim, x)
    _, vecs = eigsolve(
        f,
        A,
        1,
        updater.which,
        maxiter = updater.kryloviters,
        krylovdim = updater.krylovdim,
        ishermitian = updater.ishermitian,
        tol = updater.tol,
    )
    return vecs[1]
end

function _dmrg_applyH(optim::MPSOptimiser, A)
    Anew = similar(A) .= 0
    for projψ in optim.projψs
        Anew .+= product(projψ, A, optim.dir)
    end
    return Anew
end

### DMRG Objective function
struct MPSObjectiveDMRG <: MPSObjective end
function measure(::MPSObjectiveDMRG, optim::MPSOptimiser)
    cost = 0.0
    for projψ in optim.projψs
        cost += inner(projψ)
    end
    return real(cost)
end
