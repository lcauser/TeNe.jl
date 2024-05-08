#=
    The MPSOptimiser is responible for changing the parameters of an MPS to achieve a desired task.
    It sweeps through the MPS, left-to-right and then right-to-left to peform some local update that 
    aims to minimise some objective function. It aims to be quite flexible: e.g., the energy could be 
    minimised using DMRG, the energy variance using gradient descent, or time-evolution using TDVP...
=#

abstract type MPSUpdate end
abstract type MPSObjective end
abstract type MPSObserver end

mutable struct MPSOptimiser
    # Tensor networks to store 
    ψ::MPS
    projψs::Vector{<:MPSProjection}

    # The update 
    update::MPSUpdate
    dir::Bool
    sweeps::Int
    normalize::Bool
    nsites::Int
    cutoff::Float64 
    mindim::Int
    maxdim::Int 
    tol::Float64

    # Objective function 
    objective::MPSObjective
    costs::Vector{<:Number}
    relativecheck::Bool
    verbose::Bool

    # Observers 
    observers::Vector{<:MPSObserver}
end
export MPSOptimiser 

### Iniating an MPS optimiser 
"""
    MPSOptimiser(ψ::MPS, projψs::Vector{<:MPSProjection}, update::MPSUpdate,
        objective::MPSObjective, [observers::Vector{<:MPSObserver}])

Create an optimiser for an MPS `ψ` with the projections `projψs` used in the optimisation.
The `update` defines the type of local update done for the MPS, and the `objective` is the
objective function.
"""
function MPSOptimiser(ψ::MPS, projψs::Vector{<:MPSProjection}, update::MPSUpdate,
    objective::MPSObjective, observers::Vector{<:MPSObserver}=MPSObserver[];
    normalize::Bool=true, relativecheck::Bool=true, verbose::Bool=true, tol::Float64=1e-6,
    nsites::Int=2, cutoff::Float64=_TeNe_cutoff, mindim::Int=0, maxdim::Int=0)
    
    optim = MPSOptimiser(ψ, projψs, update, false, 0, normalize, nsites, cutoff, mindim,
        maxdim, tol, objective, eltype(objective)[], relativecheck, verbose, observers)
    push!(optim.costs, measure(objective, optim))
    return optim
end


### Updating the hyperparamters 
function update_hyperparameters(optim::MPSOptimiser; nsites=nothing, mindim=nothing,
    maxdim=nothing, cutoff=nothing)

    if !isnothing(nsites)
        optim.nsites = nsites
    end
    if !isnothing(maxdim)
        optim.maxdim = maxdim
    end
    if !isnothing(mindim)
        optim.mindim = mindim
    end
    if !isnothing(cutoff)
        optim.cutoff = cutoff
    end

    return optim.nsites, optim.mindim, optim.maxdim, optim.cutoff
end

### Moving canonical center of all objects 
function movecenter!(optim::MPSOptimiser, site::Int)
    movecenter!(optim.ψ, site)
    for projψ in optim.projψs
        movecenter!(projψ, site)
    end
end

### Sweeping procedure 
export sweep!
"""
    sweep!(optim::MPSOptimiser, sweeps::Int=1, tol::Float64=1e-6; kwargs...)

Perform a number of iterations `sweeps` of optimisation on the MPS using the MPSOptimiser `optim`.
Set `sweeps=0` to do until convergence is met, according to `tol`.

Optimisation hyperparameters, such as truncation criteria, will default to the last used, 
but can be changed using keyword arguments.

# Optional Keyword Arguments 

    - `nsites=nothing`: The number of sites to update in the optimisation.
    - `minsweeps::Int=0`: The minimum number of iterations to do.
    - `cutoff::Float64=$(_TeNe_cutoff)`: Truncation criteria to reduce the bond dimension.
      Good values range from 1e-8 to 1e-14.
    - `mindim=nothing`: Mininum dimension for truncated.
    - `maxdim=nothing`: Maximum bond dimension for truncation. Set to 0 to have
      no limit.
"""
function sweep!(optim::MPSOptimiser, sweeps::Int=1, tol::Float64=1e-6; nsites=nothing,
    mindim=nothing, maxdim=nothing, cutoff=nothing, minsweeps::Int=0)
    # Update hyperparameters
    nsites, mindim, maxdim, cutoff = update_hyperparameters(optim; nsites=nsites,
        mindim=mindim, maxdim=maxdim, cutoff=cutoff)

    # Loop through the number of sweeps 
    iters = 0
    mbd = maxbonddim(optim.ψ)
    movecenter!(optim, optim.dir ? lastindex(optim.ψ) : firstindex(optim.ψ))
    while true
        # Loop through the sites & update
        sites = optim.dir ? Base.range(lastindex(optim.ψ), firstindex(optim.ψ) + nsites - 1, step = -1) :
            Base.range(firstindex(optim.ψ), lastindex(optim.ψ) + 1 - nsites)
        for site in sites
            # Move the canonical center 
            movecenter!(optim, site)

            # Combine tensors
            firstsite = optim.dir ? site - nsites + 1 : site
            lastsite = optim.dir ? site : site + nsites - 1
            A = optim.ψ[firstsite]
            for i in Base.range(firstsite+1, lastsite)
                A = contract(A, optim.ψ[i], ndims(A), 1; tocache=!(i==lastsite))
            end

            # Find the updated composite tensor 
            A = update(optim.update, optim, A)

            # Restore the MPS 
            replacesites!(optim.ψ, A, site, optim.dir; normalize=optim.normalize, cutoff=cutoff,
                maxdim=maxdim, mindim=mindim)
        end
        movecenter!(optim, optim.dir ? firstindex(optim.ψ) : lastindex(optim.ψ))
        optim.dir = !optim.dir

        # Measure the cost and observables 
        cost = measure(optim.objective, optim)
        push!(optim.costs, cost)

        # Output information
        iters += 1
        optim.sweeps += 1 # Increment an internal counter 
        if optim.verbose 
            @printf("iter=%d, objective=%.8E, maxbonddim=%d \n", optim.sweeps,
                   real(cost), maxbonddim(optim.ψ))
        end

        # Check for convergence 
        if iters >= minsweeps 
            # Check sweeps 
            if sweeps != 0 && iters >= sweeps
                break
            end

            # Check cost 
            diff = optim.costs[end-1] - optim.costs[end]
            if optim.relativecheck
                diff = 2*diff / (abs(optim.costs[end]) + abs(optim.costs[end-1]))
            end
            if diff < tol && mbd == maxbonddim(optim.ψ)
                break 
            end
            mbd = maxbonddim(optim.ψ)
        end
    end
end

### MPS Update 
function update(::MPSUpdate, ::MPSOptimiser, A)
    return A
end

### MPS Objective 
function measure(::MPSObjective, ::MPSOptimiser)
    return 0.0
end
Base.eltype(::MPSObjective) = Number

### MPS Observer 
function measure(::MPSObserver, ::MPSOptimiser)
    nothing
end
