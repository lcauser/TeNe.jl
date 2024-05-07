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
function MPSOptimiser(ψ::MPS, projψs::Vector{<:MPSProjection}, update::MPSUpdate,
    objective::MPSObjective, observers::Vector{<:MPSObserver}=MPSObserver[];
    normalize::Bool=true, relativecheck::Bool=true, verbose::Bool=true)
    
    return MPSOptimiser(ψ, projψs, update, false, 0, normalize, objective, eltype(objective)[],
        relativecheck, verbose, observers)
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
function sweep!(optim::MPSOptimiser, sweeps::Int=1, tol::Float64=1e-6; nsites::Int=2,
    mindim::Int=0, maxdim::Int=0, cutoff::Float64=1e-12, minsweeps::Int=0)
    
    # Loop through the number of sweeps 
    iters = 0
    movecenter!(optim, !optim.dir ? lastindex(optim.ψ) : firstindex(optim.ψ))
    while true
        # Loop through the sites & update
        sites = optim.dir ? Base.range(lastindex(optim.ψ), firstindex(optim.ψ) + nsites - 1, step = -1) :
            Base.range(firstindex(optim.ψ), lastindex(optim.ψ) + 1 - nsites)
        for site in sites
            # Move the canonical center 
            movecenter!(optim, site)

            # Combine tensors
            firstsite = optim.dir ? site - nsites + 1 : site
            lastsite = optim.dir ? site : site - nsites + 1
            A = optim.ψ[firstsite]
            for i in Base.range(firstsite+1, lastsite)
                A = contract(A, optim.ψ[i], ndims(A), 1; tocache=!(i==lastsite))
            end

            # Find the updated composite tensor 
            A = update(optim.update, optim, A)

            # Restore the MPS 
            replacesites!(optim.ψ, A, firstsite, optim.normalize; cutoff=cutoff, maxdim=maxdim, mindim=mindim)
        end
        movecenter!(optim, optim.dir ? lastindex(optim.ψ) : firstindex(optim.ψ))
        optim.dir = !optim.dir

        # Measure the cost and observables 
        cost = measure(optim.objective, optim)
        push!(optim.costs, cost)

        # Output information 
        if optim.verbose 
            @printf("iter=%d, energy=%.8E, maxbonddim=%d \n", optim.sweeps,
                   real(cost), maxbonddim(optim.ψ))
        end

        # Increment counters & check for convergence 
        iters += 1
        optim.sweeps += 1 # Increment an internal counter
        if iters >= minsweeps 
            # Check sweeps 
            if sweeps != 0 && iters >= sweeps
                break
            end

            # Check cost 
            diff = optim.costs[end] - optim.cost[end-1]
            if optim.relativecheck
                diff = 2*diff / (abs(optim.costs[end]) + abs(optim.cost[end-1]))
            end
            if diff < tol
                break 
            end
        end
    end
end

### MPS Update 
function update(::MPSUpdate, ::MPSOptimiser, A::AbstractArray)
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


struct MPSUpdateDMRG <: MPSUpdate end
struct MPSObjectiveDMRG <: MPSObjective end