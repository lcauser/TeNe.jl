#=
    The state preparation algorithm; given an initial MPS, a circuit ansatz
    and a final state, the algorithm will use the a polar decompositon
    to optimise the circuit sequentially.
=#

mutable struct StateOptimiser
    # Tensor networks to store 
    ϕ::MPS # Target state 
    circuit::Circuit
    ψ::MPS # Initial state 
    proj::ProjMPSCircuit

    # Update properties 
    sweeps::Int
    tol::Float64
    dir::Bool
    verbose::Int 
    costs::Vector{Float64}
end

export StateOptimiser

"""
    StateOptimiser(ϕ::MPS, circuit::Circuit, ψ::MPS; kwargs...)

Optimise a unitary circuit that evolves `ψ` to `ϕ`.

# Optional keyword arguments 
    - `minsweeps::Int=1`: The minimum number of sweeps to perform.
    - `maxsweeps::Int=0`: The maximum number of sweeps to do. Set to zero for no limit.
    - `tol::Float64=1e-8`: The convergence tolerance.
    - `cutoff::Float64=$(_TeNe_cutoff)`: The SVD cutoff to use when estimating the environment.
    - `maxdim::Int=0`: The minimum bond dimension of the MPSs for estimating the environment.
    - `verbose::Int=1`: Output optimisation information? Set to 0 for no, 1 for after each sweep, 2 after each row.
"""
function StateOptimiser(
    ϕ::MPS, circuit::Circuit, ψ::MPS;
    minsweeps::Int=1,
    maxsweeps::Int=0,
    tol::Float64=1e-8,
    cutoff::Float64=_TeNe_cutoff,
    maxdim::Int=0,
    verbose::Int=1
)   
    proj = ProjMPSCircuit(ϕ, circuit, ψ; cutoff=cutoff, maxdim=maxdim)
    optim = StateOptimiser(ϕ, circuit, ψ, proj, 0, tol, false, verbose, Float64[abs(product(proj))^2])
    sweep!(optim, maxsweeps; minsweeps=minsweeps)
    return optim
end

function sweep!(optim::StateOptimiser, sweeps::Int=0; minsweeps::Int=1)
    # Fetch properties 
    circuit = optim.circuit
    proj = optim.proj

    # Sweep through 
    iters = 0
    while true
        # Fetch the rows to sweep through
        rows = Base.range(firstindex(circuit.layers), lastindex(circuit.layers))
        rows = optim.dir ? reverse(rows) : rows
        
        for row in rows 
            # Move the center to the correct row 
            movecenter!(proj, row)

            # Loop through the gates 
            layer = getlayer(proj, proj.depth_center)
            for g in Base.OneTo(length(layer.gates))
                # Move center 
                movecenter!(proj, row, g)

                # Find the optimal gate the place there 
                gate = creategate(conj(project(proj)))
                makeunitary!(gate)
                layer.gates[g].gate .= gate.gate
            end
            if optim.verbose == 2
                cost = product(proj) 
                @printf("iter=%d, row=%d overlap=%.8E \n", optim.sweeps+1, row, real(cost))
            end
        end

        # Check convergence
        cost = abs(product(proj))^2
        push!(optim.costs, cost)
        iters += 1
        optim.sweeps += 1

        if optim.verbose == 1
            @printf("iter=%d, overlap=%.8E \n", optim.sweeps, real(cost))
        end

        if iters > minsweeps
            if sweeps > 0 && iters >= sweeps
                break
            end
            lastcost = optim.costs[end-1]
            if 2*abs(cost - lastcost) / abs(cost + lastcost) < optim.tol
                break
            end
        end
    end
end