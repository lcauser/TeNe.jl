#=
    The time-evolving block decimation method for an MPS ansatz.
=#

export tebd 

function tebd(ψ::MPS, H::OpList, δt::Number, Δt::Number, t::Number, observers::MPSObserver...; kwargs...)
    # Create the gates 
    U = trotterize(H, δt;
        order=get(kwargs, :order, 2),
        type=get(kwargs, :type, :compressed)
    )

    # Keyword arguments for truncation
    cutoff::Float64 = get(kwargs, :cutoff, _TeNe_cutoff)
    maxdim::Int = get(kwargs, :maxdim, 0)
    mindim::Int = get(kwargs, :mindim, 0)
    normalize::Bool = get(kwargs, :normalize, true)

    # Determine time steps 
    num_small = Int(ceil(real(Δt / δt)))
    num_large = Int(ceil(real(t / Δt)))
    if !isapprox(num_small * δt, Δt)
        throw(ArgumentError("Δt must be divisible by δt."))
    end
    if !isapprox(num_large * Δt, t)
        throw(ArgumentError("t must be divisible by Δt."))
    end

    # Take initial measurements 
    _tebd_measurements!(ψ, observers...)

    # Check for energy convergence?
    energycheck::Bool = get(kwargs, :energy_convergence, false)
    energytol::Number = get(kwargs, :energy_tol, 1e-5)
    if energycheck
        lastE = sum(inner(ψ, H, ψ))
    end

    # Apply the gates 
    verbose::Bool = get(kwargs, :verbose, true)
    for i = 1:num_large
        # Apply gates multiple times
        for _ = 1:num_small 
            applygates!(U, ψ; cutoff=cutoff, maxdim=maxdim, mindim=mindim, normalize=normalize)
        end

        # Take measurements
        _tebd_measurements!(ψ, observers...)
        
        # Output
        E = sum(inner(ψ, H, ψ))
        if verbose 
            @printf("time=%.4E/%.4E, maxbonddim=%d, energy=%.6E \n", Δt*i, t, maxbonddim(ψ), real(E))
        end

        # Check convergence 
        if energycheck
            diff = abs(E) > 1e-10 ? 2*abs(lastE-E)/abs(lastE+E) : abs(lastE-E)
            if diff < energytol 
                break 
            end
            lastE = E
        end
    end
end

function _tebd_measurements!(ψ::MPS, observers::MPSObserver...)
    for observer in observers
        measure!(observer, ψ) 
    end
end