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

    # Apply the gates 
    verbose::Bool = get(kwargs, :verbose, true)
    for i = 1:num_large
        for _ = 1:num_small 
            applygates!(U, ψ; cutoff=cutoff, maxdim=maxdim, mindim=mindim, normalize=normalize)
        end
        
        if verbose 
            E = sum(inner(ψ, H, ψ))
            @printf("time=%.4E/%.4E, maxbonddim=%d, energy=%.6E \n", Δt*i, t, maxbonddim(ψ), real(E))
        end
    end
end