function benchmark_stateoptimiser()
    # Parameters for check 
    N = 20
    J = 0.3
    h = 1.0

    # Create Hamiltonian 
    qu = Qubits()
    H = OpList(qu, N)
    for i = 1:N
        add!(H, "x", i, -h)
    end
    for i = 1:(N-1)
        add!(H, ["z", "z"], [i, i+1], -J)
    end
    H = MPO(H)

    # Do DMRG 
    ψ = randommps(2, N, 1)
    energy, optim = dmrg(ψ, H; cutoff = 1e-12, verbose = false)

    # Optimise circuit 
    ϕ = productmps(qu, ["dn" for _ = 1:N])
    U = randombwcircuit(2, N, 3)
    optim = StateOptimiser(ψ, U, ϕ; tol = 1e-5, verbose = 0)
end
