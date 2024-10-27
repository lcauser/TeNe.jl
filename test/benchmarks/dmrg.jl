function benchmark_dmrg()
    # Parameters for check 
    N = 20
    J = 0.9
    h = 1.0

    # Create Hamiltonian 
    qu = Qubits()
    H = OpList(qu, N)
    for i = 1:N
        add!(H, "x", i, -h)
    end
    for i = 1:N-1
        add!(H, ["z", "z"], [i, i+1], -J)
    end
    H = MPO(H)

    # Do DMRG 
    ψ = productmps(N, [sqrt(0.5), sqrt(0.5)], normalize=true)
    energy, optim = dmrg(ψ, H; cutoff=1e-12, verbose=false)
end