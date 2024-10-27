@testset "DMRG" begin
    # Parameters for check 
    tol = 1e-3
    N = 20
    J = 0.9
    h = 1.0
    E = -24.065959296525854
    E2 = -23.772631597776236

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
    ψ = randommps(2, N, 1)
    energy, optim = dmrg(ψ, H; cutoff=1e-12, verbose=false)
    @test abs(energy - E) / abs(E) < tol
    ψ2 = randommps(2, N, 1)
    energy2, optim = dmrg(ψ2, H, ψ; cutoff=1e-12, verbose=false)
    @test abs(energy2 - E2) / abs(E2) < tol
end