@testset "Circuits-qubitgates" begin
    # Single qubit gates 
    @test isapprox(tensor(QubitXGate()), [0 1; 1 0])
    @test isapprox(tensor(QubitZGate()), [1 0; 0 -1])
    @test isapprox(tensor(QubitYGate()), [0 -1im; 1im 0])
    @test isapprox(tensor(QubitIdGate()), [1 0; 0 1])
    @test isapprox(tensor(QubitSGate()), [1 0; 0 1im])
    @test isapprox(tensor(QubitSxGate()), 0.5*[1+1im 1-1im; 1-1im 1+1im])
    @test isapprox(tensor(QubitHadamardGate()), sqrt(1/2)*[1 1; 1 -1])

    # Two-qubit gates 
    @test isapprox(
        tensor(QubitCNOTGate()),
        tensorproduct([1 0; 0 0], [0 1; 1 0]; tocache = false) +
        tensorproduct([0 0; 0 1], [1 0; 0 1]; tocache = false),
    )
    @test isapprox(
        tensor(QubitCNOTReverseGate()),
        tensorproduct([0 1; 1 0], [1 0; 0 0]; tocache = false) +
        tensorproduct([1 0; 0 1], [0 0; 0 1]; tocache = false),
    )
    @test isapprox(
        tensor(QubitCZGate()),
        tensorproduct([1 0; 0 0], [1 0; 0 -1]; tocache = false) +
        tensorproduct([0 0; 0 1], [1 0; 0 1]; tocache = false),
    )
    swap_gate = tensorproduct(ComplexF64[1 0; 0 1], [1 0; 0 1]; tocache = false)
    swap_gate = permutedims(swap_gate, (1, 4, 3, 2))
    @test isapprox(tensor(QubitSWAPGate()), swap_gate)
    swap_gate[1, 2, 2, 1] = swap_gate[2, 1, 1, 2] = 1im
    @test isapprox(tensor(QubitiSWAPGate()), swap_gate)

    # Single qubit rotation gates
    Rx = QubitRxGate(π)
    @test isapprox(tensor(Rx), [cos(π/2) -1im*sin(π/2); -1im*sin(π/2) cos(π/2)])
    setparams!(Rx, π/2)
    @test isapprox(tensor(Rx), [cos(π/4) -1im*sin(π/4); -1im*sin(π/4) cos(π/4)])

    Ry = QubitRyGate(π)
    @test isapprox(tensor(Ry), [cos(π/2) -sin(π/2); sin(π/2) cos(π/2)])
    setparams!(Ry, π/2)
    @test isapprox(tensor(Ry), [cos(π/4) -sin(π/4); sin(π/4) cos(π/4)])

    Rz = QubitRzGate(π)
    @test isapprox(tensor(Rz), [exp(-1im*π/2) 0; 0 exp(1im*π/2)])
    setparams!(Rz, π/2)
    @test isapprox(tensor(Rz), [exp(-1im*π/4) 0; 0 exp(1im*π/4)])

    phase = QubitPhaseGate(π)
    @test isapprox(tensor(phase), [1 0; 0 exp(1im*π)])
    setparams!(phase, π/2)
    @test isapprox(tensor(phase), [1 0; 0 exp(1im*π/2)])

    rot = QubitRotationGate(π, π/2, π/4)
    @test isapprox(
        tensor(rot),
        [cos(π/2) -exp(1im*π/2)*sin(π/2); exp(1im*π/4)*sin(π/2) exp(1im*(π/2+π/4))cos(π/2)],
    )
    setparams!(rot, (π/2, π/4, π))
    @test isapprox(
        tensor(rot),
        [cos(π/4) -exp(1im*π/4)*sin(π/4); exp(1im*π)*sin(π/4) exp(1im*(π/4+π))cos(π/4)],
    )
end
