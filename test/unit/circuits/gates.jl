@testset "circuitgate" begin
    # Test the validation 
    @test_throws ArgumentError("The gates must have an even number of dimensions.") creategate(
        rand(4, 4, 4, 4, 4),
    )
    @test_throws ArgumentError("The gate must have the same dimensions sizes.") creategate(
        rand(4, 4, 4, 2),
    )

    # Test the polar decomposition
    U = creategate(randn(ComplexF64, 2, 2))
    makeunitary!(U)
    U2 = contract(tensor(U), tensor(U), 2, 2, false, true)
    @test isapprox(U2, [1 0; 0 1])

    U = creategate(randn(ComplexF64, 2, 2, 2, 2))
    makeunitary!(U)
    U2 = contract(tensor(U), tensor(U), (2, 4), (2, 4), false, true)
    U2 = permutedims(U2, (1, 3, 2, 4))
    @test isapprox(U2, tensorproduct([1 0; 0 1], [1 0; 0 1]))

    U = creategate(randn(ComplexF64, 2, 2, 2, 2, 2, 2))
    makeunitary!(U)
    U2 = contract(tensor(U), tensor(U), (2, 4, 6), (2, 4, 6), false, true)
    U2 = permutedims(U2, (1, 4, 2, 5, 3, 6))
    @test isapprox(U2, tensorproduct(tensorproduct([1 0; 0 1], [1 0; 0 1]), [1 0; 0 1]))
end


@testset "applygate-statevector" begin
    # Test the validation 
    ψ = randomsv(2, 10)
    U = makeunitary(creategate(randn(ComplexF64, 3, 3, 3, 3)))
    @test_throws DomainError("The length of sites does not match the length of the gate.") applygate!(
        ψ,
        U,
        (1, 2, 3),
    )
    @test_throws DomainError(
        "The list of sites (2, 11) does not fall between 1 and $(length(ψ)).",
    ) applygate!(ψ, U, (2, 11))
    @test_throws ArgumentError(
        "The $(typeof(ψ)) and $(typeof(U)) have incomptible dimensions.",
    ) applygate!(ψ, U, (2, 9))

    # Test application 
    ψ = randomsv(2, 10)
    U = makeunitary(creategate(randn(ComplexF64, 2, 2)))
    applygate!(U, ψ, 5)
    @test isapprox(ψ*ψ, 1)
    applygate!(ψ, U, 8)
    @test isapprox(ψ*ψ, 1)

    ψ = randomsv(2, 10)
    U = makeunitary(creategate(randn(ComplexF64, 2, 2, 2, 2)))
    applygate!(U, ψ, (3, 7))
    @test isapprox(ψ*ψ, 1)
    applygate!(ψ, U, (10, 2))
    @test isapprox(ψ*ψ, 1)

    ψ = randomsv(2, 10)
    U = makeunitary(creategate(randn(ComplexF64, 2, 2, 2, 2, 2, 2)))
    applygate!(U, ψ, (1, 5, 3))
    @test isapprox(ψ*ψ, 1)
    applygate!(ψ, U, (3, 9, 10))
    @test isapprox(ψ*ψ, 1)

    ψ = randomsv(2, 10)
    ψ′ = deepcopy(ψ)
    U = makeunitary(creategate(randn(ComplexF64, 2, 2, 2, 2)))
    applygate!(U, ψ, (4, 5))
    tensor(U) .= conj.(tensor(U))
    applygate!(ψ, U, (4, 5))
    @test isapprox(ψ*ψ′, 1)
end

@testset "applygate-stateoperator" begin
    # Test the validation 
    O = randomso(2, 10)
    U = makeunitary(creategate(randn(ComplexF64, 3, 3, 3, 3)))
    @test_throws DomainError("The length of sites does not match the length of the gate.") applygate!(
        O,
        U,
        (1, 2, 3),
    )
    @test_throws DomainError(
        "The list of sites (2, 11) does not fall between 1 and $(length(O)).",
    ) applygate!(O, U, (2, 11))
    @test_throws ArgumentError(
        "The $(typeof(O)) and $(typeof(U)) have incomptible dimensions.",
    ) applygate!(O, U, (2, 9))
    @test_throws DomainError("The length of sites does not match the length of the gate.") applygate!(
        U,
        O,
        (1, 2, 3),
    )
    @test_throws DomainError(
        "The list of sites (2, 11) does not fall between 1 and $(length(O)).",
    ) applygate!(U, O, (2, 11))
    @test_throws ArgumentError(
        "The $(typeof(O)) and $(typeof(U)) have incomptible dimensions.",
    ) applygate!(U, O, (2, 9))

    # Test the application
    O = randomso(2, 10)
    O′ = deepcopy(O)
    U = makeunitary(creategate(randn(ComplexF64, 2, 2)))
    applygate!(U, O, 5)
    tensor(U) .= adjoint(tensor(U))
    applygate!(U, O, 5)
    @test isapprox(tensor(O), tensor(O′))

    O = randomso(2, 10)
    O′ = deepcopy(O)
    U = makeunitary(creategate(randn(ComplexF64, 2, 2)))
    applygate!(O, U, 3)
    tensor(U) .= adjoint(tensor(U))
    applygate!(O, U, 3)
    @test isapprox(tensor(O), tensor(O′))

    O = randomso(2, 10)
    O′ = deepcopy(O)
    U = makeunitary(creategate(randn(ComplexF64, 2, 2, 2, 2)))
    applygate!(U, O, (3, 7))
    tensor(U) .= conj.(tensor(U))
    tensor(U) .= permutedims(tensor(U), (2, 1, 4, 3))
    applygate!(U, O, (3, 7))
    @test isapprox(tensor(O), tensor(O′))

    O = randomso(2, 10)
    O′ = deepcopy(O)
    U = makeunitary(creategate(randn(ComplexF64, 2, 2, 2, 2)))
    applygate!(O, U, (9, 1))
    tensor(U) .= conj.(tensor(U))
    tensor(U) .= permutedims(tensor(U), (2, 1, 4, 3))
    applygate!(O, U, (9, 1))
    @test isapprox(tensor(O), tensor(O′))
end

@testset "applygate-mps" begin
    # Test the validation 
    ψ = randommps(2, 20, 3)
    U = makeunitary(creategate(randn(ComplexF64, 3, 3, 3, 3)))
    @test_throws DomainError(
        "The list of sites (20, 21) does not fall between 1 and $(length(ψ)).",
    ) applygate!(U, ψ, 20)
    @test_throws ArgumentError(
        "The $(typeof(ψ)) and $(typeof(U)) have incomptible dimensions.",
    ) applygate!(U, ψ, 5)
    @test_throws DomainError(
        "The list of sites (20, 21) does not fall between 1 and $(length(ψ)).",
    ) applygate!(U, ψ, 21, true)
    @test_throws ArgumentError(
        "The $(typeof(ψ)) and $(typeof(U)) have incomptible dimensions.",
    ) applygate!(U, ψ, 5, true)

    # Test the application 
    ψ = randommps(2, 20, 3)
    U = makeunitary(creategate(randn(ComplexF64, 2, 2)))
    applygate!(U, ψ, 5)
    @test isapprox(ψ*ψ, 1)
    applygate!(ψ, U, 13)
    @test isapprox(ψ*ψ, 1)

    ψ = randommps(2, 20, 3)
    U = makeunitary(creategate(randn(ComplexF64, 2, 2, 2, 2)))
    applygate!(U, ψ, 3)
    @test isapprox(ψ*ψ, 1)
    applygate!(ψ, U, 10)
    @test isapprox(ψ*ψ, 1)

    ψ = randommps(2, 20, 3)
    U = makeunitary(creategate(randn(ComplexF64, 2, 2, 2, 2, 2, 2)))
    applygate!(U, ψ, 1)
    @test isapprox(ψ*ψ, 1)
    applygate!(ψ, U, 4)
    @test isapprox(ψ*ψ, 1)

    ψ = randommps(2, 20, 3)
    ψ′ = deepcopy(ψ)
    U = makeunitary(creategate(randn(ComplexF64, 2, 2, 2, 2)))
    applygate!(U, ψ, 7)
    tensor(U) .= conj.(tensor(U))
    applygate!(ψ, U, 7)
    @test isapprox(ψ*ψ′, 1)
end

@testset "gate-manipulations" begin
    # Test conjugate 
    U = makeunitary(creategate(randn(ComplexF64, 2, 2, 2, 2)))
    U2 = conj(U)
    @test isapprox(conj(U.gate), U2.gate)

    # Test transpose 
    U = makeunitary(creategate(randn(ComplexF64, 2, 2, 2, 2)))
    U2 = transpose(U)
    @test isapprox(permutedims(U.gate, (2, 1, 4, 3)), U2.gate)

    # Test adjoint 
    U = makeunitary(creategate(randn(ComplexF64, 2, 2, 2, 2)))
    U2 = adjoint(U)
    U3 = conj(transpose(U))
    @test isapprox(conj(permutedims(U.gate, (2, 1, 4, 3))), U2.gate)
    @test isapprox(U3.gate, U2.gate)
end
