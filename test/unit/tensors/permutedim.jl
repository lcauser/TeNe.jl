@testset "permutedim" begin
    A = randn(ComplexF64, 4, 6, 3, 2)
    try
        B = permutedim(A, 2, 5)
    catch e
        @test (e isa ArgumentError)
    end
    try
        B = permutedim(A, 2, -2)
    catch e
        @test (e isa ArgumentError)
    end

    B = permutedim(A, 1, 3)
    C = permutedims(A, (2, 3, 1, 4))
    @test isapprox(B, C)

    B = permutedim(A, 1, -1)
    C = permutedims(A, (2, 3, 4, 1))
    @test isapprox(B, C)

    B = zeros(ComplexF64, 3, 4, 6, 2)
    permutedim!(B, A, 3, 1)
    C = permutedims(A, (3, 1, 2, 4))
    @test isapprox(B, C)
end
