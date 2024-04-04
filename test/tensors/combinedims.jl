@testset "combinedims" begin
    A = randn(ComplexF64, 4, 6, 3, 2)
    B = zeros(ComplexF64, 3, 2, 24)
    combinedims!(B, A, (1, 2))
    C = reshape(permutedims(A, (3, 4, 1, 2)), (3, 2, 24))
    @test isapprox(B, C)

    B, cmb = combinedims(A, (1, 2))
    @test isapprox(B, C)

    B, cmb = combinedims(A, (3, 2))
    C = reshape(permutedims(A, (1, 4, 3, 2)), (4, 2, 18))
    @test isapprox(B, C)

    B, cmb = combinedims(A, (3, 4); return_copy=true)
    C = reshape(A, (4, 6, 6))
    @test isapprox(B, C)

    B, cmb = combinedims(A, (3, 4); return_copy=false)
    C = reshape(A, (4, 6, 6))
    @test isapprox(B, C)
end