@testset "exp" begin
    A = rand(ComplexF64, 4, 3, 3, 4, 3, 3)
    B = TeNe.exp(A, (1, 5, 6), prefactor=-1im)
    C = -1im * permutedims(A, (2, 3, 4, 1, 5, 6))
    C = LinearAlgebra.exp(reshape(C, (36, 36)))
    C = reshape(C, (3, 3, 4, 4, 3, 3))
    C = permutedims(C, (4, 1, 2, 3, 5, 6))
    @test isapprox(B, C)

    B = TeNe.exp(A, (1, 2, 3), (4, 5, 6))
    C = reshape(LinearAlgebra.exp(reshape(A, (36, 36))), (4, 3, 3, 4, 3, 3))
    @test isapprox(B, C)
end