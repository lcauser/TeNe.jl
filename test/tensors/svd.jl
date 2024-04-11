@testset "svd" begin
    A = rand(ComplexF64, 16, 20, 24)
    U, S, V = tsvd(A, 1)
    B = reshape(permutedims(A, (2, 3, 1)), (20*24, 16))
    t = LinearAlgebra.svd(B)
    @test isapprox(U, reshape(t.U, 20, 24, 16))
    @test isapprox(S, Diagonal(t.S))
    @test isapprox(V, t.Vt)

    A = rand(ComplexF64, 16, 20, 24)
    U, S, V = tsvd(A, 2)
    B = reshape(permutedims(A, (1, 3, 2)), (16*24, 20))
    t = LinearAlgebra.svd(B)
    @test isapprox((U, S, V), (reshape(t.U, 16, 24, 20), Diagonal(t.S), t.Vt))

    A = rand(ComplexF64, 16, 20, 24)
    U, S, V = tsvd(A, -1)
    B = reshape(A, (16*20, 24))
    t = LinearAlgebra.svd(B)
    @test isapprox((U, S, V), (reshape(t.U, 16, 20, 24), Diagonal(t.S), t.Vt))

    A = rand(ComplexF64, 16, 20, 24)
    U, S, V = tsvd(A, (1, 3))
    B = reshape(permutedims(A, (2, 1, 3)), (20, 16*24))
    t = LinearAlgebra.svd(B)
    @test isapprox((U, S, V), (t.U, Diagonal(t.S), reshape(t.Vt, 20, 16, 24)))
end