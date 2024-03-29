@testset "tensorproduct" begin
    ### Test matrix-vector multiplication
    A = randn(ComplexF64, 5, 6)
    B = randn(ComplexF64, 3, 4)
    C = tensorproduct(A, B, true, false)
    @test size(C) == (size(A)..., size(B)...)
end