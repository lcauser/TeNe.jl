@testset "tensorproduct" begin
    # Vectors 
    A = randn(ComplexF64, 5)
    B = randn(ComplexF64, 3)
    C = tensorproduct(A, B, false, false)
    @test size(C) == (size(A)..., size(B)...)
    D = zeros(ComplexF64, 5, 3)
    for i = 1:5
        for j = 1:3
            D[i, j] = A[i] * B[j]
        end
    end
    @test isapprox(C, D)
    E = zeros(ComplexF64, 5, 3)
    tensorproduct!(E, A, B)
    @test isapprox(E, D)

    # Matrices
    A = randn(ComplexF64, 5, 6)
    B = randn(ComplexF64, 3, 4)
    C = tensorproduct(A, B, true, false)
    @test size(C) == (size(A)..., size(B)...)
    D = zeros(ComplexF64, 5, 6, 3, 4)
    for i = 1:5
        for j = 1:6
            for k = 1:3
                for m = 1:4
                    D[i, j, k, m] = conj(A[i, j]) * B[k, m]
                end
            end
        end
    end
    @test isapprox(C, D)
end
