@testset "trace" begin
    ### Trace over one
    A = randn(ComplexF64, 5, 6, 7)
    B = trace(A, 3)
    @test size(B) == (5, 6)
    C = zeros(ComplexF64, 5, 6)
    for i = 1:5
        for j = 1:6
            for k = 1:7
                C[i, j] += A[i, j, k]
            end
        end
    end
    @test C == B

    ### Trace over two 
    A = randn(ComplexF64, 5, 6, 7, 6)
    try
        B = trace(A, 1, 2; conj = true)
    catch e
        @test (e isa ArgumentError)
    end
    try
        B = trace(A, 1, 2, 3)
    catch e
        @test (e isa ArgumentError)
    end
    B = trace(A, 2, 4)
    @test size(B) == (5, 7)
    C = zeros(ComplexF64, 5, 7)
    for i = 1:5
        for j = 1:7
            for k = 1:6
                C[i, j] += A[i, k, j, k]
            end
        end
    end
    @test C == B

    ### Trace over 3
    A = randn(ComplexF64, 3, 3, 4, 5, 3, 6)
    try
        B = trace(A, 1, 3; conj = true)
    catch e
        @test (e isa ArgumentError)
    end
    try
        B = trace(A, 1, 2, 4)
    catch e
        @test (e isa ArgumentError)
    end
    B = trace(A, 1, 2, 5; conj = true)
    @test size(B) == (4, 5, 6)
    C = zeros(ComplexF64, 4, 5, 6)
    for i = 1:4
        for j = 1:5
            for k = 1:6
                for n = 1:3
                    C[i, j, k] += conj(A[n, n, i, j, n, k])
                end
            end
        end
    end
    @test C == B

end
