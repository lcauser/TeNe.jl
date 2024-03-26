using SimpleTensors
using Test

@testset "contract.jl" begin
    ### Test matrix-vector multiplication
    A = randn(ComplexF64, 5, 5)
    B = randn(ComplexF64, 4)
    try C = contract(A, B, 1, 2) catch e  @test (e isa ArgumentError) end
    try C = contract(A, B, 1, 1) catch e  @test (e isa ArgumentError) end

    A = randn(ComplexF64, 5, 4)
    check = false
    try
        C = contract(A, B, 2, 1)
        check = true
    catch e end
    @test check

    A = randn(ComplexF64, 4, 4)
    @test contract(A, B, 2, 1) == A * B
    @test contract(A, B, 1, 1) == transpose(A) * B
    @test contract(A, B, 1, 1, true) == adjoint(A) * B
    @test contract(A, B, 1, 1, true, true) == adjoint(A) * conj(B)


    ### Test matrix-matrix multiplications
    A = randn(ComplexF64, 5, 5)
    B = randn(ComplexF64, 4, 4)
    try C = contract(A, B, 1, 3) catch e  @test (e isa ArgumentError) end
    try C = contract(A, B, 1, 1) catch e  @test (e isa ArgumentError) end

    B = randn(ComplexF64, 5, 4)
    check = false
    try
        C = contract(A, B, 2, 1)
        check = true
    catch e end
    @test check

    A = randn(ComplexF64, 3, 4)
    B = randn(ComplexF64, 4, 5)
    @test contract(A, B, 2, 1) == A * B
    @test contract(A, B, 2, 1, true) == conj(A) * B
    @test contract(A, B, 2, 1, true, true) == conj(A) * conj(B)
    B = randn(ComplexF64, 5, 4)
    @test contract(A, B, 2, 2) == A * transpose(B)
    B = randn(ComplexF64, 4, 3)
    @test contract(A, B, 1, 2) == transpose(A) * transpose(B)

    # General tensor contractions
    A = rand(ComplexF64, 4, 5, 6)
    B = rand(ComplexF64, 3, 6, 5)
    @test size(contract(A, B, 2, 3)) == (4, 6, 3, 6)
    @test size(contract(A, B, (2, 3), (3, 2))) == (4, 3)
    try C = contract(A, B, 1, 1) catch e  @test (e isa ArgumentError) end
    try C = contract(A, B, [1, 3], [1, 2]) catch e  @test (e isa ArgumentError) end

    A = rand(ComplexF64, 2, 3, 4, 5)
    B = rand(ComplexF64, 4, 5, 6, 7)
    C = reshape(reshape(A, (6, 20)) * reshape(B, (20, 42)), (2, 3, 6, 7))
    @test C == contract(A, B, (3, 4), (1, 2))
end