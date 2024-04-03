@testset "combinedims" begin
    A = randn(ComplexF64, 4, 6, 3, 2)
    B = zeros(ComplexF64, 3, 2, 24)
    combine
end