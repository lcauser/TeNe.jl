@testset "statevector" begin 
    @test begin
        ψ = randomsv(3, 6)
        isapprox(ψ*ψ, 1)
    end

    @test begin
        ψ = productsv(8, [sqrt(0.2), sqrt(0.8)])
        ϕ = productsv(8, [1, 0])
        isapprox(ψ*ϕ, sqrt(0.2)^8)
    end
end