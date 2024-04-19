@testset "mps" begin 
    @test begin
        ψ = randommps(3, 21, 7)
        isapprox(ψ*ψ, 1)
    end

    @test begin
        ψ = productmps(16, [sqrt(0.2), sqrt(0.8)])
        ϕ = productmps(16, [1, 0])
        isapprox(ψ*ϕ, sqrt(0.2)^16)
    end
end
