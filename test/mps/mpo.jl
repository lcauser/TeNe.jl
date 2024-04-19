@testset "mpo" begin 
    @test begin
        ψ = productmps(20, [0, 1])
        ϕ = productmps(20, [1, 0])
        O = productmpo(20, [0 1; 0 0])
        isapprox(inner(ϕ, O, ψ), 1)
    end

    @test begin
        ψ = productmps(20, [0, 1])
        ϕ = productmps(20, [1, 0])
        O = productmpo(20, [0 1; 0 0])
        isapprox(inner(ψ, O, ϕ), 0)
    end

    @test begin
        ψ = productmps(20, [0, 1])
        ϕ = productmps(20, [1, 0])
        O = productmpo(20, [0 1; 0 0])
        isapprox(inner(ϕ, adjoint(O), ψ), 0)
    end

    @test begin
        O1 = randommpo(30, 3, 4)
        O2 = randommpo(30, 3, 5)
        O3 = randommpo(30, 3, 6)
        trace(conj(O1), adjoint(O2), transpose(O3))
        true
    end

    @test begin
        ψ = randommps(5, 11, 4)
        O = randommpo(5, 11, 7)
        ψ′ = applympo(O, ψ; alg=:zipup)
        true
    end

    @test begin
        ψ = randommps(5, 11, 4)
        O = randommpo(5, 11, 7)
        ψ′ = applympo(O, ψ; alg=:naive)
        true
    end

    @test begin
        ψ = productmps(20, [0, 1])
        ϕ = productmps(20, [1, 0])
        O = productmpo(20, [0 1; 0 0])
        ψ′ = applympo(O, ψ; alg=:naive)
        isapprox(inner(ψ′, ϕ), 1)
    end

    @test begin
        ψ = productmps(20, [0, 1])
        ϕ = productmps(20, [1, 0])
        O = productmpo(20, [0 1; 0 0])
        ψ′ = applympo(O, ψ; alg=:zipup)
        isapprox(inner(ψ′, ϕ), 1)
    end

end
