@testset "mps-applympo" begin
    @test begin
        ψ = randommps(5, 11, 4)
        O = randommpo(5, 11, 7)
        ψ′ = applympo(O, ψ; alg = :zipup)
        true
    end

    @test begin
        ψ = randommps(5, 11, 4)
        O = randommpo(5, 11, 7)
        ψ′ = applympo(O, ψ; alg = :naive)
        true
    end

    @test begin
        ψ = productmps(20, [0, 1])
        ϕ = productmps(20, [1, 0])
        O = productmpo(20, [0 1; 0 0])
        ψ′ = applympo(O, ψ; alg = :naive)
        isapprox(inner(ψ′, ϕ), 1)
    end

    @test begin
        ψ = productmps(20, [0, 1])
        ϕ = productmps(20, [1, 0])
        O = productmpo(20, [0 1; 0 0])
        ψ′ = applympo(O, ψ; alg = :zipup)
        isapprox(inner(ψ′, ϕ), 1)
    end

    @test begin
        ψ = productsv(8, [0, 1])
        ϕ = productsv(8, [1, 0])
        O = productmpo(8, [0 1; 0 0])
        ψ′ = applympo(O, ψ)
        isapprox(inner(ψ′, ϕ), 1)
    end

    @test begin
        lt = Qubits()
        O = productmpo(lt, ["x" for _ = 1:20])
        ψ = productmps(lt, ["up" for _ = 1:20])
        ϕ = productmps(lt, ["dn" for _ = 1:20])
        ψ′ = applympo(O, ψ)
        isapprox(inner(ψ′, ϕ), 1)
    end
end
