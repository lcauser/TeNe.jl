@testset "stateoperator" begin
    lt = Qubits()
    H = OpList(lt, 10)
    for i = 1:9
        add!(H, ["z", "z"], [i, i+1])
    end
    for i = 1:10
        add!(H, "x", i, 0.5)
    end
    H = StateOperator(H)
    ψ = productsv(lt, ["up" for _ = 1:10])
    @test isapprox(inner(ψ, H, ψ), 9.0)
    ψ = productsv(lt, ["s" for _ = 1:10])
    @test isapprox(inner(ψ, H, ψ), 5.0)
end