@testset "creatempo" begin
    lt = Qubits()
    H = OpList(lt, 20)
    for i = 1:19
        add!(H, ["z", "z"], [i, i+1], 1)
    end
    for i = 1:20
        add!(H, ["x"], [i], 1)
    end
    H = MPO(H)
    @test maxbonddim(H) == 3

    ψ = productmps(lt, ["up" for _ = 1:20])
    @test isapprox(inner(ψ, H, ψ), 19)
    ψ = productmps(lt, ["s" for _ = 1:20])
    @test isapprox(inner(ψ, H, ψ), 20)
    ψ = productmps(lt, ["as" for _ = 1:20])
    @test isapprox(inner(ψ, H, ψ), -20)
end