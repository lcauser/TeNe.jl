@testset "mps-inner" begin 
    @test begin
        ψ = randommps(3, 21, 7)
        isapprox(ψ*ψ, 1)
    end

    @test begin
        ψ = productmps(16, [sqrt(0.2), sqrt(0.8)])
        ϕ = productmps(16, [1, 0])
        isapprox(ψ*ϕ, sqrt(0.2)^16)
    end

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
        lt = Qubits()
        O = productmpo(lt, ["x" for _ = 1:20])
        ψ = productmps(lt, ["up" for _ = 1:20])
        ϕ = productmps(lt, ["dn" for _ = 1:20])
        isapprox(inner(ψ, O, ϕ), 1)
    end
end