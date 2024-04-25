@testset "stateoperator" begin 
    @test begin
        ψ = productsv(8, [0, 1])
        ϕ = productsv(8, [1, 0])
        O = productso(8, [0 1; 0 0])
        isapprox(inner(ϕ, O, ψ), 1)
    end

    @test begin
        ψ = productsv(8, [0, 1])
        ϕ = productsv(8, [1, 0])
        O = productso(8, [0 1; 0 0])
        isapprox(inner(ψ, O, ϕ), 0)
    end

    @test begin
        ψ = productsv(8, [1, 0])
        O = productso(8, [0 1; 0 0])
        isapprox(inner(ψ, O, adjoint(O), ψ), 1)
    end
    @test begin
        ψ = productsv(8, [1, 0])
        ϕ = productsv(8, [0, 1])
        O = productso(8, [0 1; 0 0])
        isapprox(inner(ψ, O, adjoint(O), O, ϕ), 1)
    end

    @test begin
        ψ = productsv(8, [0, 1])
        ϕ = productsv(8, [1, 0])
        O = productso(8, [0 1; 0 0])
        isapprox(inner(ϕ, adjoint(O), ψ), 0)
    end

    @test begin
        O1 = randomso(8, 3)
        O2 = randomso(8, 3)
        O3 = randomso(8, 3)
        trace(conj(O1), adjoint(O2), transpose(O3))
        true
    end

    @test begin 
        ψ = productsv(8, [0, 1])
        ϕ = productsv(8, [1, 0])
        O = productso(8, [0 1; 0 0])
        ψ2 = O * ψ
        isapprox(ϕ.tensor, ψ2.tensor)
    end

    @test begin 
        O1 = productso(8, [0 1; 0 0])
        O2 = productso(8, [0 1; 0 0] * adjoint([0 1; 0 0]))
        O3 = O1 * adjoint(O1)
        isapprox(O2.tensor, O3.tensor)
    end

    @test begin 
        O1 = productso(8, [0 1; 0 0])
        O2 = productso(8, [0 0; 1 0])
        O3 = productso(8, [0 1; 0 0] * adjoint([0 1; 0 0]))
        O4 = O1 * O2
        isapprox(O3.tensor, O4.tensor)
    end
end