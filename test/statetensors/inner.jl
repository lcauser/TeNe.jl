@testset "statetensor-inner" begin 
    @test begin
        ψ = randomsv(3, 6)
        isapprox(ψ*ψ, 1)
    end

    @test begin
        ψ = productsv(8, [sqrt(0.2), sqrt(0.8)])
        ϕ = productsv(8, [1, 0])
        isapprox(ψ*ϕ, sqrt(0.2)^8)
    end

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
end