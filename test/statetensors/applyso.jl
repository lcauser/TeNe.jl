@testset "statetensor-applyso" begin 
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