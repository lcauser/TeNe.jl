@testset "ProjMPS" begin
    @testset "ψ ⋅ ψ" begin 
        ψ = randommps(3, 35, 6)
        inner1 = inner(ψ, ψ)
        projψ = ProjMPS(ψ, ψ; sym=true)
        movecenter!(projψ, 9)
        @test isapprox(inner1, inner(projψ))
        A = ψ[9]
        @test isapprox(inner1, inner(projψ, A))
        @test isapprox(inner1, inner(projψ, A, true))
        A = contract(A, ψ[10], 3, 1; tocache=false)
        @test isapprox(inner1, inner(projψ, A))
        movecenter!(projψ, 10)
        @test isapprox(inner1, inner(projψ, A, true))
        movecenter!(projψ, 9)
        A = contract(A, ψ[11], 4, 1; tocache=false)
        @test isapprox(inner1, inner(projψ, A))
        movecenter!(projψ, 11)
        @test isapprox(inner1, inner(projψ, A, true))
    end

    @testset "ψ ⋅ O ⋅ ψ" begin 
        ψ = randommps(2, 11, 6)
        O = randommpo(2, 11, 5)
        inner1 = inner(ψ, O, ψ)
        projψ = ProjMPS(ψ, O, ψ; sym=true)
        movecenter!(projψ, 9)
        @test isapprox(inner1, inner(projψ))
        A = ψ[9]
        @test isapprox(inner1, inner(projψ, A))
        @test isapprox(inner1, inner(projψ, A, true))
        A = contract(A, ψ[10], 3, 1; tocache=false)
        @test isapprox(inner1, inner(projψ, A))
        movecenter!(projψ, 10)
        @test isapprox(inner1, inner(projψ, A, true))
        movecenter!(projψ, 9)
        A = contract(A, ψ[11], 4, 1; tocache=false)
        @test isapprox(inner1, inner(projψ, A))
        movecenter!(projψ, 11)
        @test isapprox(inner1, inner(projψ, A, true))
    end

    @testset "ψ ⋅ O ⋅ P ⋅ ψ" begin 
        ψ = randommps(2, 15, 6)
        O = randommpo(2, 15, 5)
        P = randommpo(2, 15, 5)
        inner1 = inner(ψ, O, P, ψ)
        projψ = ProjMPS(ψ, O, P, ψ; sym=true)
        movecenter!(projψ, 9)
        @test isapprox(inner1, inner(projψ))
        A = ψ[9]
        @test isapprox(inner1, inner(projψ, A))
        @test isapprox(inner1, inner(projψ, A, true))
        A = contract(A, ψ[10], 3, 1; tocache=false)
        @test isapprox(inner1, inner(projψ, A))
        movecenter!(projψ, 10)
        @test isapprox(inner1, inner(projψ, A, true))
        movecenter!(projψ, 9)
        A = contract(A, ψ[11], 4, 1; tocache=false)
        @test isapprox(inner1, inner(projψ, A))
        movecenter!(projψ, 11)
        @test isapprox(inner1, inner(projψ, A, true))
    end

    @testset "ϕ . ψ" begin 
        ϕ = randommps(2, 31, 3)
        ψ = randommps(2, 31, 3)
        inner1 = inner(ϕ, ψ)
        projψ = ProjMPS(ϕ, ψ; sym=false)
        movecenter!(projψ, 9)
        @test isapprox(inner1, inner(projψ))
        A = ψ[9]
        @test isapprox(inner1, inner(projψ, A))
        @test isapprox(inner1, inner(projψ, A, true))
        A = contract(A, ψ[10], 3, 1; tocache=false)
        @test isapprox(inner1, inner(projψ, A))
        movecenter!(projψ, 10)
        @test isapprox(inner1, inner(projψ, A, true))
        movecenter!(projψ, 9)
        A = contract(A, ψ[11], 4, 1; tocache=false)
        @test isapprox(inner1, inner(projψ, A))
        movecenter!(projψ, 11)
        @test isapprox(inner1, inner(projψ, A, true))
    end
end