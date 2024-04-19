@testset "gmps" begin 
    @test begin
        ψ = randomgmps(1, 2, 30, 16)
        movecenter!(ψ, 30)
        movecenter!(ψ, 1)
        true
    end

    @test begin
        ψ = randomgmps(2, 2, 30, 16)
        truncate!(ψ; maxdim=8)
        maxbonddim(ψ) == 8
    end

    @test begin
        ψ = randomgmps(3, 2, 30, 16)
        expand!(ψ, 32, 0.0)
        maxbonddim(ψ) == 32
    end

    @test begin 
        ψ = randomgmps(2, 2, 30, 16)
        ψconj = conj(ψ)
        Base.mightalias(ψconj[10], ψ[10])
    end

end
