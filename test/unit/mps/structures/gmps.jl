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
        expand!(ψ, 32, 0.01)
        maxbonddim(ψ) == 32
    end

    @test begin 
        ψ = randomgmps(2, 2, 30, 16)
        ψconj = conj(ψ)
        Base.mightalias(ψconj[10], ψ[10])
    end

    # Tests for checking replacesites! functionality 
    @test begin
        ψ = randommps(2, 20, 9)
        ψ2 = deepcopy(ψ)
        site = 7
        range = 3
        movecenter!(ψ, site)
        ten = ψ[site]
        for i = 1:range-1
            ten = contract(ten, ψ[site+i], ndims(ten), 1)
        end
        TeNe.replacesites!(ψ, ten, site; cutoff=1e-16)
        isapprox(inner(ψ, ψ2), 1.0)
    end

    @test begin
        ψ = randommps(2, 20, 9)
        ψ2 = deepcopy(ψ)
        site = 14
        range = 4
        movecenter!(ψ, site)
        ten = ψ[site-3]
        for i = 1:range-1
            ten = contract(ten, ψ[site-3+i], ndims(ten), 1)
        end
        TeNe.replacesites!(ψ, ten, site, true; cutoff=1e-16)
        isapprox(inner(ψ, ψ2), 1.0)
    end

    @test begin
        ψ = randommpo(2, 30, 13)
        ψ2 = deepcopy(ψ)
        site = 11
        range = 3
        movecenter!(ψ, site)
        ten = ψ[site]
        for i = 1:range-1
            ten = contract(ten, ψ[site+i], ndims(ten), 1)
        end
        TeNe.replacesites!(ψ, ten, site; cutoff=1e-16)
        isapprox(trace(ψ), trace(ψ2))
    end

    @test begin
        ψ = randommpo(2, 35, 4)
        ψ2 = deepcopy(ψ)
        site = 18
        range = 4
        movecenter!(ψ, site)
        ten = ψ[site-3]
        for i = 1:range-1
            ten = contract(ten, ψ[site-3+i], ndims(ten), 1)
        end
        TeNe.replacesites!(ψ, ten, site, true; cutoff=1e-16)
        isapprox(trace(ψ), trace(ψ2))
    end

end
