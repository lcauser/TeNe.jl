@testset "latticesites-bosons" begin 
    for dim = 2:10
        lt = Bosons(dim)
        for i = Base.OneTo(dim)
            state1 = state(lt, string(i-1))
            state2 = zeros(eltype(state1), dim)
            state2[i] = 1
            @test isapprox(state1, state2)
        end

        @test isapprox(op(lt, "n"), op(lt, "adag")*op(lt, "a"))
    end
end