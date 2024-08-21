@testset "projmpscircuit" begin
    @testset "mps-connector" begin
        # The first test will create a state and apply a circuit.
        # It will test that the projection creates the correct MPSs.
        ψ = randommps(2, 10, 1)
        U = randombwcircuit(2, 10, 5)
        ϕ = randommps(2, 10, 1)
        ϕ′ = deepcopy(ϕ)
        ψ′ = applygates(U, ψ)
        ol1 = inner(ϕ, ψ′)
        projU = ProjMPSCircuit(ϕ, U, ψ; cutoff=0.0)
        ol2s = []
        for i = 2:5
            movecenter!(projU, i)
            ol2 = inner(conj(TeNe.topblock(projU, i-1)), TeNe.bottomblock(projU, i))
            push!(ol2s, ol2)
        end
        @test all(isapprox.(ol2s, ol1))

        ## Test for building left blocks 
        ol2s = []
        for i = 1:5
            movecenter!(projU, i)
            left = ones(Float64, 1, 1)
            for j = 1:length(U.layers[end+1-i].sites)+1
                left = copy(TeNe._buildleft(projU, j, left))
            end
            push!(ol2s, left[])
        end
        @test all(isapprox.(ol2s, ol1))
    end
    
end