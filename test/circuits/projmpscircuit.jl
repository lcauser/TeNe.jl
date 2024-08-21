@testset "projmpscircuit" begin
    @testset "mps-connector" begin
        # The first test will create a state and apply a circuit.
        # It will test that the projection creates the correct MPSs.
        ψ = randommps(2, 10, 1)
        U = randombwcircuit(2, 10, 5)
        ϕ = randommps(2, 10, 1)
        ψ′ = applympo(U, ψ)
        ol1 = inner(ϕ, ψ′)
        projU = ProjMPSCircuit(ϕ, U, ψ)
        ol2s = []
        for i = 2:5
            movecenter!(projU, i)
            ol2 = inner(TeNe.topblock(projU, i-1), TeNe.bottomblock(projU, i))
            push!(ol2s, ol2)
        end
        println(ol1)
        println(ol2s)
    end
    
end