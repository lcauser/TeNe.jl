@testset "mps-trace" begin 
    @test begin
        O1 = randommpo(30, 3, 4)
        O2 = randommpo(30, 3, 5)
        O3 = randommpo(30, 3, 6)
        trace(conj(O1), adjoint(O2), transpose(O3))
        true
    end
end