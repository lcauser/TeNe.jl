@testset "statetensor-trace" begin
    @test begin
        O1 = randomso(8, 3)
        O2 = randomso(8, 3)
        O3 = randomso(8, 3)
        trace(conj(O1), adjoint(O2), transpose(O3))
        true
    end
end
