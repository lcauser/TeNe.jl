@testset "latticesites-qubits" begin
    lt = Qubits()
    @test isapprox(state(lt, "up"), [1, 0])
    @test isapprox(op(lt, "x"), [0 1; 1 0])
end
