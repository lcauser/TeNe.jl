@testset "circuit" begin
    @testset "mps-connector" begin
        circuit = CircuitMPS(2, 20)
        for i = 1:10
            add!(circuit, TeNe._unitary_close_to_id(2, 2, 0.5), (2*i-1, 2*i))
        end

        @test length(circuit.layers) == 1
        add!(circuit, TeNe._unitary_close_to_id(2, 2, 0.5), (4, 5))
        @test length(circuit.layers) == 2
        add!(circuit, TeNe._unitary_close_to_id(2, 2, 0.5), (4, 5))
        @test length(circuit.layers) == 3

        circuit = randombwcircuit(2, 20, 3)
        @test length(circuit.layers) == 3
        circuit = randombwcircuit(2, 20, 4)
        @test length(circuit.layers) == 4
    end

end
