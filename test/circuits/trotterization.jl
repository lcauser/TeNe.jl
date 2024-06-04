@testset "Trotterization" begin
    # Testing _create_trotter_gate_compressed
    @testset "_create_trotter_gate_compressed" begin
        @testset "two sites" begin
            H = OpList(Qubits(), 20)
            for i = 1:20
                add!(H, "x", i, 0.34)
            end
            for i = 1:19
                add!(H, ["z", "z"], [i, i+1], 0.56)
            end

            @test begin
                op = TeNe._create_trotter_gate_compressed(H, 0.1, 5)
                op2 = 0.56*tensorproduct([1 0; 0 -1], [1 0; 0 -1]; tocache=false)
                op2 += 0.5*0.34*tensorproduct([0 1; 1 0], [1 0; 0 1]; tocache=false)
                op2 += 0.5*0.34*tensorproduct([1 0; 0 1], [0 1; 1 0]; tocache=false)
                op2 = TeNe.exp(0.1*op2, (2, 4))
                isapprox(op, op2)
            end

            @test begin
                op = TeNe._create_trotter_gate_compressed(H, 0.1, 1)
                op2 = 0.56*tensorproduct([1 0; 0 -1], [1 0; 0 -1]; tocache=false)
                op2 += 0.34*tensorproduct([0 1; 1 0], [1 0; 0 1]; tocache=false)
                op2 += 0.5*0.34*tensorproduct([1 0; 0 1], [0 1; 1 0]; tocache=false)
                op2 = TeNe.exp(0.1*op2, (2, 4))
                isapprox(op, op2)
            end

            @test begin
                op = TeNe._create_trotter_gate_compressed(H, 0.1, 19)
                op2 = 0.56*tensorproduct([1 0; 0 -1], [1 0; 0 -1]; tocache=false)
                op2 += 0.5*0.34*tensorproduct([0 1; 1 0], [1 0; 0 1]; tocache=false)
                op2 += 0.34*tensorproduct([1 0; 0 1], [0 1; 1 0]; tocache=false)
                op2 = TeNe.exp(0.1*op2, (2, 4))
                isapprox(op, op2)
            end
        end

        @testset "three sites" begin
            H = OpList(Qubits(), 20)
            for i = 1:20
                add!(H, "x", i, 0.34)
            end
            for i = 1:19
                add!(H, ["z", "z"], [i, i+1], 0.56)
            end
            for i = 1:18
                add!(H, ["z", "x", "z"], [i, i+1, i+2], -1.4)
            end

            @test begin
                op = TeNe._create_trotter_gate_compressed(H, 0.1, 5)
                op2 = -1.4*tensorproduct(tensorproduct([1 0; 0 -1], [0 1; 1 0]; tocache=false), [1 0; 0 -1]; tocache=false)
                op2 += 0.56/2*tensorproduct(tensorproduct([1 0; 0 -1], [1 0; 0 -1]; tocache=false), [1 0; 0 1]; tocache=false)
                op2 += 0.56/2*tensorproduct(tensorproduct([1 0; 0 1], [1 0; 0 -1]; tocache=false), [1 0; 0 -1]; tocache=false)
                op2 += 0.34/3*tensorproduct(tensorproduct([0 1; 1 0], [1 0; 0 1]; tocache=false), [1 0; 0 1]; tocache=false)
                op2 += 0.34/3*tensorproduct(tensorproduct([1 0; 0 1], [0 1; 1 0]; tocache=false), [1 0; 0 1]; tocache=false)
                op2 += 0.34/3*tensorproduct(tensorproduct([1 0; 0 1], [1 0; 0 1]; tocache=false), [0 1; 1 0]; tocache=false)
                op2 = TeNe.exp(0.1*op2, (2, 4, 6))
                isapprox(op, op2)
            end

            @test begin
                op = TeNe._create_trotter_gate_compressed(H, 0.1, 1)
                op2 = -1.4*tensorproduct(tensorproduct([1 0; 0 -1], [0 1; 1 0]; tocache=false), [1 0; 0 -1]; tocache=false)
                op2 += 0.56*tensorproduct(tensorproduct([1 0; 0 -1], [1 0; 0 -1]; tocache=false), [1 0; 0 1]; tocache=false)
                op2 += 0.56/2*tensorproduct(tensorproduct([1 0; 0 1], [1 0; 0 -1]; tocache=false), [1 0; 0 -1]; tocache=false)
                op2 += 0.34*tensorproduct(tensorproduct([0 1; 1 0], [1 0; 0 1]; tocache=false), [1 0; 0 1]; tocache=false)
                op2 += 0.34/2*tensorproduct(tensorproduct([1 0; 0 1], [0 1; 1 0]; tocache=false), [1 0; 0 1]; tocache=false)
                op2 += 0.34/3*tensorproduct(tensorproduct([1 0; 0 1], [1 0; 0 1]; tocache=false), [0 1; 1 0]; tocache=false)
                op2 = TeNe.exp(0.1*op2, (2, 4, 6))
                isapprox(op, op2)
            end

            @test begin
                op = TeNe._create_trotter_gate_compressed(H, 0.1, 2)
                op2 = -1.4*tensorproduct(tensorproduct([1 0; 0 -1], [0 1; 1 0]; tocache=false), [1 0; 0 -1]; tocache=false)
                op2 += 0.56/2*tensorproduct(tensorproduct([1 0; 0 -1], [1 0; 0 -1]; tocache=false), [1 0; 0 1]; tocache=false)
                op2 += 0.56/2*tensorproduct(tensorproduct([1 0; 0 1], [1 0; 0 -1]; tocache=false), [1 0; 0 -1]; tocache=false)
                op2 += 0.34/2*tensorproduct(tensorproduct([0 1; 1 0], [1 0; 0 1]; tocache=false), [1 0; 0 1]; tocache=false)
                op2 += 0.34/3*tensorproduct(tensorproduct([1 0; 0 1], [0 1; 1 0]; tocache=false), [1 0; 0 1]; tocache=false)
                op2 += 0.34/3*tensorproduct(tensorproduct([1 0; 0 1], [1 0; 0 1]; tocache=false), [0 1; 1 0]; tocache=false)
                op2 = TeNe.exp(0.1*op2, (2, 4, 6))
                isapprox(op, op2)
            end

            @test begin
                op = TeNe._create_trotter_gate_compressed(H, 0.1, 18)
                op2 = -1.4*tensorproduct(tensorproduct([1 0; 0 -1], [0 1; 1 0]; tocache=false), [1 0; 0 -1]; tocache=false)
                op2 += 0.56/2*tensorproduct(tensorproduct([1 0; 0 -1], [1 0; 0 -1]; tocache=false), [1 0; 0 1]; tocache=false)
                op2 += 0.56*tensorproduct(tensorproduct([1 0; 0 1], [1 0; 0 -1]; tocache=false), [1 0; 0 -1]; tocache=false)
                op2 += 0.34/3*tensorproduct(tensorproduct([0 1; 1 0], [1 0; 0 1]; tocache=false), [1 0; 0 1]; tocache=false)
                op2 += 0.34/2*tensorproduct(tensorproduct([1 0; 0 1], [0 1; 1 0]; tocache=false), [1 0; 0 1]; tocache=false)
                op2 += 0.34*tensorproduct(tensorproduct([1 0; 0 1], [1 0; 0 1]; tocache=false), [0 1; 1 0]; tocache=false)
                op2 = TeNe.exp(0.1*op2, (2, 4, 6))
                isapprox(op, op2)
            end

            @test begin
                op = TeNe._create_trotter_gate_compressed(H, 0.1, 17)
                op2 = -1.4*tensorproduct(tensorproduct([1 0; 0 -1], [0 1; 1 0]; tocache=false), [1 0; 0 -1]; tocache=false)
                op2 += 0.56/2*tensorproduct(tensorproduct([1 0; 0 -1], [1 0; 0 -1]; tocache=false), [1 0; 0 1]; tocache=false)
                op2 += 0.56/2*tensorproduct(tensorproduct([1 0; 0 1], [1 0; 0 -1]; tocache=false), [1 0; 0 -1]; tocache=false)
                op2 += 0.34/3*tensorproduct(tensorproduct([0 1; 1 0], [1 0; 0 1]; tocache=false), [1 0; 0 1]; tocache=false)
                op2 += 0.34/3*tensorproduct(tensorproduct([1 0; 0 1], [0 1; 1 0]; tocache=false), [1 0; 0 1]; tocache=false)
                op2 += 0.34/2*tensorproduct(tensorproduct([1 0; 0 1], [1 0; 0 1]; tocache=false), [0 1; 1 0]; tocache=false)
                op2 = TeNe.exp(0.1*op2, (2, 4, 6))
                isapprox(op, op2)
            end
        end
    end

    # Uncompressed gates 
    @testset "_create_trotter_gate_uncompressed" begin 
        @testset "two sites" begin 
            H = OpList(Qubits(), 20)
            for i = 1:20
                add!(H, "x", i, 0.34)
            end
            for i = 1:19
                add!(H, ["z", "z"], [i, i+1], 0.56)
            end

            @test begin
                op = TeNe._create_trotter_gate_uncompressed(H, 0.1, 5, 1)
                op2 = TeNe.exp(0.34*0.1*[0.0 1; 1 0], (2))
                isapprox(op, op2)
                true
            end
        end
    end
end