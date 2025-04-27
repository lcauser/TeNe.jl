@testset "oplist" begin
    lt = Qubits()
    @test begin
        ops = OpList(lt, 20)
        true
    end

    @test begin
        ops = OpList(lt, 20)
        ops2 = deepcopy(ops)
        Base.mightalias(ops.ops, ops2.ops) == false
    end

    @test begin
        ops = OpList(lt, 20)
        for i = 1:19
            add!(ops, ["z", "z"], [i, i+1], rand())
        end
        add!(ops, "z", 1, rand())
        true
    end

    @test begin
        ops1 = OpList(Qubits(), 20)
        for i = 1:19
            add!(ops1, ["z", "z"], [i, i+1], rand())
        end

        ops2 = OpList(Qubits(), 10)
        for i = 1:10
            add!(ops2, "x", i, rand())
        end
        ops = 2*ops1 + ops2 / 2
        true
    end

    @test begin
        ops = OpList(Qubits(), 20)
        for i = 1:19
            add!(ops, ["z", "z"], [i, i+1], rand())
        end
        siterange(ops) == 2
    end

    @test begin
        ops = OpList(Qubits(), 20)
        for i = 1:19
            add!(ops, ["z", "z"], [i, i+1], 1)
        end
        ten1 = totensor(ops, 2)
        ten2 = tensorproduct([1 0; 0 -1], [1 0; 0 -1])
        isapprox(ten1, ten2)
    end

    # TODO: sitetensor currently broken, fix later
    #=
    @test begin
        ops = OpList(Qubits(), 20)
        for i = 1:19
            add!(ops, ["z", "z"], [i, i+1], 1)
        end
        for i = 1:20
            add!(ops, ["x"], [i], 1)
        end
        ten1 = sitetensor(ops, 2)
        ten2 = tensorproduct([1 0; 0 -1], [1 0; 0 -1]; tocache=false)
        ten2 += tensorproduct([0 1; 1 0], [1 0; 0 1])
        isapprox(ten1, ten2)
    end
    =#

end
