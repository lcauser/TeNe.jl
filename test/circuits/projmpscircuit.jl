@testset "projmpscircuit" begin
    @testset "bw-circuit" begin
        # The first test will create a state and apply a circuit.
        # It will test that the projection creates the correct MPSs.
        ψ = randommps(2, 10, 1)
        U = randombwcircuit(2, 10, 5; ϵ=1.0)
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

        # Test for building left blocks 
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

        # Test for building right blocks
        ol2s = []
        for i = 1:5
            movecenter!(projU, i)
            right = ones(Float64, 1, 1)
            for j = length(U.layers[end+1-i].sites):-1:0
                right = copy(TeNe._buildright(projU, j, right))
            end
            push!(ol2s, right[])
        end
        @test all(isapprox.(ol2s, ol1))

        # Test for product of gate 
        ol2s = []
        for i = 1:5
            for j = 1:length(U.layers[end+1-i].sites)
                movecenter!(projU, i, j)
                push!(ol2s, product(projU))
            end
        end
        @test all(isapprox.(ol2s, ol1))

        # Test projection 
        projU = ProjMPSCircuit(ψ′, U, ψ; cutoff=0.0)
        checks = []
        for i = 1:5
            for j = 1:length(U.layers[end+1-i].sites)
                movecenter!(projU, i, j)
                gate = creategate(conj(TeNe.project(projU)))
                makeunitary!(gate)
                gate2 = U.layers[end+1-i].gates[j]             
                U.layers[end+1-i].gates[j] = gate
                push!(checks, isapprox(1.0, product(projU)))
            end
        end
        @test all(checks)
    end

    @testset "3-bw-circuit" begin
        # The first test will create a state and apply a circuit.
        # It will test that the projection creates the correct MPSs.
        ψ = randommps(2, 10, 1)
        U = randombwcircuit(2, 10, 6, 3; ϵ=1.0)
        ϕ = randommps(2, 10, 4)
        ϕ′ = deepcopy(ϕ)
        ψ′ = applygates(U, ψ)
        ol1 = inner(ϕ, ψ′)
        projU = ProjMPSCircuit(ϕ, U, ψ; cutoff=0.0)
        ol2s = []
        for i = 2:6
            movecenter!(projU, i)
            ol2 = inner(conj(TeNe.topblock(projU, i-1)), TeNe.bottomblock(projU, i))
            push!(ol2s, ol2)
        end
        @test all(isapprox.(ol2s, ol1))

        # Test for building left blocks 
        ol2s = []
        for i = 1:6
            movecenter!(projU, i)
            left = ones(Float64, 1, 1)
            for j = 1:length(U.layers[end+1-i].sites)+1
                left = copy(TeNe._buildleft(projU, j, left))
            end
            push!(ol2s, left[])
        end
        @test all(isapprox.(ol2s, ol1))

        # Test for building right blocks
        ol2s = []
        for i = 1:6
            movecenter!(projU, i)
            right = ones(Float64, 1, 1)
            for j = length(U.layers[end+1-i].sites):-1:0
                right = copy(TeNe._buildright(projU, j, right))
            end
            push!(ol2s, right[])
        end
        @test all(isapprox.(ol2s, ol1))

        # Test for product of gate 
        ol2s = []
        for i = 1:6
            for j = 1:length(U.layers[end+1-i].sites)
                movecenter!(projU, i, j)
                push!(ol2s, product(projU))
            end
        end
        @test all(isapprox.(ol2s, ol1))

        # Test projection 
        projU = ProjMPSCircuit(ψ′, U, ψ; cutoff=0.0)
        checks = []
        for i = 1:6
            for j = 1:length(U.layers[end+1-i].sites)
                movecenter!(projU, i, j)
                gate = creategate(conj(TeNe.project(projU)))
                makeunitary!(gate)
                gate2 = U.layers[end+1-i].gates[j]             
                U.layers[end+1-i].gates[j] = gate
                push!(checks, isapprox(1.0, product(projU)))
            end
        end
        @test all(checks)
    end

    @testset "gap-circuit" begin
        # The first test will create a state and apply a circuit.
        # It will test that the projection creates the correct MPSs.
        ψ = randommps(2, 10, 1)
        ϕ = randommps(2, 10, 4)

        U = randombwcircuit(2, 10, 2)
        for i = range(1, 8, step=3)
            add!(U, TeNe._unitary_close_to_id(2, 2, 1.0), (i, i+2))
        end
        for i = range(2, 8, step=3)
            add!(U, TeNe._unitary_close_to_id(2, 2, 1.0), (i, i+2))
        end
        for i = range(3, 8, step=3)
            add!(U, TeNe._unitary_close_to_id(2, 2, 1.0), (i, i+2))
        end
        println(length(U.layers))

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

        # Test for building left blocks 
        ol2s = []
        println(ol1)
        for i = 1:5
            movecenter!(projU, i)
            left = ones(Float64, 1, 1)
            for j = 1:length(U.layers[end+1-i].sites)+1
                left = copy(TeNe._buildleft(projU, j, left))
            end
            println(left)
            push!(ol2s, left[])
        end
        @test all(isapprox.(ol2s, ol1))

        """
        # Test for building right blocks
        ol2s = []
        for i = 1:5
            movecenter!(projU, i)
            right = ones(Float64, 1, 1)
            for j = length(U.layers[end+1-i].sites):-1:0
                right = copy(TeNe._buildright(projU, j, right))
            end
            push!(ol2s, right[])
        end
        @test all(isapprox.(ol2s, ol1))

        # Test for product of gate 
        ol2s = []
        for i = 1:5
            for j = 1:length(U.layers[end+1-i].sites)
                movecenter!(projU, i, j)
                push!(ol2s, product(projU))
            end
        end
        @test all(isapprox.(ol2s, ol1))

        # Test projection 
        projU = ProjMPSCircuit(ψ′, U, ψ; cutoff=0.0)
        checks = []
        for i = 1:5
            for j = 1:length(U.layers[end+1-i].sites)
                movecenter!(projU, i, j)
                gate = creategate(conj(TeNe.project(projU)))
                makeunitary!(gate)
                gate2 = U.layers[end+1-i].gates[j]             
                U.layers[end+1-i].gates[j] = gate
                push!(checks, isapprox(1.0, product(projU)))
            end
        end
        @test all(checks)
        """
    end
end