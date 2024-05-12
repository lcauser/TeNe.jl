@testset "circuitlayers" begin 
    ψ = randommps(2, 20, 1)
    layer = CircuitLayer(2, 20)
    for i = 1:10
        add!(layer, makeunitary(creategate(rand(2, 2, 2, 2))), (2*i-1, 2*i))
    end
    applygates!(layer, ψ; cutoff=1e-16)
    println([bonddim(ψ, i) for i=1:19])
    println(center(ψ))

    layer2 = CircuitLayer(2, 20)
    for i = 1:9
        add!(layer2, makeunitary(creategate(rand(2, 2, 2, 2))), (2*i, 2*i+1))
    end
    applygates!(layer2, ψ; cutoff=1e-16)
    println([bonddim(ψ, i) for i=1:19])
    println(center(ψ))
    println(inner(ψ, ψ))
end