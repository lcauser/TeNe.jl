@testset "circuitlayer" begin
    @testset "mps-connector" begin
        # Test for nearest neighbour additions
        @test begin
            layer = CircuitLayer(2, 20)
            for i = 20:-2:2
                add!(layer, TeNe._unitary_close_to_id(2, 2, 0.5), (i, i-1), CircuitMPS())
            end

            j = 1
            check = true
            for i in eachindex(layer.sites)
                if layer.sites[i][1] != j
                    check = false
                end
                j += 2
            end
            check == true
        end
    end

end
