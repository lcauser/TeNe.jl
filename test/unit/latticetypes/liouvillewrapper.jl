@testset "latticesites-liouvillewrapper" begin 
    for lt in [Qubits(), Bosons(3)]
        lt2 = LiouvilleWrapper(lt)
        @test dim(lt2) == dim(lt)^2
        @test Base.eltype(lt2) == ComplexF64
        for name in lt.statenames
            st = state(lt, name)
            @test state(lt2, name) == kron(st, st)
        end
        for name in lt.opnames
            oper = op(lt, name)
            @test op(lt2, name) == kron(oper, transpose(oper))
        end
    end
end