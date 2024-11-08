@testset "latticesites-liouvillewrapper" begin 
    for lt in [Qubits(), Bosons(3)]
        lt2 = LiouvilleWrapper(lt)
        @test dim(lt2) == dim(lt)^2
        @test Base.eltype(lt2) == ComplexF64
        for name in lt.statenames
            st = state(lt, name)
            @test state(lt2, name) == kron(st, st)
        end
        for name1 in lt.opnames
            for name2 in lt.opnames
                oper1 = op(lt, name1)
                oper2 = op(lt, name2)
                @test op(lt2, name1*"_"*name2) == kron(oper1, transpose(oper2))
            end
        end
    end
end