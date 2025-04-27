#=
    A wrapper for that maps a LatticeType space onto a Liouville space.
=#

export LiouvilleWrapper
"""
    LiouvilleWrapper(lt::LatticeTypes)

Takes a LatticeType for a Hilbert space as an argument, and constructs a new LatticeType
that lives in the associate Liouville space.

The mapping is done according to Choi's isomorphism:
    - States: ρ = |ψ><ψ| → |ψ> ⊗ |ψ>
    - Operators: O1 ⋅ O2 -> O1 ⊗ transpose(O2)

This will take the product of all states defined in lt, and the kronecker product of
all combinations of operators.

# Examples

```jldoctest
julia> qubits = LiouvilleWrapper(Qubits())
Lattice Type
Dimension: 4
States: ["up", "dn", "1", "0", "+", "-", "s", "as"]
Operators: ["id_id", "id", "id_x", "id_y", "id_z", "id_pu", "id_pd", "id_n", "id_s+", "id_s-", "x_id", "x_x", "x_y", "x_z", "x_pu", "x_pd", "x_n", "x_s+", "x_s-", "y_id", "y_x", "y_y", "y_z", "y_pu", "y_pd", "y_n", "y_s+", "y_s-", "z_id", "z_x", "z_y", "z_z", "z_pu", "z_pd", "z_n", "z_s+", "z_s-", "pu_id", "pu_x", "pu_y", "pu_z", "pu_pu", "pu_pd", "pu_n", "pu_s+", "pu_s-", "pd_id", "pd_x", "pd_y", "pd_z", "pd_pu", "pd_pd", "pd_n", "pd_s+", "pd_s-", "n_id", "n_x", "n_y", "n_z", "n_pu", "n_pd", "n_n", "n_s+", "n_s-", "s+_id", "s+_x", "s+_y", "s+_z", "s+_pu", "s+_pd", "s+_n", "s+_s+", "s+_s-", "s-_id", "s-_x", "s-_y", "s-_z", "s-_pu", "s-_pd", "s-_n", "s-_s+", "s-_s-"]
```
"""
function LiouvilleWrapper(lt::LatticeTypes)
    d = dim(lt)^2
    T = Base.promote_type(Base.eltype(lt), ComplexF64)
    new_lt = LatticeTypes(d; T = T)

    # Kronecker product of all states
    for name in lt.statenames
        st = kron(state(lt, name), state(lt, name))
        add!(new_lt, name, st)
    end

    # Kronecker product of all operators 
    for name1 in lt.opnames
        for name2 in lt.opnames
            oper = kron(op(lt, name1), transpose(op(lt, name2)))
            name = name1*"_"*name2
            add!(new_lt, name, oper)
            if name == "id_id"
                add!(new_lt, "id", oper)
            end
        end
    end

    return new_lt
end
