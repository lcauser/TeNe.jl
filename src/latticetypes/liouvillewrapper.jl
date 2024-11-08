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
Operators: ["id", "x", "y", "z", "pu", "pd", "n", "s+", "s-"]
```
"""
function LiouvilleWrapper(lt::LatticeTypes)
    d = dim(lt)^2
    T = Base.promote_type(Base.eltype(lt), ComplexF64)
    new_lt = LatticeTypes(d; T=T)

    # Kronecker product of all states
    for name in lt.statenames 
        st = kron(state(lt, name), state(lt, name))
        add!(new_lt, name, st)
    end

    # Kronecker product of all operators 
    for name1 in lt.opnames 
        for name2 in lt.opnames
            oper = kron(op(lt, name1), transpose(op(lt, name2)))
            add!(new_lt, name1*"_"*name2, oper)
        end
    end

    return new_lt
end