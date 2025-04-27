#=
    The Qubits type is a general purpose lattice type for subsyltems with two degress
    of freedom.
=#

export Qubits
"""
    Qubits()

A lattice site with the information for qubits.

```jldoctest
julia> qubits = Qubits()
Lattice Type
Dimension: 2
States: ["up", "dn", "1", "0", "+", "-", "s", "as"]
Operators: ["id", "x", "y", "z", "pu", "pd", "n", "s+", "s-"]
```
"""
function Qubits()
    # Create the sitetype
    lt = LatticeTypes(2)

    # Add the ltate
    add!(lt, "up", [1, 0])
    add!(lt, "dn", [0, 1])
    add!(lt, "1", [1, 0])
    add!(lt, "0", [0, 1])
    add!(lt, "+", [0.5^0.5, 0.5^0.5])
    add!(lt, "-", [0.5^0.5, -0.5^0.5])
    add!(lt, "s", [0.5^0.5, 0.5^0.5])
    add!(lt, "as", [0.5^0.5, -0.5^0.5])

    # Add operators
    add!(lt, "id", [1 0; 0 1])
    add!(lt, "x", [0 1; 1 0])
    add!(lt, "y", [0 -1im; 1im 0])
    add!(lt, "z", [1 0; 0 -1])
    add!(lt, "pu", [1 0; 0 0])
    add!(lt, "pd", [0 0; 0 1])
    add!(lt, "n", [1 0; 0 0])
    add!(lt, "s+", [0 1; 0 0])
    add!(lt, "s-", [0 0; 1 0])

    return lt
end
