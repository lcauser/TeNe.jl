#=
    Truncated bosonic modes
=#

export Bosons
"""
    Bosons(dim::int)

Create a lattice of truncated bosonic modes with maximum dimension `dim`.

```jldoctest
julia> bosons = Bosons(3)
Lattice Type
Dimension: 3
States: ["0", "1", "2"]
Operators: ["id", "n", "0", "1", "2", "a", "adag"]
```
"""
function Bosons(dim::Int)
    # Create the sitetype
    T = ComplexF64
    lt = LatticeTypes(dim; T=T)

    # Add the states
    for i in Base.OneTo(dim)
        state = zeros(T, dim)
        state[i] = 1
        add!(lt, string(i-1), state)
    end

    # Add operators 
    add!(lt, "id", diagm(ones(T, dim)))
    add!(lt, "n", diagm(Vector{T}(Base.range(0, dim-1))))
    for i in Base.OneTo(dim)
        oper = zeros(T, dim, dim)
        oper[i, i] = 1
        add!(lt, string(i-1), oper)
    end
    a = zeros(T, dim, dim)
    adag = zeros(T, dim, dim)
    for i in Base.OneTo(dim-1)
        a[i, i+1] = sqrt(i)
        adag[i+1, i] = sqrt(i)
    end
    add!(lt, "a", a)
    add!(lt, "adag", adag)
    add!(lt, "n2", adag*adag*a*a)

    return lt
end