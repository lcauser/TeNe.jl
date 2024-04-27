#=
    Operator lists are a high-level interface for strings of operators,
    such as Hamiltonians.
=#

mutable struct OpList{Q<:Number}
    lt::LatticeTypes
    length::Int
    ops::Vector{Vector{String}}
    sites::Vector{Vector{Int}}
    coeffs::Vector{Q}
end

### Initiate
export OpList

"""
    OpList(lt::LatticeTypes, length::Int)

Create an operator list from operators in a lattice type.

# Examples 
```julia-repl
julia> ops = OpList(Qubits(), 10);
```
"""
function OpList(lt::LatticeTypes, length::Int)
    return OpList{eltype(lt)}(lt, length, Vector{String}[], Vector{Int}[], eltype(lt)[])
end

### Copys 
function Base.copy(ops::OpList)
    ops2 = OpList(ops.lt, ops.length)
    ops2.sites = copy(ops.sites)
    ops2.ops = copy(ops.ops)
    ops2.coeffs = copy(ops.coeffs)
    return ops2
end

### Add to the list
export add!
"""
    add!(ops::OpList, ops::Vector{String}, sites::Vector{Int},
         coeff<:Number = 1)
    add!(ops::OpList, op::String, site::Int, coeff<:Number = 1)

Add an operator to the list defined by local operators at given sites.
"""
function add!(ops::OpList, opers::Vector{String}, sites::Vector{Int},
              coeff::Q = 1.0) where {Q<:Number}
    # Validate the data
    if length(opers) != length(sites)
        throw(ArgumentError("Mismatch in the length of opers and sites."))
    end

    # Order the operators and sites
    perms = sortperm(sites)
    sites = sites[perms]
    opers = opers[perms]

    # More validations
    for i in eachindex(sites)
        if !(opers[i] in ops.lt.opnames)
            throw(ArgumentError("The operator $(opers[i]) is undefined."))
        end
        if 0 > sites[i] || sites[i] > ops.length
            throw(ArgumentError("The sites must be between 1 and $(ops.length)."))
        end
        if sum([sites[i] == site2 for site2 in sites]) > 1
            throw(ArgumentError("There are two or more operators on the same site."))
        end
    end

    # Add to list
    push!(ops.ops, opers)
    push!(ops.sites, sites)
    push!(ops.coeffs, coeff)
end

function add!(ops::OpList, oper::String, site::Int, coeff::Q = 1.0) where {Q<:Number}
    add!(ops, [oper], [site], coeff)
end