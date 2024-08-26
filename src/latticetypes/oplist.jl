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

### Some basic properties 
TeNe.rank(::OpList) = 2
TeNe.dim(H::OpList) = dim(H.lt)
Base.length(H::OpList) = H.length
innerdims(H::OpList) = ntuple(x->dim(H), Val(length(H)))
outerdims(H::OpList) = ntuple(x->dim(H), Val(length(H)))
Base.eltype(H::OpList) = eltype(H.lt)

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

### Adding operator lists
"""
    add(ops1::OpList, ops2::OpList)
    +(ops1::OpList, ops2::OpList)

Join two oplists.
"""
function add(ops1::OpList, ops2::OpList)
    # Validation 
    if !(ops1.lt == ops2.lt)
        for i in eachindex(ops2.ops)
            for j in eachindex(ops2.ops[i])
                if !opin(ops1.lt, ops2.ops[i][j])
                    throw(ArgumentError("LatticeTypes do not match."))
                end
            end
        end
    end

    # Create new list 
    ops = OpList(ops1.lt, max(ops1.length, ops2.length))

    # Add the operators
    ops.sites = copy(ops1.sites)
    append!(ops.sites, ops2.sites)
    ops.ops = copy(ops1.ops)
    append!(ops.ops, ops2.ops)
    ops.coeffs = copy(ops1.coeffs)
    append!(ops.coeffs, ops2.coeffs)
    return ops
end
+(ops1::OpList, ops2::OpList) = add(ops1, ops2)

### Multiplication of opeartor lists 
function *(x::Number, y::OpList)
    y = deepcopy(y)
    for i = 1:length(y.coeffs)
        y.coeffs[i] *= x
    end
    return y
end

function *(x::OpList, y::Number)
    return *(y, x)
end

function /(x::OpList, y::Number)
    return *((1/y), x)
end

### Determine properties about the operator list
export siterange
"""
    siterange(oplist::OpList)

Determine the maximal interaction range within an operator list.
"""
function siterange(ops::OpList)
    rng = 1
    for sites in ops.sites
        rng = max(rng, sites[end]-sites[1]+1)
    end
    return rng
end

export siteindexs
"""
    siteindexs(ops::OpList, site::Int)

Determine the indexs of operators in a list which start at a given site.
"""
function siteindexs(ops::OpList, site::Int)
    idxs = []
    for i = 1:length(ops.sites)
        if ops.sites[i][1] == site
            push!(idxs, i)
        end
    end
    return idxs
end

### Conversions to actual operators
export totensor
"""
    totensor(ops::OpList, idx::Int, [paddingleft::Int, paddingright::Int]; kwargs...)

Construct the tensor from an operator in an OpList.
Use `paddingleft` and `paddingright` to add identity paddings to either side. 

# Optional Keyword Arguments 

    - `tocache::Bool=false`: Store the result in the cache?
"""
function totensor(ops::OpList, idx::Int, paddingleft::Int=0, paddingright::Int=0; tocache::Bool=false)
    # Fetch the relevent information
    opers = ops.ops[idx]
    sites = ops.sites[idx]
    #rng = min(siterange(ops), ops.length - sites[1] + 1)
    rng = sites[end] - sites[begin] + 1

    # Create the tensor through a tensor product
    prod = ones(eltype(ops.lt), )
    i = 1
    for site = Base.range(1-paddingleft, rng+paddingright)
        if i <= sites[end] && site > 0 && site <= rng && (sites[i] == sites[1] + site - 1)
            oper = opers[i]
            i += 1
        else
            oper = "id"
        end
        prod = tensorproduct(prod, op(ops.lt, oper); tocache = site==rng+paddingright ? tocache : true)
    end
    return ops.coeffs[idx]*prod
end

export sitetensor
"""
    sitetensor(ops::OpList, idx::Int)

Return the tensor for all operators starting at a site.
"""
function sitetensor(ops::OpList,  idx::Int)
    # Validate the idx
    (idx < 0 || idx > ops.length) && error("Site index is out of range.")

    # Get the indexs which start at the site
    idxs = siteindexs(ops, idx)
    length(idxs) == 0 && return false

    # Loop through adding them
    ten = 1
    for i = 1:length(idxs)
        if i == 1
            ten = totensor(ops, idxs[i])
        else
            ten += totensor(ops, idxs[i])
        end
    end

    return ten
end