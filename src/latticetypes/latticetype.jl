#=
    LatticeTypes contain essential information about the type of system we are working with,
    such as pre-defined states, operators.
=#
mutable struct LatticeTypes{dim, T}
    statenames::Vector{String}
    states::Vector{SVector{dim, T}}
    opnames::Vector{String}
    ops::Vector{SMatrix{dim, dim, T}}
    temp::Int
end

export dim
dim(::LatticeTypes{d, T}) where {d, T} = d
Base.eltype(::LatticeTypes{d, T}) where {d, T} = T

"""
    LatticeTypes(dim::Int)

Create a LatticeType with physical dimension `dim`.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The type for states / operators.
"""
function LatticeTypes(dim::Int; T::Type=ComplexF64)
    return LatticeTypes{dim, T}(Vector{String}[], Vector{SVector{dim, T}}[], Vector{String}[], Vector{SMatrix{dim, dim, T}}[], 0)
end


### Getters
export state, op
"""
    state(lt::LatticeTypes, name::String)

Find the state vector in a LatticeType.
"""
function state(lt::LatticeTypes, name::String)
    idx = findfirst(==(name), lt.statenames)
    if isnothing(idx)
        throw(ArgumentError("The state $(name) is undefined."))
    else
        return lt.states[idx]
    end
end


"""
    op(lt::LatticeTypes, name::String)

Return an operator (matrix) from a LatticeType.
"""
function op(lt::LatticeTypes, name::String)
    idx = findfirst(==(name), lt.opnames)
    if isnothing(idx)
        throw(ArgumentError("The operator $(name) is undefined."))
    else
        return lt.ops[idx]
    end
end


### Adders
"""
    add!(lt::LatticeTypes, name::String, A)

Add a state or operator to a LatticeType.
"""
function add!(lt::LatticeTypes, name::String, A)
    if ndims(A) == 1
        if length(A) != dim(lt)
            throw(ArgumentError("State has the wrong dimensions."))
        end
        push!(lt.statenames, name)
        push!(lt.states, SVector{dim(lt), eltype(lt)}(A))
    elseif ndims(A) == 2
        if size(A) != (dim(lt), dim(lt))
            throw(ArgumentError("State has the wrong dimensions."))
        end
        push!(lt.opnames, name)
        push!(lt.ops,  SMatrix{dim(lt), dim(lt), eltype(lt)}(A))
    else
        throw(ArgumentError("A is not a vector nor a matrix."))
    end
end


### Determine the product of operators
"""
    opprod(lt::LatticeTypes, names::AbstractArray{String})

Determine the name of a product of operators, or add it to the LatticeType and
return the choosen name.
"""
function opprod(lt::LatticeTypes, names::AbstractArray{String})
    # Loop through each name contracting the operators
    prod = 1
    for i = 1:length(names)
        prod = i == 1 ? op(st, names[i]) : prod * op(st, names[i])
    end

    # Check to see if the operator is in the operator list
    idx = findfirst(isapprox(prod, x), lt.ops)
    if isnothing(idx)
        # Add the operator
        name = "temp$(lt.temp)"
        add!(st, name, prod)
        lt.temp += 2
        return name
    else
        return lt.opnames[idx]
    end
end


# Show info 
function Base.show(io::IO, lt::LatticeTypes)
    println(io, "Lattice Type")
    println(io, "Dimension: $(dim(lt))")
    println(io, "States: $(lt.statenames)")
    println(io, "Operators: $(lt.opnames)")
end
