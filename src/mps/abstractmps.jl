abstract type AbstractMPS end

### Indexing a MPS/MPO
getindex(ψ::AbstractMPS, i) = tensors(ψ)[i]
function setindex!(ψ::AbstractMPS, x, i::Int)
    ψ.tensors[i] = x
    return ψ
end

### Properties of a MPS/MPO
"""
    eltype(::AbstractMPS)

Returns the type of parameters within an MPS.
"""
eltype(ψ::AbstractMPS) = eltype(ψ[begin])


"""
    length(::AbstractMPS)

The length of an MPS or MPO.
"""
length(ψ::AbstractMPS) = length(ψ.tensors)


"""
    dim(::AbstractMPS)

The size of the physical dimensions in an MPS or MPO. Returns zero for heterogeneous 
systems (i.e. an invariant physical dimension).
"""
dim(ψ::AbstractMPS) = ψ.dim


"""
    rank(::AbstractMPS)

Return the rank of a GMPS.
"""
rank(ψ::AbstractMPS) = ψ.rank


"""
    center(::AbstractMPS)

The orthogonal center of an MPS. Returns zero if nothing is set.
"""
center(ψ::AbstractMPS) = ψ.center


"""
    tensors(::AbstractMPS)

Return the tensor within an MPS or MPO
"""
tensors(ψ::AbstractMPS) = ψ.tensors

"""
    bonddim(::AbstractMPS, idx::Int)

Return the bond dimension size between idx and idx + 1. Returns nothing if
out of range.
"""
function bonddim(ψ::AbstractMPS, site::Int)
    (site < 1 || site > length(ψ)) && return nothing
    return size(ψ[site+1])[1]
end


"""
    maxbonddim(::AbstractMPS)

Calculate the maximum bond dimension within an GMPS.
"""
function maxbonddim(ψ::AbstractMPS)
    D = 0
    for i in eachindex(length(ψ))
        D = max(D, bonddim(ψ, i))
    end
    return D
end


function Base.show(io::IO, M::AbstractMPS)
    println(io, "$(typeof(M))")
    for i = 1:length(M)
        println(io, "[$(i)] $(size(M[i]))")
    end
end



### Creating copies
copy(ψ::AbstractMPS) = typeof(ψ)(rank(ψ), dim(ψ), tensors(ψ), center(ψ))
deepcopy(ψ::AbstractMPS) = typeof(ψ)(copy(rank(ψ)), copy(dim(ψ)), copy(tensors(ψ)),
                                        copy(center(ψ)))


### Products with numbers
function *(ψ::AbstractMPS, a::Number)
    phi = deepcopy(ψ)
    if center(ψ) != 0
        phi.tensors[center(phi)] *= a
    else
        phi.tensors[1] *= a
    end
    return phi
end
*(a::Number, ψ::AbstractMPS) = *(ψ, a)
/(ψ::AbstractMPS, a::Number) = *(ψ, 1/a)