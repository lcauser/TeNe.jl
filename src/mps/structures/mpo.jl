#=
    Matrix product operators are typically used to represent Hamiltonians or operators
    of 1D quantum many-body systems, but can more generally be used to describe a rank-2
    tensor of N discrete degrees of freedom.
=#

### Traits for MPOs 
export conj, transpose, adjoint, istranspose
# Transpose
struct TransposeMPO <: GMPSTrait
    MPS::GMPS{2}
end
TeNe.transpose(O::GMPS{2}) = TransposeMPO(O)
TeNe.transpose(O::TransposeMPO) = O.MPS
istranspose(ψ::TransposeMPO) = true

# Adjoint
struct AdjointMPO <: GMPSTrait
    MPS::GMPS{2}
end
TeNe.adjoint(O::GMPS{2}) = AdjointMPO(O)
TeNe.adjoint(O::AdjointMPO) = O.MPS
isconj(ψ::AdjointMPO) = true
istranspose(ψ::AdjointMPO) = true

# Trait rules
TeNe.conj(O::TransposeMPO) = AdjointMPO(O.MPS)
TeNe.conj(O::AdjointMPO) = TransposeMPO(O.MPS)
TeNe.transpose(O::ConjGMPS{2}) = AdjointMPO(O.MPS)
TeNe.transpose(O::AdjointMPO) = ConjGMPS(O.MPS)
TeNe.adjoint(O::ConjGMPS{2}) = TransposeMPO(O.MPS)
TeNe.adjoint(O::TransposeMPO) = ConjGMPS(O.MPS)


const MPO = Union{GMPS{2},ConjGMPS{2},TransposeMPO,AdjointMPO}
export MPO


export ismpo
"""
    ismpo(O)

Check to see if an object is an MPO.
"""
function ismpo(O)
    return typeof(O) <: MPO
end

# Get the size of physical dimensions
export innerdim, outerdim
innerdim(O::Union{GMPS{2},ConjGMPS{2}}, site::Int) = size(O[site], 2)
innerdim(O::Union{AdjointMPO,TransposeMPO}, site::Int) = size(O[site], 3)
outerdim(O::Union{GMPS{2},ConjGMPS{2}}, site::Int) = size(O[site], 3)
outerdim(O::Union{AdjointMPO,TransposeMPO}, site::Int) = size(O[site], 2)
innerdims(O::Union{GMPS{2},ConjGMPS{2}}) = dims(O, 1)
innerdims(O::Union{AdjointMPO,TransposeMPO}) = dims(O, 2)
outerdims(O::Union{GMPS{2},ConjGMPS{2}}) = dims(O, 2)
outerdims(O::Union{AdjointMPO,TransposeMPO}) = dims(O, 1)

# Return the correct indices
innerind(::Union{GMPS{2},ConjGMPS{2}}, ::Int = 0) = 2
innerind(::Union{AdjointMPO,TransposeMPO}, ::Int = 0) = 3
outerind(::Union{GMPS{2},ConjGMPS{2}}, ::Int = 0) = 3
outerind(::Union{AdjointMPO,TransposeMPO}, ::Int = 0) = 2

### Initalising MPOs 
export randommpo, productmpo
"""
    MPO(dim::Int, length::Int)

Create an MPO with physical dimension `dim` and `length` sites.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensors.
"""
function MPO(dim::Int, length::Int; kwargs...)
    return GMPS(2, dim, length; kwargs...)
end


"""
    randommpo(dim::Int, length::Int, bonddim::Int; kwargs...)

Create an MPO with dimensions `dim`, size `length` and bond dimension `bonddim`,
with random tensors.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensors.
"""
function randommpo(dim::Int, length::Int, bonddim::Int; kwargs...)
    return randomgmps(2, dim, length, bonddim; kwargs...)
end


"""
    productmpo(N::Int, A::AbstractArray; kwargs...)

Create a product MPO of size `N`, composed of array `A`. 
`A` can be a matrix or rank-4 tensor.

# Optional Keyword Arguments

    - `T::Type=ComplexF64`: The element type for the tensors.
"""
function productmpo(N::Int, A::AbstractArray; T::Type = ComplexF64)
    if ndims(A) == 2
        if size(A, 1) != size(A, 2)
            throw(ArgumentError("A must be a square matrix!"))
        end
        O = MPO(size(A, 1), N; T = T)
        for i in Base.OneTo(N)
            O[i] = reshape(copy(A), (1, size(A)..., 1))
        end
    elseif ndims(A) == 4
        if size(A, 2) != size(A, 3)
            throw(ArgumentError("Array must have physical dimensions of the same length!"))
        end
        if size(A, 1) != size(A, 4)
            throw(ArgumentError("Array must have equal bond dimensions!"))
        end
        O = MPO(size(A, 2), N; T = T)
        O[1] = A[1:1, :, :, :]
        for i in Base.range(1+firstindex(O), lastindex(O))
            O[i] = copy(A)
        end
        O[lastindex(O)] = A[:, :, :, end:end]
    else
        throw(ArgumentError("You must provide an array with just two or four dimensions."))
    end
    return O
end

"""
    productmpo(lt::LatticeTypes, ops::AbstractVector{String})

Create an MPO from a string of operators.

# Example 

```julia-repl
julia> lt = Qubits();
julia> O = productmpo(lt, ["x" for _ = 1:20]);
```
"""
function productmpo(lt::LatticeTypes, ops::AbstractVector{String})
    O = MPO(dim(lt), length(ops); T = eltype(lt))
    for i in eachindex(ops)
        O[i][1, :, :, 1] .= op(lt, ops[i])
    end
    movecenter!(O, firstindex(O))
    return O
end

"""
    MPO(O::StateOperator; kwargs...)

Write a StateOperator as a MPO.

# Optional Keyword Arguments
    
    - `cutoff::Float64=0.0`: Truncation criteria to reduce the bond dimension.
    Good values range from 1e-8 to 1e-14.
    - `mindim::Int=1`: Minimum dimension for the truncation.
    - `maxdim::Int=0`: Maximum bond dimension for truncation. Set to 0 to have
    no limit.

# Examples
```julia-repl
julia> ψ = randomsv(2, 10);
julia> ψ = MPS(ψ; cutoff=1e-12);
```
"""
function MPO(O::GStateTensor{2}; kwargs...)
    return GMPS(O; kwargs...)
end
MPO(O::ConjGStateTensor{2}; kwargs...) = conj(GMPS(O.StateTensor; kwargs...))
MPO(O::TransposeStateOperator; kwargs...) = transpose(GMPS(O.StateTensor; kwargs...))
MPO(O::AdjointStateOperator; kwargs...) = adjoint(GMPS(O.StateTensor; kwargs...))
