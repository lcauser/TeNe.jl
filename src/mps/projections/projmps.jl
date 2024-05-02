#=
    Projections for inner products of MPS
=#

mutable struct ProjMPS{Q<:Number}
    objects::Vector{Union{MPS, MPO}}
    lefts::Vector{<:AbstractArray}
    rights::Vector{<:AbstractArray}
    center::Int 
    λ::Q
end

function ProjMPS(ψ::MPS, args::Union{MPS, MPO}; λ::Number=1.0, center::Int=1)
    _inner_validation(ψ, args...)
    _vec_op_vec_validation(ψ, args[end], args[begin:end-1]...)

end

function buildleft!(projψ::ProjMPS, site::Int)
    # Build on the previous left block
    ten = leftblock(projψ, site-1)
    ten = contract(ten, projψ.objects[begin][site], 1, 1, false, !isconj(projψ.objects[begin]))
    for i = Base.OneTo(length(projψ.objects)-2)
        ten = contract(ten, projψ.objects[begin+i][site], (1, ndims(ten)-1),
            (1, istranspose(projψ.objects[begin+i]) ? 3 : 2), false, isconj(projψ.objects[begin+i]))
    end
    ten = contract(ten, projψ.objects[end][site], (1, ndims(ten)-1), (1, 2), false,
        isconj(projψ.objects[end]))

    # Save to memory
    if size(ten) == size(leftblock(projψ, site))
        projψ.lefts[site] .= ten
    else
        projψ.lefts[site] = copy(ten)
    end
end