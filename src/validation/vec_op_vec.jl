#=
    Validation for expectations of operators.
=#

function _vec_op_vec_validation(ψ::TensorNetworkVector, ϕ::TensorNetworkVector, Os::TensorNetworkOperator...)
    lengths = (length(ψ), length.(Os)..., length(ϕ))
    if !all(lengths .== lengths[begin])
        throw(ArgumentError("Length mismatch for inner product with lengths $(lengths)."))
    end
    outers = (dims(ψ), outerdims.(Os)...)
    inners = (innerdims.(Os)..., dims(ϕ))
    
    for i in eachindex(outers)
        check = inners[i] .== outers[i]
        if !all(check)
            throw(ArgumentError("Dimensions mistmatch for contraction number $(i) in the inner product."))
        end
    end
    return true
end

function _inner_validation(args...)
    if rank(args[begin]) != 1
        throw(ArgumentError("The first term in the inner product must be rank-1."))
    end
    for i in Base.OneTo(length(args)-2)
        if rank(args[begin+i]) != 2
            throw(ArgumentError("The middle terms in the inner product must be rank-2."))
        end
    end
    if rank(args[end]) != 1
        throw(ArgumentError("The last term in the inner product must be rank-1."))
    end
end