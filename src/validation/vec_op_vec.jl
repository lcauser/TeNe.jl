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