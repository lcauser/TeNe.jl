#=
    Validation for tensor network traces of operators.
=#

function _op_trace_validation(Os::TensorNetworkOperator...)
    lengths = length.(Os)
    if !all(lengths .== lengths[begin])
        throw(ArgumentError("Length mismatch for trace with lengths $(lengths)."))
    end

    for i in eachindex(Os)
        outers = outerdims(Os[i])
        inners = innerdims(Os[i == lastindex(Os) ? firstindex(Os) : i+1])
        check = inners .== outers
        if !all(check)
            throw(
                ArgumentError(
                    "Dimensions mistmatch for contraction number $(i) in the inner product.",
                ),
            )
        end
    end
    return true
end
