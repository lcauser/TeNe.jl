#=
    Validation for tensor network products for operators.
=#

function _op_op_validation(O1::TensorNetworkOperator, O2::TensorNetworkOperator)
    if length(O1) != length(O2)
        throw(ArgumentError("Length mismatch for $(typeof(O1)) and $(typeof(O2)) with lengths $(length(O1)) and $(length(O2))."))
    end
    check = outerdims(O1) .== innerdims(O2)
    if !all(check)
        check = check .== false 
        if all(check)
            throw(ArgumentError("Dimensions mistmatch for $(typeof(O1)) and $(typeof(O2))."))
        else
            throw(ArgumentError("Dimensions mistmatch for $(typeof(O1)) and $(typeof(O2)) at sites $(findall(check .== false))."))
        end
    end
    return true
end

function _op_op_validation(O::TensorNetworkOperator, O1::TensorNetworkOperator, O2::TensorNetworkOperator)
    _op_vec_validation(O1, O2)
    if length(O) != length(O1)
        throw(ArgumentError("Length mismatch for $(typeof(O)) and $(typeof(O1)) with lengths $(length(O)) and $(length(O1))."))
    end
    check = innerdims(O) .== innerdims(O1)
    if !all(check)
        check = check .== false 
        if all(check)
            throw(ArgumentError("Dimensions mistmatch for $(typeof(O)) and $(typeof(O1))."))
        else
            throw(ArgumentError("Dimensions mistmatch for $(typeof(O)) and $(typeof(O1)) at sites $(findall(check .== false))."))
        end
    end
    check = outerdims(O) .== outerdims(O2)
    if !all(check)
        check = check .== false 
        if all(check)
            throw(ArgumentError("Dimensions mistmatch for $(typeof(O)) and $(typeof(O2))."))
        else
            throw(ArgumentError("Dimensions mistmatch for $(typeof(O)) and $(typeof(O2)) at sites $(findall(check .== false))."))
        end
    end
    return true
end
