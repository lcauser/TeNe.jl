#=
    Validation for tensor network products for operators and vectors.
=#

function _op_vec_validation(O::TensorNetworkOperator, ψ::TensorNetworkVector)
    if length(ψ) != length(O)
        throw(
            ArgumentError(
                "Length mismatch for $(typeof(O)) and $(typeof(ψ)) with lengths $(length(O)) and $(length(ψ)).",
            ),
        )
    end
    check = dims(ψ) .== outerdims(O)
    if !all(check)
        check = check .== false
        if all(check)
            throw(ArgumentError("Dimensions mistmatch for $(typeof(O)) and $(typeof(ψ))."))
        else
            throw(
                ArgumentError(
                    "Dimensions mistmatch for $(typeof(O)) and $(typeof(ψ)) at sites $(findall(check .== false)).",
                ),
            )
        end
    end
    return true
end

function _op_vec_validation(
    ϕ::TensorNetworkVector,
    O::TensorNetworkOperator,
    ψ::TensorNetworkVector,
)
    _op_vec_validation(O, ψ)
    if length(ψ) != length(ϕ)
        throw(
            ArgumentError(
                "Length mismatch for $(typeof(ϕ)) and $(typeof(ψ)) with lengths $(length(ϕ)) and $(length(ψ)).",
            ),
        )
    end
    check = innerdims(O) .== dims(ϕ)
    if !all(check)
        check = check .== false
        if all(check)
            throw(ArgumentError("Dimensions mistmatch for $(typeof(ϕ)) and $(typeof(O))."))
        else
            throw(
                ArgumentError(
                    "Dimensions mistmatch for $(typeof(ϕ)) and $(typeof(O)) at sites $(findall(check .== false)).",
                ),
            )
        end
    end
    return true
end
