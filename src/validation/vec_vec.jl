#=
    Validation for tensor network inner products of Vectorial objects with
    Vectorial objects
=#

function _vec_vec_validation(ψ::TensorNetworkVector, ϕ::TensorNetworkVector)
    if length(ψ) != length(ϕ)
        throw(
            ArgumentError(
                "Length mismatch for $(typeof(ψ)) and $(typeof(ϕ)) with lengths $(length(ψ)) and $(length(ϕ)).",
            ),
        )
    end
    check = dims(ψ) .== dims(ϕ)
    if !all(check)
        check = check .== false
        if all(check)
            throw(ArgumentError("Dimensions mistmatch for $(typeof(ψ)) and $(typeof(ϕ))."))
        else
            throw(
                ArgumentError(
                    "Dimensions mistmatch for $(typeof(ψ)) and $(typeof(ϕ)) at sites $(findall(check .== false)).",
                ),
            )
        end
    end
    return true
end
