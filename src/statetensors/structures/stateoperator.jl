#=
    State operators are an exact representation of matrices for many-body systems with discrete
    degrees of freedom.
=#

### Traits for State Operators 

# Transpose
struct TransposeStateOperator <: GStateTensorTrait
    StateTensor::GStateTensor{2}
end
TeNe.transpose(O::GStateTensor{2}) = TransposeStateOperator(O)
TeNe.transpose(O::TransposeStateOperator) = O.StateTensor
istranspose(ψ::TransposeStateOperator) = true

# Adjoint
struct AdjointStateOperator <: GStateTensorTrait
    StateTensor::GStateTensor{2}
end
TeNe.adjoint(O::GStateTensor{2}) = AdjointStateOperator(O)
TeNe.adjoint(O::AdjointStateOperator) = O.StateTensor
isconj(ψ::AdjointStateOperator) = true
istranspose(ψ::AdjointStateOperator) = true

# Trait rules
TeNe.conj(O::TransposeStateOperator) = AdjointStateOperator(O.StateTensor)
TeNe.conj(O::AdjointStateOperator) = TransposeStateOperator(O.StateTensor)
TeNe.transpose(O::ConjGStateTensor{2}) = AdjointStateOperator(O.StateTensor)
TeNe.transpose(O::AdjointStateOperator) = ConjGStateTensor(O.StateTensor)
TeNe.adjoint(O::ConjGStateTensor{2}) = TransposeStateOperator(O.StateTensor)
TeNe.adjoint(O::TransposeStateOperator) = ConjGStateTensor(O.StateTensor)