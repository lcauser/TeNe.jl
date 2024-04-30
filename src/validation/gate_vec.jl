function _gate_vec_validation(U::AbstractGate, ψ::TensorNetworkVector, sites)
    # Check sites 
    ψlen = length(ψ)
    if length(sites) != length(U)
        throw(ArgumentError("The length of sites does not match the length of the gate."))
    end
    if any(map(j->(j>ψlen || j <= 0), sites))
        throw(ArgumentError("The list of sites $(sites) does not fall between 1 and $(length(ψ))."))
    end

    # Check dimensions of sites 
    if any(map(j->dim(ψ, j)!=dim(U), sites))
        throw(ArgumentError("The StateVector has the wrong physical dimensions."))
    end
end