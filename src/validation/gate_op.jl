function _gate_op_validation(U::AbstractGate, O::TensorNetworkOperator, sites)
    # Check sites 
    Olen = length(O)
    if length(sites) != length(U)
        throw(ArgumentError("The length of sites does not match the length of the gate."))
    end
    if any(map(j->(j>Olen || j <= 0), sites))
        throw(ArgumentError("The list of sites $(sites) does not fall between 1 and $(length(O))."))
    end

    # Check dimensions of sites 
    if any(map(j->innerdim(O, j)!=dim(U), sites))
        throw(ArgumentError("The StateVector has the wrong physical dimensions."))
    end
end

function _gate_op_validation(O::TensorNetworkOperator, U::AbstractGate, sites)
    # Check sites 
    Olen = length(O)
    if length(sites) != length(U)
        throw(ArgumentError("The length of sites does not match the length of the gate."))
    end
    if any(map(j->(j>Olen || j <= 0), sites))
        throw(ArgumentError("The list of sites $(sites) does not fall between 1 and $(length(O))."))
    end

    # Check dimensions of sites 
    if any(map(j->outerdim(O, j)!=dim(U), sites))
        throw(ArgumentError("The StateVector has the wrong physical dimensions."))
    end
end
