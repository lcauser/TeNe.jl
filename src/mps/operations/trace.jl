#=
    Traces of MPOs.

    Need to add contraction checks.
=#

export trace


"""
    trace(Os::MPO...)

Compute the trace of the product of many MPOs.
"""
function TeNe.trace(Os::MPO...)
    _op_trace_validation(Os...)
    if length(Os) == 1
        return _mpo_trace(Os...)
    else
        return _mpo_mpo_trace(Os...)
    end
end

function _mpo_trace(O::MPO)
     # Type info...
     T = eltype(O)
     conjO = isconj(O)

     # Do the contraction 
    block = cache(T, size(O[begin], 1), 2, 1) .= 1
    for i in eachindex(O)
        block = contract(block, trace(O[i], 2, 3), 1, 1, false, conjO)
    end
    return block[]
end

function _mpo_mpo_trace(Os::MPO...)
    # Type info...
    T = Base.promote_op(*, eltype.(Os)...)
    conjOs = map(O->isconj(O), Os)
    transOs = map(O->istranspose(O), Os)

    # Do the contraction 
    block = cache(T, map(O->size(O[begin], 1), Os), 2, 1) .= 1
    for i in eachindex(Os[begin])
        # Contract with the first MPO in the term
        block = contract(block, Os[1][i], 1, 1, false, conjOs[1])
        if transOs[1]
            block = permutedim(block, ndims(block)-1, ndims(block)-2)
        end

        # Contract with the central MPOs
        for j in range(firstindex(Os)+1, lastindex(Os)-1)
            block = contract(block, Os[j][i], (1, ndims(block)-1),
                            (1, transOs[j] ? 3 : 2), false, conjOs[j])
        end

        # Contract with final MPO 
        block = contract(block, Os[end][i], (1, ndims(block)-1, 2),
                         (1, transOs[end] ? 3 : 2, transOs[end] ? 2 : 3), false, conjOs[end])
    end

    return block[]
end