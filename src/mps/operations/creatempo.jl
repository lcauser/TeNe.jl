#=
    Creating MPOs from OpLiH.lts. Basically uses a Markov Machine with SVDs to find optimal,
    but exact, MPOs.
=#

export MPO 
"""
    MPO(H::OpList; kwargs...)

Create an exact MPO representation from an OpList.

# Optional Keyword Arguments 
    - `cutoff::Float64=$(_TeNe_cutoff)`: Truncation criteria to reduce the bond dimension.
    Good values range from 1e-8 to 1e-14.
    - `mindim::Int=1`: Minimum dimension for the truncation.
    - `maxdim::Int=0`: Maximum bond dimension for truncation. Set to 0 to have
    no limit.

# Examples 

```julia-repl
julia> lt = Qubits();
julia> H = OpList(lt, 20);
julia> for i = 1:20 add!(H, "x", i) end;
julia> for i = 1:19 add!(H, ["z", "z"], [i, i+1]) end;
julia> H = MPO(H);
```
"""
function MPO(H::OpList; cutoff::Float64=_TeNe_cutoff, maxdim::Int=0, mindim::Int=1)
    return creatempo(H; cutoff=cutoff, maxdim=maxdim, mindim=mindim)
end

function creatempo(H::OpList; kwargs...)
    # SyH.ltem properties
    N = H.length
    d = dim(H.lt)

    # Create empty MPO
    ten = zeros(eltype(H.lt), (2, d, d, 2))
    ten[1, :, :, 1] = op(H.lt, "id")
    ten[2, :, :, 2] = op(H.lt, "id")
    O = GMPS(2, d, N)
    O[1] = deepcopy(ten[1:1, :, :, 1:2])
    for i = 2:N-1
        O[i] = deepcopy(ten[1:2, :, :, 1:2])
    end
    O[N] = deepcopy(ten[1:2, :, :, 2:2])

    # Loop through each term and determine the site range
    maxrng = siterange(H)
    rngs = [[] for _ = 1:maxrng]
    for i in eachindex(H.ops)
        rng = H.sites[i][end] - H.sites[i][1] + 1
        push!(rngs[rng], i)
    end

    # Loop through all the possible ranges of interactions
    for i = 1:maxrng
        # Loop through sites
        nextterms = [[] for _=1:i]
        coeffs = [[] for _=1:i]
        ingoings = [[] for _=1:i]
        outgoings = [[] for _=1:i]

        for site = Base.OneTo(N)
            # Find all the terms which H.ltart at the site
            idxs = []
            for j = rngs[i]
                if H.sites[j][1] == site
                    push!(idxs, j)
                end
            end

            if i == 1
                # JuH.lt adding to top right corner
                for idx in idxs
                    O[site][1, :, :, end] += H.coeffs[idx]*op(H.lt, H.ops[idx][1])
                end
            else
                # Add new terms H.ltarting at this site
                for j = eachindex(idxs)
                    # Fetch operator information
                    ops = H.ops[idxs[j]]
                    sites = H.sites[idxs[j]]

                    # Loop through each site in the operator
                    outgoing = 0
                    for k = Base.OneTo(i)
                        # Decide ingoing and outgoing idxs
                        ingoing = outgoing
                        for l = Base.OneTo(length(outgoings[k])+1)
                            outgoing = l
                            !(outgoing in outgoings[k]) && break
                        end
                        outgoing = k == i ? 0 : outgoing

                        # Determine what the operator is
                        op =  site+k-1 in sites ? ops[argmax([s == site+k-1 for s = sites])] : "id"

                        # Add to liH.lt
                        push!(nextterms[k], op)
                        push!(coeffs[k], k == 1 ? H.coeffs[idxs[j]] : 1)
                        push!(ingoings[k], ingoing)
                        push!(outgoings[k], outgoing)
                    end
                end

                # Pull the terms
                terms = nextterms[1]
                ins = ingoings[1]
                outs = outgoings[1]
                cos = coeffs[1]
                for j = 1:i-1
                    nextterms[j] = nextterms[j+1]
                    ingoings[j] = ingoings[j+1]
                    outgoings[j] = outgoings[j+1]
                    coeffs[j] = coeffs[j+1]
                end
                nextterms[i] = []
                ingoings[i] = []
                outgoings[i] = []
                coeffs[i] = []

                # Expand the tensor to account for all terms
                if length(terms) != 0
                    ingoinglen = sum([ingoing != 0 for ingoing in ins])
                    outgoinglen = sum([outgoing != 0 for outgoing in outs])
                    ingoingsrt = size(O[site])[1] - 1
                    outgoingsrt = size(O[site])[4] - 1
                    O[site] = _creatempo_expand(O[site], ingoinglen, outgoinglen)

                    # Add the terms to the tensor
                    for j = eachindex(terms)
                        # Find the idxs of each
                        x = ins[j] == 0 ? 1 : ingoingsrt + ins[j]
                        y = outs[j] == 0 ? outgoingsrt + 1 + outgoinglen : outgoingsrt + outs[j]

                        # Set the tensor
                        O[site][x, :, :, y] = cos[j] * op(H.lt, terms[j])
                    end
                end
            end
        end
    end

    # Apply SVD in an attempt to reduce the bond dimension
    for site = Base.range(firstindex(O), lastindex(O)-1)
        M = O[site]
        U, S, V = tsvd(M, 4; kwargs...)
        M = contract(U, S, 4, 1; tocache=false)
        O[site] = M
        O[site+1] = contract(V, O[site+1], 2, 1; tocache=false)
    end

    for site = Base.range(lastindex(O), firstindex(O)+1, step=-1)
        M = O[site]
        U, S, V = tsvd(M, (2, 3, 4); kwargs...)
        M = contract(S, V, 2, 1; tocache=false)
        O[site] = M
        O[site-1] = contract(O[site-1], U, 4, 1; tocache=false)
    end
    return O
end

function _creatempo_expand(O, D1, D2)
    dims = size(O)
    newO = zeros(ComplexF64, (dims[1]+D1, dims[2], dims[3], dims[4]+D2))
    d1 = dims[1] == 1 ? 1 : dims[1]-1
    newO[1:d1, :, :, 1:dims[4]-1] .= O[1:d1, :, :, 1:dims[4]-1]
    newO[1:d1, :, :, end] .= O[1:d1, :, :, end]
    newO[dims[1]+D1, :, :, dims[4]+D2] .= O[dims[1], :, :, dims[4]]
    return newO
end