#=
    Functions for contracting tensors 
=#



function contract(x, y, cix, ciy, conjx=false, conjy=false)

end

function contract!(z, x, y, cix, ciy, conjx=false, conjy=false)
    # Retrieve the dimensions of x, y, z
    ix = 1:ndims(x)
    iy = 1:ndims(y)
    sx = size(x)
    sy = size(y)
    sz = size(z)
    cix = collect(cix)
    ciy = collect(ciy)

    # Permuted indices
    rix = collect(symdiff(ix, cix))
    riy = collect(symdiff(iy, ciy))
    pix = [rix; cix]
    piy = [ciy; riy]
    
    # Checks to ensure contraction can be done 
    sx[cix] != sy[ciy] && throw(ArgumentError("Dimensions $(sx[cix]) do not match $(sy[ciy])."))
    collect(sz) != [sx[rix]...; sy[riy]...] && throw(ArgumentError("Invalid dimensions in the resulting array."))

    # Permute the tensors 
    px = permutedims!(cache(eltype(x), sx[pix]), x, pix)
    py = permutedims!(cache(eltype(y), sy[piy], length(x) == length(y) ? 2 : 1), y, piy)

    # Conjugations 
    ifelse(conjx, px .= conj.(px), nothing)
    ifelse(conjy, py .= conj.(py), nothing)

    # Reshape the tensors
    px = reshape(px, (prod(sx[cix]), prod(sx[rix])))
    py = reshape(py, (prod(sy[ciy]), prod(sy[riy])))

    # Contract
    mul!(reshape(z, (prod(sx[rix]), prod(sy[riy]))),
         reshape(px, (prod(sx[rix]), prod(sx[cix]))),
         reshape(py, (prod(sy[ciy]), prod(sy[riy]))))
end