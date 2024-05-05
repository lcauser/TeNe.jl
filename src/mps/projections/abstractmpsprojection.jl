abstract type MPSProjection end 

### Properties of projections
export center, leftblock, rightblock, block
"""
    length(proj::MPSProjection)

Determine the number of sites in an MPS projection.
"""
Base.length(proj::MPSProjection) = length(proj.objects[1])


"""
    center(proj::MPSProjection)

Return the center of the block.
"""
center(proj::MPSProjection) = proj.center


### Retrieving blocks from a projection 
"""
    leftedge(proj::MPSProjection)

Return the left edge block for a projection.
"""
function leftedge(proj::MPSProjection)
    return cache_ones(Tuple(map(j->size(proj.objects[j][1], 1), Base.eachindex(proj.objects))))
end


"""
    rightedge(proj::MPSProjection)

Return the right edge block for a projection.
"""
function rightedge(proj::MPSProjection)
    return cache_ones(Tuple(map(j->size(proj.objects[j][end], ndims(proj.objects[j][end])),
        Base.eachindex(proj.objects))))
end

"""
    leftblock(proj::MPSProjection, idx::Int)

Fetch the left block at a given site.
"""
function leftblock(proj::MPSProjection, idx::Int)
    if idx < 1
        return leftedge(proj)
    elseif idx >= length(proj)
        throw(ArgumentError("Index $(idx) is out of range."))
    else
        return proj.lefts[idx]
    end
end

"""
    rightblock(proj::MPSProjection, idx::Int)

Fetch the right block at a given site.
"""
function rightblock(proj::MPSProjection, idx::Int)
    if idx > length(proj)
        return rightedge(proj)
    elseif idx <= 1
        throw(ArgumentError("Index $(idx) is out of range."))
    else
        return proj.rights[idx-1]
    end
end

"""
    block(proj::MPSProjection, idx::Int)

Fetch the block at a given site.
"""
function block(proj::MPSProjection, idx::Int)
    if idx > proj.center
        return rightblock(proj, idx)
    elseif idx < proj.center 
        return leftblock(proj, idx)
    else
        throw(ArgumentError("There is ambiguity in indexing at the centre of the projection. Use `leftblock` or `rightblock`."))
    end
end

Base.getindex(proj::MPSProjection, idx::Int) = block(proj, idx)
function Base.setindex!(proj::MPSProjection, x, idx::Int)
    proj.blocks[idx] = x
    return proj
end

### Moving the centre of a block 
export movecenter!
"""
    movecenter!(proj::MPSProjection, idx::Int)

Move the center of the projection.
"""
function movecenter!(proj::MPSProjection, idx::Int)
    if center(proj) == 0
        for i in Base.OneTo(idx-1)
            buildleft!(proj, i)
        end
        for i in Base.range(length(proj), idx+1, step=-1)
            buildright!(proj, i)
        end
    else
        if idx > center(proj)
            for i in Base.range(center(proj), idx-1)
                buildleft!(proj, i)
            end
        elseif idx < center(proj)
            for i = Base.range(center(proj), idx+1, step=-1)
                buildright!(proj, i)
            end
        end
    end
    proj.center = idx
end