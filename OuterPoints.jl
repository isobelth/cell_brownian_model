
module OuterPoints

using StaticArrays

@inline function outerPoints!(pos::MMatrix)
    hull = sp.ConvexHull(reshape(filter(!iszero, pos), (Ncells,3)))
    vertex_points = [hull.vertices.+1]
    return vertex_points
end

export outerPoints

end
