"""
    VectorFEM
Module describing math vectors in FEM. Vectors in FEM are represented as simple tuples with
3 elements (x, y and z coordinates of that vector). This module contains functions which
represent common operations with such vectors.
"""
module VectorFEM

export veclength_3d, vecfrompoints_2d, vecfrompoints_3d

"""
    vecfrompoints_2d(firstpt::Tuple{Number, Number}, 
        secondpt::Tuple{Number, Number, Number}

Construct vector in form `Tuple{Number, Number, Number}` by 2 points. Math which is used
to perform operation: ``p_x = p_2_x - p_1_x``, ``p_y = p_2_y - p_1_y``, ``p_y = 0``.

# Arguments
- `firstpt::Tuple{Number, Number}`: first point;
- `secondpt::Tuple{Number, Number}`: second point
"""
function vecfrompoints_2d(firstpt::Tuple{Number, Number}, 
    secondpt::Tuple{Number, Number}
)
    return (secondpt[1] - firstpt[1], secondpt[2] - firstpt[2], 0)
end

"""
    vecfrompoints_3d(firstpt::Tuple{Number, Number, Number}, 
        secondpt::Tuple{Number, Number, Number}

Construct vector in form `Tuple{Number, Number, Number}` by 2 points. Math which is used
to perform operation: ``p_x = p_2_x - p_1_x``, ``p_y = p_2_y - p_1_y``, ``p_y = p_2_y - 
p_1_z``.

# Arguments
- `firstpt::Tuple{Number, Number, Number}`: first point;
- `secondpt::Tuple{Number, Number, Number}`: second point
"""
function vecfrompoints_3d(firstpt::Tuple{Number, Number, Number}, 
    secondpt::Tuple{Number, Number, Number}
)
    return (secondpt[1] - firstpt[1], secondpt[2] - firstpt[2], secondpt[3] - firstpt[3])
end

"""
    veclength(vec::Tuple{Number, Number, Number})

Compute euclidian length of vector `vec`.

# Arguments
- `vec::Tuple{Number, Number, Number}`: vector which length needs to be calculated.
"""
function veclength_3d(vec::Tuple{Number, Number, Number})
    return sqrt(vec[1]^2 + vec[2]^2 + vec[3]^2)
end

end  # VectorFEM