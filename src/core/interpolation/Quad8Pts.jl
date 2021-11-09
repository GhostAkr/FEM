module Quad8Pts

using ElementTypes

export FiniteElement, Quad8Type
export jacGlobToLoc, DetJs, gradMatr, displInterpMatr, nodesFromDirection, getRSFromNode

struct Quad8Type <: FiniteElement
    name::String
end

h1(r, s) = 0.25 * (1 + r) * (1 + s) - 0.5 * h5(r, s) - 0.5 * h8(r, s)
h2(r, s) = 0.25 * (1 - r) * (1 + s) - 0.5 * h5(r, s) - 0.5 * h6(r, s)
h3(r, s) = 0.25 * (1 - r) * (1 - s) - 0.5 * h6(r, s) - 0.5 * h7(r, s)
h4(r, s) = 0.25 * (1 + r) * (1 - s) - 0.5 * h7(r, s) - 0.5 * h8(r, s)
h5(r, s) = 0.5 * (1 - r^2) * (1 + s)
h6(r, s) = 0.5 * (1 - s^2) * (1 - r)
h7(r, s) = 0.5 * (1 - r^2) * (1 - s)
h8(r, s) = 0.5 * (1 - s^2) * (1 + r)

function dxr(r, s, xCoords::Array{Float64})
    c1 = (1 + s) / 4 + 0.5 * r * (1 + s) + 0.25 * (-1 + s^2)
    c2 = 0.25 * (-1 - s) + 0.5 * r * (1 + s) + 0.25 * (1 - s^2)
    c3 = (0.5 * r * (1 - s) + 0.25 * (-1 + s) + 0.25 * (1 - s^2))
    c4 = (1 - s) / 4 + 0.5 * r * (1 - s) + 0.25 * (-1 + s^2)
    c5 = -r * (1 + s)
    c6 = -0.5 * (1 - s^2)
    c7 = -r * (1 - s)
    c8 = 0.5 * (1 - s^2)
    return c1 * xCoords[1] + c2 * xCoords[2] + c3 * xCoords[3] + c4 * xCoords[4] + c5 * xCoords[5] + c6 * xCoords[6] + c7 * xCoords[7] + c8 * xCoords[8]
end  # dxr

function dyr(r, s, yCoords::Array{Float64})
    c1 = (1 + s) / 4 + 0.5 * r * (1 + s) + 0.25 * (-1 + s^2)
    c2 = 0.25 * (-1 - s) + 0.5 * r * (1 + s) + 0.25 * (1 - s^2)
    c3 = (0.5 * r * (1 - s) + 0.25 * (-1 + s) + 0.25 * (1 - s^2))
    c4 = (1 - s) / 4 + 0.5 * r * (1 - s) + 0.25 * (-1 + s^2)
    c5 = -r * (1 + s)
    c6 = -0.5 * (1 - s^2)
    c7 = -r * (1 - s)
    c8 = 0.5 * (1 - s^2)
    return c1 * yCoords[1] + c2 * yCoords[2] + c3 * yCoords[3] + c4 * yCoords[4] + c5 * yCoords[5] + c6 * yCoords[6] + c7 * yCoords[7] + c8 * yCoords[8]
end  # dyr

function dxs(r, s, xCoords::Array{Float64})
    c1 = (1 + r) / 4 + 0.25 * (-1 + r^2) + 0.5 * (1 + r) * s
    c2 = (1 - r) / 4 + 0.25 * (-1 + r^2) + 0.5 * (1 - r) * s
    c3 = 0.25 * (-1 + r) + 0.25 * (1 - r^2) + 0.5 * (1 - r) * s
    c4 = 0.25 * (-1 - r) + 0.25 * (1 - r^2) + 0.5 * (1 + r) * s
    c5 = 0.5 * (1 - r^2)
    c6 = -(1 - r) * s
    c7 = -0.5 * (1 - r^2)
    c8 = -(1 + r) * s
    return c1 * xCoords[1] + c2 * xCoords[2] + c3 * xCoords[3] + c4 * xCoords[4] + c5 * xCoords[5] + c6 * xCoords[6] + c7 * xCoords[7] + c8 * xCoords[8]
end  # dxs

function dys(r, s, yCoords::Array{Float64})
    c1 = (1 + r) / 4 + 0.25 * (-1 + r^2) + 0.5 * (1 + r) * s
    c2 = (1 - r) / 4 + 0.25 * (-1 + r^2) + 0.5 * (1 - r) * s
    c3 = 0.25 * (-1 + r) + 0.25 * (1 - r^2) + 0.5 * (1 - r) * s
    c4 = 0.25 * (-1 - r) + 0.25 * (1 - r^2) + 0.5 * (1 + r) * s
    c5 = 0.5 * (1 - r^2)
    c6 = -(1 - r) * s
    c7 = -0.5 * (1 - r^2)
    c8 = -(1 + r) * s
    return c1 * yCoords[1] + c2 * yCoords[2] + c3 * yCoords[3] + c4 * yCoords[4] + c5 * yCoords[5] + c6 * yCoords[6] + c7 * yCoords[7] + c8 * yCoords[8]
end  # dys

    """
    jacGlobToLoc(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})

Jacobi matrix for conversion between different coordinate systems.

# Arguments
- `r`: r-coordinate;
- `s`: s-coordinate;
- `xCoords::Array{Float64}`: x coordinates of each node in current element;
- `yCoords::Array{Float64}`: y coordinates of each node in current element.
"""
ElementTypes.jacGlobToLoc(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}, elemTypeInd::Quad8Type) = [dxr(r, s, xCoords) dyr(r, s, yCoords); dxs(r, s, xCoords) dys(r, s, yCoords)]  # Jacobi's matrix for conversion from local coordinates to global

"""
    jacGlobToLocInv(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})

Inversed Jacobi matrix for conversion between different coordinate systems.

# Arguments
- `r`: r-coordinate;
- `s`: s-coordinate;
- `xCoords::Array{Float64}`: x coordinates of each node in current element;
- `yCoords::Array{Float64}`: y coordinates of each node in current element.
"""
function jacGlobToLocInv(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})
    elemTypeInd = Quad8Type("Quad8Type")
    inv(jacGlobToLoc(r, s, xCoords, yCoords, elemTypeInd))
end

"""
    DetJs(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})

Determinant of appropriate Jacobi matrix.

# Arguments
- `r`: r-coordinate;
- `s`: s-coordinate;
- `xCoords::Array{Float64}`: x coordinates of each node in current element;
- `yCoords::Array{Float64}`: y coordinates of each node in current element;
- `load_direction::Int`: direction of load;
- `elemTypeInd::Quad4Type`: type of finite element indicator.
"""
function ElementTypes.DetJs(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}, load_direction::Int, elemTypeInd::Quad8Type)
    if load_direction == 1 || load_direction == 3
        return sqrt(dxr(r, s, xCoords)^2 + dyr(r, s, yCoords)^2)
    elseif load_direction == 2 || load_direction == 4
        return sqrt(dxs(r, s, xCoords)^2 + dys(r, s, yCoords)^2)
    else
        @error("Incorrect direction while calculating \"surface\" Jacobian")
        return nothing
    end
end  # DetJs

# Supporting matrices for calculating gradient matrix
function TU(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})
    # Coefficients for dur (rc1 * u1 + rc2 * u2 + ... + rc8 * u8)
    rc1 = (1 + s) / 4 + 0.5 * r * (1 + s) + 0.25 * (-1 + s^2)
    rc2 = 0.25 * (-1 - s) + 0.5 * r * (1 + s) + 0.25 * (1 - s^2)
    rc3 = (0.5 * r * (1 - s) + 0.25 * (-1 + s) + 0.25 * (1 - s^2))
    rc4 = (1 - s) / 4 + 0.5 * r * (1 - s) + 0.25 * (-1 + s^2)
    rc5 = -r * (1 + s)
    rc6 = -0.5 * (1 - s^2)
    rc7 = -r * (1 - s)
    rc8 = 0.5 * (1 - s^2)
    # Coefficients for dus (sc1 * u1 + sc2 * u2 + ... + sc8 * u8)
    sc1 = (1 + r) / 4 + 0.25 * (-1 + r^2) + 0.5 * (1 + r) * s
    sc2 = (1 - r) / 4 + 0.25 * (-1 + r^2) + 0.5 * (1 - r) * s
    sc3 = 0.25 * (-1 + r) + 0.25 * (1 - r^2) + 0.5 * (1 - r) * s
    sc4 = 0.25 * (-1 - r) + 0.25 * (1 - r^2) + 0.5 * (1 + r) * s
    sc5 = 0.5 * (1 - r^2)
    sc6 = -(1 - r) * s
    sc7 = -0.5 * (1 - r^2)
    sc8 = -(1 + r) * s

    uxy = [rc1 0 rc2 0 rc3 0 rc4 0 rc5 0 rc6 0 rc7 0 rc8 0
           sc1 0 sc2 0 sc3 0 sc4 0 sc5 0 sc6 0 sc7 0 sc8 0]
    return jacGlobToLocInv(r, s, xCoords, yCoords) * uxy
end  # TU

function TV(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})
    # Coefficients for dvr (rc1 * v1 + rc2 * v2 + ... + rc8 * v8)
    rc1 = (1 + s) / 4 + 0.5 * r * (1 + s) + 0.25 * (-1 + s^2)
    rc2 = 0.25 * (-1 - s) + 0.5 * r * (1 + s) + 0.25 * (1 - s^2)
    rc3 = (0.5 * r * (1 - s) + 0.25 * (-1 + s) + 0.25 * (1 - s^2))
    rc4 = (1 - s) / 4 + 0.5 * r * (1 - s) + 0.25 * (-1 + s^2)
    rc5 = -r * (1 + s)
    rc6 = -0.5 * (1 - s^2)
    rc7 = -r * (1 - s)
    rc8 = 0.5 * (1 - s^2)
    # Coefficients for dvs (sc1 * v1 + sc2 * v2 + ... + sc8 * v8)
    sc1 = (1 + r) / 4 + 0.25 * (-1 + r^2) + 0.5 * (1 + r) * s
    sc2 = (1 - r) / 4 + 0.25 * (-1 + r^2) + 0.5 * (1 - r) * s
    sc3 = 0.25 * (-1 + r) + 0.25 * (1 - r^2) + 0.5 * (1 - r) * s
    sc4 = 0.25 * (-1 - r) + 0.25 * (1 - r^2) + 0.5 * (1 + r) * s
    sc5 = 0.5 * (1 - r^2)
    sc6 = -(1 - r) * s
    sc7 = -0.5 * (1 - r^2)
    sc8 = -(1 + r) * s
    vxy = [0 rc1 0 rc2 0 rc3 0 rc4 0 rc5 0 rc6 0 rc7 0 rc8
           0 sc1 0 sc2 0 sc3 0 sc4 0 sc5 0 sc6 0 sc7 0 sc8]
    return jacGlobToLocInv(r, s, xCoords, yCoords) * vxy
end  # TU

"""
    gradMatr(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})

Gradient matrix ``B`` for element.

# Arguments
- `r`: r-coordinate;
- `s`: s-coordinate;
- `xCoords::Array{Float64}`: x coordinates of each node in current element;
- `yCoords::Array{Float64}`: y coordinates of each node in current element.
"""
function ElementTypes.gradMatr(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}, elemTypeInd::Quad8Type)
    if size(xCoords)[1] != size(yCoords)[1]
        println("Incorrect nodes while calculating gradient matrix")
        return nothing
    end
    nodesPerElement = size(xCoords)[1]
    resultMatrix = Array{Float64, 2}(undef, 3, nodesPerElement * 2)
    Uxy = TU(r, s, xCoords, yCoords)
    Vxy = TV(r, s, xCoords, yCoords)
    if (size(Uxy)[2] != size(Vxy)[2] != size(resultMatrix)[2])
        println("Error while calculating gradient matrix")
        return resultMatrix
    end
    cols = size(Uxy)[2]
    for i in 1:cols
        resultMatrix[1, i] = Uxy[1, i]
        resultMatrix[2, i] = Vxy[2, i]
        resultMatrix[3, i] = Uxy[2, i] + Vxy[1, i]
    end
    return resultMatrix
end

"""
    displInterpMatr(r, s)

Displacements interpolation ``H`` matrix.

# Arguments
- `r`: r-coordinate;
- `s`: s-coordinate.
"""
function ElementTypes.displInterpMatr(r, s, elemTypeInd::Quad8Type)
    result = [h1(r, s) 0 h2(r, s) 0 h3(r, s) 0 h4(r, s) 0 h5(r, s) 0 h6(r, s) 0 h7(r, s) 0 h8(r, s) 0
              0 h1(r, s) 0 h2(r, s) 0 h3(r, s) 0 h4(r, s) 0 h5(r, s) 0 h6(r, s) 0 h7(r, s) 0 h8(r, s)]
    return result
end  # displInterpMatr

"""
    nodesFromDirection(direction::Int)

Return nodes of element associated with given load direction.

# Arguments
- `direction::Int`: given load direction.
"""
function ElementTypes.nodesFromDirection(direction::Int, elemTypeInd::Quad8Type)
    if direction == 1  # Top
        return [1, 5, 2]
    elseif direction == 2  # Left
        return [2, 6, 3]
    elseif direction == 3  # Bottom
        return [3, 7, 4]
    elseif direction == 4  # Right
        return [1, 8, 4]
    else
        println("Given load direction is not supported")
        return nothing
    end
end  # nodesFromDirection

"""
    directionFromNodes(nodes::Array, elemTypeInd::Quad8Type)

Return direction from given nodes.

# Arguments
- `nodes::Array`: nodes according to which direction should be defined.
- `elemTypeInd::Quad8Type`: element type indicator.
"""
function ElementTypes.directionFromNodes(nodes::Array, elemTypeInd::Quad8Type)
    if issubset([1, 5, 2], nodes)
        return 1  # Top
    elseif issubset([2, 6, 3], nodes)
        return 2  # Left
    elseif issubset([3, 7, 4], nodes)
        return 3  # Bottom
    elseif issubset([1, 8, 4], nodes)
        return 4  # Right
    else
        @error("Can't define direction from given nodes")
        return nothing
    end
end  # directionFromNodes

function ElementTypes.getRSFromNode(nodeIndex::Int, elemTypeInd::Quad8Type)
    if nodeIndex == 1
        return (1, 1)
    elseif nodeIndex == 2
        return (-1, 1)
    elseif nodeIndex == 3
        return (-1, -1)
    elseif nodeIndex == 4
        return (1, -1)
    elseif nodeIndex == 5
        return (0, 1)
    elseif nodeIndex == 6
        return (-1, 0)
    elseif nodeIndex == 7
        return (0, -1)
    elseif nodeIndex == 8
        return (1, 0)
    else
        println("Invalid node index while getting (r, s) coordinates from node")
        return nothing
    end
end  # getRSFromNode

end  # Quad8Pts