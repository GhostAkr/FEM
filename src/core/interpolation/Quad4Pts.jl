"""
    Quad4Pts

Module describing model with bilinear quadrilateral finite elements.
"""
module Quad4Pts

using ElementTypes

export FiniteElement, Quad4Type
export jacGlobToLoc, DetJs, gradMatr, displInterpMatr, nodesFromDirection, directionFromNodes, getRSFromNode

struct Quad4Type <: FiniteElement
    name::String
end

include("../LoadVars.jl")

h1(r, s) = 0.25 * (1 + r) * (1 + s)
h2(r, s) = 0.25 * (1 - r) * (1 + s)
h3(r, s) = 0.25 * (1 - r) * (1 - s)
h4(r, s) = 0.25 * (1 + r) * (1 - s)

x(r, s, xCoords::Array{Float64}) = h1(r, s) * xCoords[1] + h2(r, s) * xCoords[2] + h3(r, s) * xCoords[3] + h4(r, s) * xCoords[4]
y(r, s, yCoords::Array{Float64}) = h1(r, s) * yCoords[1] + h2(r, s) * yCoords[2] + h3(r, s) * yCoords[3] + h4(r, s) * yCoords[4]

u(r, s, uCoords::Array{Float64}) = h1(r, s) * uCoords[1] + h2(r, s) * uCoords[2] + h3(r, s) * uCoords[3] + h4(r, s) * uCoords[4]
v(r, s, vCoords::Array{Float64}) = h1(r, s) * vCoords[1] + h2(r, s) * vCoords[2] + h3(r, s) * vCoords[3] + h4(r, s) * vCoords[4]

dxr(r, s, xCoords::Array{Float64}) = 0.25 * (1 + s) * xCoords[1] - 0.25 * (1 + s) * xCoords[2] - 0.25 * (1 - s) * xCoords[3] + 0.25 * (1 - s) * xCoords[4]
dxs(r, s, xCoords::Array{Float64}) = 0.25 * (1 + r) * xCoords[1] + 0.25 * (1 - r) * xCoords[2] - 0.25 * (1 - r) * xCoords[3] - 0.25 * (1 + r) * xCoords[4]
dyr(r, s, yCoords::Array{Float64}) = 0.25 * (1 + s) * yCoords[1] - 0.25 * (1 + s) * yCoords[2] - 0.25 * (1 - s) * yCoords[3] + 0.25 * (1 - s) * yCoords[4]
dys(r, s, yCoords::Array{Float64}) = 0.25 * (1 + r) * yCoords[1] + 0.25 * (1 - r) * yCoords[2] - 0.25 * (1 - r) * yCoords[3] - 0.25 * (1 + r) * yCoords[4]

dur(r, s, uCoords::Array{Float64}) = 0.25 * (1 + s) * uCoords[1] - 0.25 * (1 + s) * uCoords[2] - 0.25 * (1 - s) * uCoords[3] + 0.25 * (1 - s) * uCoords[4]
dus(r, s, uCoords::Array{Float64}) = 0.25 * (1 + r) * uCoords[1] + 0.25 * (1 - r) * uCoords[2] - 0.25 * (1 - r) * uCoords[3] - 0.25 * (1 + r) * uCoords[4]
dvr(r, s, vCoords::Array{Float64}) = 0.25 * (1 + s) * vCoords[1] - 0.25 * (1 + s) * vCoords[2] - 0.25 * (1 - s) * vCoords[3] + 0.25 * (1 - s) * vCoords[4]
dvs(r, s, vCoords::Array{Float64}) = 0.25 * (1 + r) * vCoords[1] + 0.25 * (1 - r) * vCoords[2] - 0.25 * (1 - r) * vCoords[3] - 0.25 * (1 + r) * vCoords[4]

dh1r(r, s) = (1 + s) / 4
dh1s(r, s) = (1 + r) / 4
dh2r(r, s) = (-1 - s) / 4
dh2s(r, s) = (1 - r) / 4
dh3r(r, s) = (-1 + s) / 4
dh3s(r, s) = (-1 + r) / 4
dh4r(r, s) = (1 - s) / 4
dh4s(r, s) = (-1 - r) / 4

"""
    jacGlobToLoc(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})

Jacobi matrix for conversion between different coordinate systems.

# Arguments
- `r`: r-coordinate;
- `s`: s-coordinate;
- `xCoords::Array{Float64}`: x coordinates of each node in current element;
- `yCoords::Array{Float64}`: y coordinates of each node in current element.
"""
ElementTypes.jacGlobToLoc(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}, elemTypeInd::Quad4Type) = [dxr(r, s, xCoords) dyr(r, s, yCoords); dxs(r, s, xCoords) dys(r, s, yCoords)]  # Jacobi's matrix for conversion from local coordinates to global

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
    elemTypeInd = Quad4Type("Quad4Type")
    inv(jacGlobToLoc(r, s, xCoords, yCoords, elemTypeInd))
end

# Supporting matrices for calculating gradient matrix
function TU(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})
    uxy = [1 + s    0   -(1 + s)    0   -(1 - s)    0   1 - s       0
           1 + r    0   1 - r       0   -(1 - r)    0   -(1 + r)    0]
    return 0.25 * jacGlobToLocInv(r, s, xCoords, yCoords) * uxy
end  # TU

function TV(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})
    vxy = [0    1 + s   0   -(1 + s)    0   -(1 - s)    0   1 - s
           0    1 + r   0   1 - r       0   -(1 - r)    0   -(1 + r)]
    return 0.25 * jacGlobToLocInv(r, s, xCoords, yCoords) * vxy
end  # TU

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
function ElementTypes.DetJs(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}, load_direction::Int, elemTypeInd::Quad4Type)
    if load_direction == 1 || load_direction == 3
        @show(sqrt(dxr(r, s, xCoords)^2 + dyr(r, s, yCoords)^2))
        return sqrt(dxr(r, s, xCoords)^2 + dyr(r, s, yCoords)^2)
    elseif load_direction == 2 || load_direction == 4
        @show(sqrt(dxs(r, s, xCoords)^2 + dys(r, s, yCoords)^2))
        return sqrt(dxs(r, s, xCoords)^2 + dys(r, s, yCoords)^2)
    else
        @error("Incorrect direction while calculating \"surface\" Jacobian")
        return nothing
    end
end


dh1x(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = jacGlobToLocInv(r, s, xCoords, yCoords)[1, 1] * dh1r(r, s) + jacGlobToLocInv(r, s, xCoords, yCoords)[1, 2] * dh1s(r, s)
dh1y(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = jacGlobToLocInv(r, s, xCoords, yCoords)[2, 1] * dh1r(r, s) + jacGlobToLocInv(r, s, xCoords, yCoords)[2, 2] * dh1s(r, s)
dh2x(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = jacGlobToLocInv(r, s, xCoords, yCoords)[1, 1] * dh2r(r, s) + jacGlobToLocInv(r, s, xCoords, yCoords)[1, 2] * dh2s(r, s)
dh2y(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = jacGlobToLocInv(r, s, xCoords, yCoords)[2, 1] * dh2r(r, s) + jacGlobToLocInv(r, s, xCoords, yCoords)[2, 2] * dh2s(r, s)
dh3x(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = jacGlobToLocInv(r, s, xCoords, yCoords)[1, 1] * dh3r(r, s) + jacGlobToLocInv(r, s, xCoords, yCoords)[1, 2] * dh3s(r, s)
dh3y(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = jacGlobToLocInv(r, s, xCoords, yCoords)[2, 1] * dh3r(r, s) + jacGlobToLocInv(r, s, xCoords, yCoords)[2, 2] * dh3s(r, s)
dh4x(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = jacGlobToLocInv(r, s, xCoords, yCoords)[1, 1] * dh4r(r, s) + jacGlobToLocInv(r, s, xCoords, yCoords)[1, 2] * dh4s(r, s)
dh4y(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = jacGlobToLocInv(r, s, xCoords, yCoords)[2, 1] * dh4r(r, s) + jacGlobToLocInv(r, s, xCoords, yCoords)[2, 2] * dh4s(r, s)

"""
    gradMatr(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})

Gradient matrix ``B`` for element.

# Arguments
- `r`: r-coordinate;
- `s`: s-coordinate;
- `xCoords::Array{Float64}`: x coordinates of each node in current element;
- `yCoords::Array{Float64}`: y coordinates of each node in current element.
"""
function ElementTypes.gradMatr(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}, elemTypeInd::Quad4Type)
    resultMatrix = Array{Float64, 2}(undef, 3, 8)
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
- `s`: s-coordinate;
- `elemTypeInd::Quad4Type`: type of finite element indicator.
"""
function ElementTypes.displInterpMatr(r, s, elemTypeInd::Quad4Type)
    result = [  h1(r, s)   0   h2(r, s)    0   h3(r, s)    0    h4(r, s)    0
                0       h1(r, s)    0   h2(r, s)    0   h3(r, s)    0       h4(r, s)]

    return result
end  # displInterpMatr

"""
    nodesFromDirection(direction::Int)

Return nodes of element associated with given load direction.

# Arguments
- `direction::Int`: given load direction.
"""
function ElementTypes.nodesFromDirection(direction::Int, elemTypeInd::Quad4Type)
    if direction == 1
        return [1, 2]
    elseif direction == 2
        return [2, 3]
    elseif direction == 3
        return [3, 4]
    elseif direction == 4
        return [1, 4]
    else
        println("Given load direction is not supported")
        return nothing
    end
end  # nodesFromDirection

"""
    directionFromNodes(nodes::Array, elType:: FiniteElement)

Return direction from given nodes.

# Arguments
- `nodes::Array`: nodes according to which direction should be defined.
- `elType:: FiniteElement`: element type indicator.
"""
function ElementTypes.directionFromNodes(nodes::Array, elemTypeInd::Quad4Type)
    if issubset([1, 2], nodes)
        return 1  # Top
    elseif issubset([3, 2], nodes)  # TODO: provide order-insensetive way to define direction
        return 2  # Left            # Now it only works for [3, 2] case.
    elseif issubset([3, 4], nodes)
        return 3  # Bottom
    elseif issubset([1, 4], nodes)
        return 4  # Right
    else
        @error("Can't define direction from given nodes")
        return nothing
    end
end  # directionFromNodes

function ElementTypes.getRSFromNode(nodeIndex::Int, elemTypeInd::Quad4Type)
    if nodeIndex == 1
        return (1, 1)
    elseif nodeIndex == 2
        return (-1, 1)
    elseif nodeIndex == 3
        return (-1, -1)
    elseif nodeIndex == 4
        return (1, -1)
    else
        println("Invalid node index while getting (r, s) coordinates from node")
        return nothing
    end
end  # getRSFromNode

# interFunc = [h1, h2, h3, h4]  # Array of interpolation functions

end  # Quad4Pts