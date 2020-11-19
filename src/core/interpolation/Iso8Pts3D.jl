"""
    Iso8Pts3D

Module describing isoparametric 8-points 3D finite elements.
"""
module Iso8Pts3D

using ElementTypes

export FiniteElement, Iso8Pts3DType


struct Iso8Pts3DType <: FiniteElement
    name::String
end

include("../LoadVars.jl")

# Interpolation functions

h1(r, s, t) = 0.125 * (1 + r) * (1 + s) * (1 + t)
h2(r, s, t) = 0.125 * (1 - r) * (1 + s) * (1 + t)
h3(r, s, t) = 0.125 * (1 - r) * (1 - s) * (1 + t)
h4(r, s, t) = 0.125 * (1 + r) * (1 - s) * (1 + t)
h5(r, s, t) = 0.125 * (1 + r) * (1 + s) * (1 - t)
h6(r, s, t) = 0.125 * (1 - r) * (1 + s) * (1 - t)
h7(r, s, t) = 0.125 * (1 - r) * (1 - s) * (1 - t)
h8(r, s, t) = 0.125 * (1 + r) * (1 - s) * (1 - t)

# Coordinates derivatives

function dxr(r, s, t, xCoords::Array{Float64})
    # c_i --- coefficients of coordinate derivative
    c1 = 0.125 * (1 + s) * (1 + t)
    c2 = -0.125 * (1 + s) * (1 + t)
    c3 = -0.125 * (1 - s) * (1 + t)
    c4 = 0.125 * (1 - s) * (1 + t)
    c5 = 0.125 * (1 + s) * (1 - t)
    c6 = -0.125 * (1 + s) * (1 - t)
    c7 = -0.125 * (1 - s) * (1 - t)
    c8 = 0.125 * (1 - s) * (1 - t)
    return c1 * xCoords[1] + c2 * xCoords[2] + c3 * xCoords[3] + c4 * xCoords[4] + c5 * xCoords[5] + c6 * xCoords[6] + c7 * xCoords[7] + c8 * xCoords[8]
end

function dyr(r, s, t, yCoords::Array{Float64})
    # c_i --- coefficients of coordinate derivative
    c1 = 0.125 * (1 + s) * (1 + t)
    c2 = -0.125 * (1 + s) * (1 + t)
    c3 = -0.125 * (1 - s) * (1 + t)
    c4 = 0.125 * (1 - s) * (1 + t)
    c5 = 0.125 * (1 + s) * (1 - t)
    c6 = -0.125 * (1 + s) * (1 - t)
    c7 = -0.125 * (1 - s) * (1 - t)
    c8 = 0.125 * (1 - s) * (1 - t)
    return c1 * yCoords[1] + c2 * yCoords[2] + c3 * yCoords[3] + c4 * yCoords[4] + c5 * yCoords[5] + c6 * yCoords[6] + c7 * yCoords[7] + c8 * yCoords[8]
end

function dzr(r, s, t, zCoords::Array{Float64})
    # c_i --- coefficients of coordinate derivative
    c1 = 0.125 * (1 + s) * (1 + t)
    c2 = -0.125 * (1 + s) * (1 + t)
    c3 = -0.125 * (1 - s) * (1 + t)
    c4 = 0.125 * (1 - s) * (1 + t)
    c5 = 0.125 * (1 + s) * (1 - t)
    c6 = -0.125 * (1 + s) * (1 - t)
    c7 = -0.125 * (1 - s) * (1 - t)
    c8 = 0.125 * (1 - s) * (1 - t)
    return c1 * zCoords[1] + c2 * zCoords[2] + c3 * zCoords[3] + c4 * zCoords[4] + c5 * zCoords[5] + c6 * zCoords[6] + c7 * zCoords[7] + c8 * zCoords[8]
end

function dxs(r, s, t, xCoords::Array{Float64})
    # c_i --- coefficients of coordinate derivative
    c1 = 0.125 * (1 + r) * (1 + t)
    c2 = 0.125 * (1 - r) * (1 + t) 
    c3 = -0.125 * (1 - r) * (1 + t)
    c4 = -0.125 * (1 + r) * (1 + t) 
    c5 = 0.125 * (1 + r) * (1 - t)
    c6 = 0.125 * (1 - r) * (1 - t)
    c7 = -0.125 * (1 - r) * (1 - t)
    c8 = -0.125 * (1 + r) * (1 - t) 
    return c1 * xCoords[1] + c2 * xCoords[2] + c3 * xCoords[3] + c4 * xCoords[4] + c5 * xCoords[5] + c6 * xCoords[6] + c7 * xCoords[7] + c8 * xCoords[8]
end

function dys(r, s, t, yCoords::Array{Float64})
    # c_i --- coefficients of coordinate derivative
    c1 = 0.125 * (1 + r) * (1 + t)
    c2 = 0.125 * (1 - r) * (1 + t) 
    c3 = -0.125 * (1 - r) * (1 + t)
    c4 = -0.125 * (1 + r) * (1 + t) 
    c5 = 0.125 * (1 + r) * (1 - t)
    c6 = 0.125 * (1 - r) * (1 - t)
    c7 = -0.125 * (1 - r) * (1 - t)
    c8 = -0.125 * (1 + r) * (1 - t) 
    return c1 * yCoords[1] + c2 * yCoords[2] + c3 * yCoords[3] + c4 * yCoords[4] + c5 * yCoords[5] + c6 * yCoords[6] + c7 * yCoords[7] + c8 * yCoords[8]
end

function dzs(r, s, t, zCoords::Array{Float64})
    # c_i --- coefficients of coordinate derivative
    c1 = 0.125 * (1 + r) * (1 + t)
    c2 = 0.125 * (1 - r) * (1 + t) 
    c3 = -0.125 * (1 - r) * (1 + t)
    c4 = -0.125 * (1 + r) * (1 + t) 
    c5 = 0.125 * (1 + r) * (1 - t)
    c6 = 0.125 * (1 - r) * (1 - t)
    c7 = -0.125 * (1 - r) * (1 - t)
    c8 = -0.125 * (1 + r) * (1 - t) 
    return c1 * zCoords[1] + c2 * zCoords[2] + c3 * zCoords[3] + c4 * zCoords[4] + c5 * zCoords[5] + c6 * zCoords[6] + c7 * zCoords[7] + c8 * zCoords[8]
end

function dxt(r, s, t, xCoords::Array{Float64})
    # c_i --- coefficients of coordinate derivative
    c1 = 0.125 * (1 + r) * (1 + s)
    c2 = 0.125 * (1 - r) * (1 + s)
    c3 = 0.125 * (1 - r) * (1 - s)
    c4 = 0.125 * (1 + r) * (1 - s)
    c5 = -0.125 * (1 + r) * (1 + s) 
    c6 = -0.125 * (1 - r) * (1 + s)
    c7 = -0.125 * (1 - r) * (1 - s)
    c8 = -0.125 * (1 + r) * (1 - s)
    return c1 * xCoords[1] + c2 * xCoords[2] + c3 * xCoords[3] + c4 * xCoords[4] + c5 * xCoords[5] + c6 * xCoords[6] + c7 * xCoords[7] + c8 * xCoords[8]
end

function dyt(r, s, t, yCoords::Array{Float64})
    # c_i --- coefficients of coordinate derivative
    c1 = 0.125 * (1 + r) * (1 + s)
    c2 = 0.125 * (1 - r) * (1 + s)
    c3 = 0.125 * (1 - r) * (1 - s)
    c4 = 0.125 * (1 + r) * (1 - s)
    c5 = -0.125 * (1 + r) * (1 + s) 
    c6 = -0.125 * (1 - r) * (1 + s)
    c7 = -0.125 * (1 - r) * (1 - s)
    c8 = -0.125 * (1 + r) * (1 - s)
    return c1 * yCoords[1] + c2 * yCoords[2] + c3 * yCoords[3] + c4 * yCoords[4] + c5 * yCoords[5] + c6 * yCoords[6] + c7 * yCoords[7] + c8 * yCoords[8]
end

function dzt(r, s, t, zCoords::Array{Float64})
    # c_i --- coefficients of coordinate derivative
    c1 = 0.125 * (1 + r) * (1 + s)
    c2 = 0.125 * (1 - r) * (1 + s)
    c3 = 0.125 * (1 - r) * (1 - s)
    c4 = 0.125 * (1 + r) * (1 - s)
    c5 = -0.125 * (1 + r) * (1 + s) 
    c6 = -0.125 * (1 - r) * (1 + s)
    c7 = -0.125 * (1 - r) * (1 - s)
    c8 = -0.125 * (1 + r) * (1 - s)
    return c1 * zCoords[1] + c2 * zCoords[2] + c3 * zCoords[3] + c4 * zCoords[4] + c5 * zCoords[5] + c6 * zCoords[6] + c7 * zCoords[7] + c8 * zCoords[8]
end

# Coordinate system trasformation

function ElementTypes.jacGlobToLoc(r, s, t, xCoords::Array{Float64}, yCoords::Array{Float64}, zCoords::Array{Float64}, elemTypeInd::Iso8Pts3DType)
    return [dxr(r, s, t, xCoords) dyr(r, s, t, yCoords) dzr(r, s, t, zCoords)
            dxs(r, s, t, xCoords) dys(r, s, t, yCoords) dzs(r, s, t, zCoords)
            dxt(r, s, t, xCoords) dyt(r, s, t, yCoords) dzt(r, s, t, zCoords)]
end

function jacGlobToLocInv(r, s, t, xCoords::Array{Float64}, yCoords::Array{Float64}, zCoords::Array{Float64})
    elemTypeInd = Quad8Type("Iso8Pts3DType")
    inv(jacGlobToLoc(r, s, t, xCoords, yCoords, zCoords, elemTypeInd))
end

function ElementTypes.DetJs(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}, elemTypeInd::Iso8Pts3DType)
    return sqrt(dxs(r, s, t, xCoords)^2 + dys(r, s, t, yCoords)^2 + dzs(r, s, t, zCoords))
end

# Supporting matrices for calculating gradient matrix

function TU(r, s, t, xCoords::Array{Float64}, yCoords::Array{Float64}, zCoords::Array{Float64})
    rc1 = 0.125 * (1 + s) * (1 + t)
    rc2 = -0.125 * (1 + s) * (1 + t)
    rc3 = -0.125 * (1 - s) * (1 + t)
    rc4 = 0.125 * (1 - s) * (1 + t) 
    rc5 = 0.125 * (1 + s) * (1 - t) 
    rc6 = -0.125 * (1 + s) * (1 - t) 
    rc7 = -0.125 * (1 - s) * (1 - t)
    rc8 = 0.125 * (1 - s) * (1 - t)

    sc1 = 0.125 * (1 + r) * (1 + t)
    sc2 = 0.125 * (1 - r) * (1 + t)
    sc3 = -0.125 * (1 - r) * (1 + t) 
    sc4 = -0.125 * (1 + r) * (1 + t)
    sc5 = 0.125 * (1 + r) * (1 - t)
    sc6 = 0.125 * (1 - r) * (1 - t)
    sc7 = -0.125 * (1 - r) * (1 - t)
    sc8 = -0.125 * (1 + r) * (1 - t)

    tc1 = 0.125 * (1 + r) * (1 + s)
    tc2 = 0.125 * (1 - r) * (1 + s)
    tc3 = 0.125 * (1 - r) * (1 - s)
    tc4 = 0.125 * (1 + r) * (1 - s)
    tc5 = -0.125 * (1 + r) * (1 + s)
    tc6 = -0.125 * (1 - r) * (1 + s)
    tc7 = -0.125 * (1 - r) * (1 - s)
    tc8 = -0.125 * (1 + r) * (1 - s)

    uxyz = [rc1 0 0 rc2 0 0 rc3 0 0 rc4 0 0 rc5 0 0 rc6 0 0 rc7 0 0 rc8 0 0
            sc1 0 0 sc2 0 0 sc3 0 0 sc4 0 0 sc5 0 0 sc6 0 0 sc7 0 0 sc8 0 0
            tc1 0 0 tc2 0 0 tc3 0 0 tc4 0 0 tc5 0 0 tc6 0 0 tc7 0 0 tc8 0 0]

    return jacGlobToLocInv(r, s, t, xCoords, yCoords, zCoords) * uxyz
end

function TV(r, s, t, xCoords::Array{Float64}, yCoords::Array{Float64}, zCoords::Array{Float64})
    rc1 = 0.125 * (1 + s) * (1 + t)
    rc2 = -0.125 * (1 + s) * (1 + t)
    rc3 = -0.125 * (1 - s) * (1 + t)
    rc4 = 0.125 * (1 - s) * (1 + t) 
    rc5 = 0.125 * (1 + s) * (1 - t) 
    rc6 = -0.125 * (1 + s) * (1 - t) 
    rc7 = -0.125 * (1 - s) * (1 - t)
    rc8 = 0.125 * (1 - s) * (1 - t)

    sc1 = 0.125 * (1 + r) * (1 + t)
    sc2 = 0.125 * (1 - r) * (1 + t)
    sc3 = -0.125 * (1 - r) * (1 + t) 
    sc4 = -0.125 * (1 + r) * (1 + t)
    sc5 = 0.125 * (1 + r) * (1 - t)
    sc6 = 0.125 * (1 - r) * (1 - t)
    sc7 = -0.125 * (1 - r) * (1 - t)
    sc8 = -0.125 * (1 + r) * (1 - t)

    tc1 = 0.125 * (1 + r) * (1 + s)
    tc2 = 0.125 * (1 - r) * (1 + s)
    tc3 = 0.125 * (1 - r) * (1 - s)
    tc4 = 0.125 * (1 + r) * (1 - s)
    tc5 = -0.125 * (1 + r) * (1 + s)
    tc6 = -0.125 * (1 - r) * (1 + s)
    tc7 = -0.125 * (1 - r) * (1 - s)
    tc8 = -0.125 * (1 + r) * (1 - s)

    vxyz = [0 rc1 0 0 rc2 0 0 rc3 0 0 rc4 0 0 rc5 0 0 rc6 0 0 rc7 0 0 rc8 0
            0 sc1 0 0 sc2 0 0 sc3 0 0 sc4 0 0 sc5 0 0 sc6 0 0 sc7 0 0 sc8 0
            0 tc1 0 0 tc2 0 0 tc3 0 0 tc4 0 0 tc5 0 0 tc6 0 0 tc7 0 0 tc8 0]

    return jacGlobToLocInv(r, s, t, xCoords, yCoords, zCoords) * vxyz
end

function TW(r, s, t, xCoords::Array{Float64}, yCoords::Array{Float64}, zCoords::Array{Float64})
    rc1 = 0.125 * (1 + s) * (1 + t)
    rc2 = -0.125 * (1 + s) * (1 + t)
    rc3 = -0.125 * (1 - s) * (1 + t)
    rc4 = 0.125 * (1 - s) * (1 + t) 
    rc5 = 0.125 * (1 + s) * (1 - t) 
    rc6 = -0.125 * (1 + s) * (1 - t) 
    rc7 = -0.125 * (1 - s) * (1 - t)
    rc8 = 0.125 * (1 - s) * (1 - t)

    sc1 = 0.125 * (1 + r) * (1 + t)
    sc2 = 0.125 * (1 - r) * (1 + t)
    sc3 = -0.125 * (1 - r) * (1 + t) 
    sc4 = -0.125 * (1 + r) * (1 + t)
    sc5 = 0.125 * (1 + r) * (1 - t)
    sc6 = 0.125 * (1 - r) * (1 - t)
    sc7 = -0.125 * (1 - r) * (1 - t)
    sc8 = -0.125 * (1 + r) * (1 - t)

    tc1 = 0.125 * (1 + r) * (1 + s)
    tc2 = 0.125 * (1 - r) * (1 + s)
    tc3 = 0.125 * (1 - r) * (1 - s)
    tc4 = 0.125 * (1 + r) * (1 - s)
    tc5 = -0.125 * (1 + r) * (1 + s)
    tc6 = -0.125 * (1 - r) * (1 + s)
    tc7 = -0.125 * (1 - r) * (1 - s)
    tc8 = -0.125 * (1 + r) * (1 - s)

    wxyz = [0 0 rc1 0 0 rc2 0 0 rc3 0 0 rc4 0 0 rc5 0 0 rc6 0 0 rc7 0 0 rc8
            0 0 sc1 0 0 sc2 0 0 sc3 0 0 sc4 0 0 sc5 0 0 sc6 0 0 sc7 0 0 sc8
            0 0 tc1 0 0 tc2 0 0 tc3 0 0 tc4 0 0 tc5 0 0 tc6 0 0 tc7 0 0 tc8]

    return jacGlobToLocInv(r, s, t, xCoords, yCoords, zCoords) * wxyz
end

function ElementTypes.gradMatr(r, s, t, xCoords::Array{Float64}, yCoords::Array{Float64}, zCoords::Array{Float64}, elemTypeInd::Iso8Pts3DType)
    if size(xCoords)[1] != size(yCoords)[1] != size(zCoords)[1]
        println("Incorrect nodes while calculating gradient matrix")
        return nothing
    end
    nodesPerElement = size(xCoords)[1]
    resultMatrix = Array{Float64, 2}(undef, 6, nodesPerElement * 3)
    Uxyz = TU(r, s, t, xCoords, yCoords, zCoords)
    Vxyz = TV(r, s, t, xCoords, yCoords, zCoords)
    Wxyz = TW(r, s, t, xCoords, yCoords, zCoords)
    if (size(Uxyz)[2] != size(Vxyz)[2] != size(Wxyz) != size(resultMatrix)[2])
        println("Error while calculating gradient matrix")
        return nothing
    end
    cols = size(Uxyz)[2]
    for i in 1:cols
        resultMatrix[1, i] = Uxyz[1, i]
        resultMatrix[2, i] = Vxyz[2, i]
        resultMatrix[3, i] = Wxyz[3, i]
        resultMatrix[4, i] = Uxyz[2, i] + Vxyz[1, i]
        resultMatrix[5, i] = Vxyz[3, i] + Wxyz[2, i]
        resultMatrix[6, i] = Wxyz[1, i] + Uxyz[3, i]
    end
    return resultMatrix
end

function ElementTypes.displInterpMatr(r, s, elemTypeInd::Iso8Pts3DType)
    result = [h1(r, s, t) 0 0 h2(r, s, t) 0 0 h3(r, s, t) 0 0 h4(r, s, t) 0 0 h5(r, s, t) 0 0 h6(r, s, t) 0 0 h7(r, s, t) 0 0 h8(r, s, t) 0 0
              0 h1(r, s, t) 0 0 h2(r, s, t) 0 0 h3(r, s, t) 0 0 h4(r, s, t) 0 0 h5(r, s, t) 0 0 h6(r, s, t) 0 0 h7(r, s, t) 0 0 h8(r, s, t) 0
              0 0 h1(r, s, t) 0 0 r2(r, s, t) 0 0 r3(r, s, t) 0 0 r4(r, s, t) 0 0 r5(r, s, t) 0 0 r6(r, s, t) 0 0 r7(r, s, t) 0 0 h8(r, s, t)]
    return result
end  # displInterpMatr

function ElementTypes.nodesFromDirection(direction::Int, elemTypeInd::Iso8Pts3DType)
    if direction == 1  # Top
        return [1, 2, 3, 4]
    elseif direction == 2  # Left
        return [3, 4, 8, 7]
    elseif direction == 3  # Bottom
        return [5, 6, 7, 8]
    elseif direction == 4  # Right
        return [1, 2, 6, 5]
    elseif direction == 5  # To us
        return [1, 4, 8, 5]
    elseif direction == 6  # From us
        return [2, 3, 7, 6]
    else
        println("Given load direction is not supported")
        return nothing
    end
end  # nodesFromDirection

function ElementTypes.getRSFromNode(nodeIndex::Int, elemTypeInd::Iso8Pts3DType)
    if nodeIndex == 1
        return (1, 1, 1)
    elseif nodeIndex == 2
        return (-1, 1, 1)
    elseif nodeIndex == 3
        return (-1, -1, 1)
    elseif nodeIndex == 4
        return (1, -1, 1)
    elseif nodeIndex == 5
        return (1, 1, -1)
    elseif nodeIndex == 6
        return (-1, 1, -1)
    elseif nodeIndex == 7
        return (-1, -1, -1)
    elseif nodeIndex == 8
        return (1, -1, -1)
    else
        println("Invalid node index while getting (r, s, t) coordinates from node")
        return nothing
    end
end  # getRSFromNode

end  # Iso8Pts3D