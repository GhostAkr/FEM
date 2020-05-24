"""
    Quad4Pts

Module describing model with bilinear quadrilateral finite elements.
"""
module Quad4Pts

h1(r, s) = 0.25 * (1 - r) * (1 - s)
h2(r, s) = 0.25 * (1 + r) * (1 - s)
h3(r, s) = 0.25 * (1 + r) * (1 + s)
h4(r, s) = 0.25 * (1 - r) * (1 + s)

x(r, s, xCoords::Array{Float64}) = h1(r, s) * xCoords[1] + h2(r, s) * xCoords[2] + h3(r, s) * xCoords[3] + h4(r, s) * xCoords[4]
y(r, s, yCoords::Array{Float64}) = h1(r, s) * yCoords[1] + h2(r, s) * yCoords[2] + h3(r, s) * yCoords[3] + h4(r, s) * yCoords[4]

u(r, s, uCoords::Array{Float64}) = h1(r, s) * uCoords[1] + h2(r, s) * uCoords[2] + h3(r, s) * uCoords[3] + h4(r, s) * uCoords[4]
v(r, s, vCoords::Array{Float64}) = h1(r, s) * vCoords[1] + h2(r, s) * vCoords[2] + h3(r, s) * vCoords[3] + h4(r, s) * vCoords[4]

# TODO: provide somehow dynamic derivatives
dxr(r, s, xCoords::Array{Float64}) = -0.25 * (1 - s) * xCoords[1] + 0.25 * (1 - s) * xCoords[2] + 0.25 * (1 + s) * xCoords[3] - 0.25 * (1 + s) * xCoords[4]
dxs(r, s, xCoords::Array{Float64}) = -0.25 * (1 - r) * xCoords[1] - 0.25 * (1 + r) * xCoords[2] + 0.25 * (1 + r) * xCoords[3] - 0.25 * (1 - r) * xCoords[4]
dyr(r, s, yCoords::Array{Float64}) = -0.25 * (1 - s) * yCoords[1] + 0.25 * (1 - s) * yCoords[2] + 0.25 * (1 + s) * yCoords[3] - 0.25 * (1 + s) * yCoords[4]
dys(r, s, yCoords::Array{Float64}) = -0.25 * (1 - r) * yCoords[1] - 0.25 * (1 + r) * yCoords[2] + 0.25 * (1 + r) * yCoords[3] + 0.25 * (1 - r) * yCoords[4]

dur(r, s, uCoords::Array{Float64}) = -0.25 * (1 - s) * uCoords[1] + 0.25 * (1 - s) * uCoords[2] + 0.25 * (1 + s) * uCoords[3] - 0.25 * (1 + s) * uCoords[4]
dus(r, s, uCoords::Array{Float64}) = -0.25 * (1 - r) * uCoords[1] - 0.25 * (1 + r) * uCoords[2] + 0.25 * (1 + r) * uCoords[3] + 0.25 * (1 - r) * uCoords[4]
dvr(r, s, vCoords::Array{Float64}) = -0.25 * (1 - s) * vCoords[1] + 0.25 * (1 - s) * vCoords[2] + 0.25 * (1 + s) * vCoords[3] - 0.25 * (1 + s) * vCoords[4]
dvs(r, s, vCoords::Array{Float64}) = -0.25 * (1 - r) * vCoords[1] - 0.25 * (1 + r) * vCoords[2] + 0.25 * (1 + r) * vCoords[3] + 0.25 * (1 - r) * vCoords[4]

# TODO: provide somehow dynamic derivatives
dh1r(r, s) = (-1 + s) / 4
dh1s(r, s) = (-1 + r) / 4
dh2r(r, s) = (1 - s) / 4
dh2s(r, s) = (-1 - r) / 4
dh3r(r, s) = (1 + s) / 4
dh3s(r, s) = (1 + r) / 4
dh4r(r, s) = (-1 - s) / 4
dh4s(r, s) = (1 - r) / 4

"""
    jacGlobToLoc(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})

Jacobi matrix for conversion between different coordinate systems.

# Arguments
- `r`: r-coordinate;
- `s`: s-coordinate;
- `xCoords::Array{Float64}`: x coordinates of each node in current element;
- `yCoords::Array{Float64}`: y coordinates of each node in current element.
"""
jacGlobToLoc(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = [dxr(r, s, xCoords) dxs(r, s, xCoords); dyr(r, s, yCoords) dys(r, s, yCoords)]  # Jacobi's matrix for conversion from local coordinates to global

"""
    jacGlobToLocInv(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})

Inversed Jacobi matrix for conversion between different coordinate systems.

# Arguments
- `r`: r-coordinate;
- `s`: s-coordinate;
- `xCoords::Array{Float64}`: x coordinates of each node in current element;
- `yCoords::Array{Float64}`: y coordinates of each node in current element.
"""
jacGlobToLocInv(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = inv(jacGlobToLoc(r, s, xCoords, yCoords))

# Supporting matrices for calculating gradient matrix
function TU(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})
    uxy = [-(1 - s)    0   1 - s    0   1 + s    0   -(1 + s)       0
           -(1 - r)    0   -(1 + r)       0   1 + r    0   1 - r    0]
    return 0.25 * jacGlobToLocInv(r, s, xCoords, yCoords) * uxy
end  # TU

function TV(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})
    vxy = [0    -(1 - s)   0   1 - s    0   1 + s    0   -(1 + s)
           0    -(1 - r)   0   -(1 + r)       0   1 + r    0   1 - r]
    return 0.25 * jacGlobToLocInv(r, s, xCoords, yCoords) * vxy
end  # TU

"""
    DetJs(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})

Determinant of appropriate Jacobi matrix.

# Arguments
- `r`: r-coordinate;
- `s`: s-coordinate;
- `xCoords::Array{Float64}`: x coordinates of each node in current element;
- `yCoords::Array{Float64}`: y coordinates of each node in current element.
"""
DetJs(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = sqrt(dxs(r, s, xCoords)^2 + dys(r, s, yCoords)^2)

dh1x(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = jacGlobToLocInv(r, s, xCoords, yCoords)[1, 1] * dh1r(r, s) + jacGlobToLocInv(r, s, xCoords, yCoords)[1, 2] * dh1s(r, s)
dh1y(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = jacGlobToLocInv(r, s, xCoords, yCoords)[2, 1] * dh1r(r, s) + jacGlobToLocInv(r, s, xCoords, yCoords)[2, 2] * dh1s(r, s)
dh2x(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = jacGlobToLocInv(r, s, xCoords, yCoords)[1, 1] * dh2r(r, s) + jacGlobToLocInv(r, s, xCoords, yCoords)[1, 2] * dh2s(r, s)
dh2y(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = jacGlobToLocInv(r, s, xCoords, yCoords)[2, 1] * dh2r(r, s) + jacGlobToLocInv(r, s, xCoords, yCoords)[2, 2] * dh2s(r, s)
dh3x(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = jacGlobToLocInv(r, s, xCoords, yCoords)[1, 1] * dh3r(r, s) + jacGlobToLocInv(r, s, xCoords, yCoords)[1, 2] * dh3s(r, s)
dh3y(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = jacGlobToLocInv(r, s, xCoords, yCoords)[2, 1] * dh3r(r, s) + jacGlobToLocInv(r, s, xCoords, yCoords)[2, 2] * dh3s(r, s)
dh4x(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = jacGlobToLocInv(r, s, xCoords, yCoords)[1, 1] * dh4r(r, s) + jacGlobToLocInv(r, s, xCoords, yCoords)[1, 2] * dh4s(r, s)
dh4y(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}) = jacGlobToLocInv(r, s, xCoords, yCoords)[2, 1] * dh4r(r, s) + jacGlobToLocInv(r, s, xCoords, yCoords)[2, 2] * dh4s(r, s)

# TODO: For now this is not effective
function gradMatr(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})
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

# TODO: calculate this through appropriate derivatives
# Precalculated version
"""
    gradMatr(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})

Gradient matrix ``B`` for element.

# Arguments
- `r`: r-coordinate;
- `s`: s-coordinate;
- `xCoords::Array{Float64}`: x coordinates of each node in current element;
- `yCoords::Array{Float64}`: y coordinates of each node in current element.
"""
function gradMatrTmp(r, s, xCoords::Array{Float64}, yCoords::Array{Float64})
    resultMatrix = Array{Float64, 2}(undef, 3, 8)
    resultMatrix[1, 1] = 1 + s
    resultMatrix[1, 2] = 0
    resultMatrix[1, 3] = -(1 + s)
    resultMatrix[1, 4] = 0
    resultMatrix[1, 5] = -(1 - s)
    resultMatrix[1, 6] = 0
    resultMatrix[1, 7] = 1 - s
    resultMatrix[1, 8] = 0

    resultMatrix[2, 1] = 0
    resultMatrix[2, 2] = 1 + r
    resultMatrix[2, 3] = 0
    resultMatrix[2, 4] = 1 - r
    resultMatrix[2, 5] = 0
    resultMatrix[2, 6] = -(1 - r)
    resultMatrix[2, 7] = 0
    resultMatrix[2, 8] = -(1 + r)

    resultMatrix[3, 1] = 1 + r
    resultMatrix[3, 2] = 1 + s
    resultMatrix[3, 3] = 1 - r
    resultMatrix[3, 4] = -(1 + s)
    resultMatrix[3, 5] = -(1 - r)
    resultMatrix[3, 6] = -(1 - s)
    resultMatrix[3, 7] = -(1 + r)
    resultMatrix[3, 8] = 1 - s

    resultMatrix *= 0.25

    return resultMatrix
end  # gradMatr

# Precalculated displacements interpolation H matrix
"""
    displInterpMatr(r, s)

Displacements interpolation ``H`` matrix.

# Arguments
- `r`: r-coordinate;
- `s`: s-coordinate.
"""
function displInterpMatr(r, s)
    result = zeros(Real, 2, 8)
    result[1, 3] = 0.5 * (1 - s)
    result[1, 5] = 0.5 * (1 + s)
    result[2, 4] = 0.5 * (1 - s)
    result[2, 6] = 0.5 * (1 + s)
    return result
end  # displInterpMatr

interFunc = [h1, h2, h3, h4]  # Array of interpolation functions

end  # Quad4Pts