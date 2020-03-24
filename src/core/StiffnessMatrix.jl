export stiffnessMatrix

using Quad4Pts
using HCubature
using multipleIntegral

function F(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}, elasticityMatrix::AbstractArray)  # F = B^T * C * B * det(J)
    B = Quad4Pts.gradMatr(r, s, xCoords, yCoords)  # Gradient matrix
    BTransp = transpose(B)
    J = Quad4Pts.jacGlobToLoc(r, s, xCoords, yCoords)  # Jacobi's matrix
    return B * elasticityMatrix * BTransp * det(J)
end

function stiffnessMatrix(elasticityMatrix::AbstractArray, parameters::processPars, elementNum::Int)
    xCoords = [parameters.mesh.nodes[parameters.mesh.elements[elementNum][i]][1] for i in 1:4]
    yCoords = [parameters.mesh.nodes[parameters.mesh.elements[elementNum][i]][2] for i in 1:4]
    thickness = 1
    rLim = Dict(lower => -1, upper => 1)
    sLim = Dict(lower => -1, upper => 1)
    M = 6  # Number of sections by r
    N = 6  # Number of sections by s
    multipleIntegral.cellMethod(F, rLim, sLim, M, N)
end

