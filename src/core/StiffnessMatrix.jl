export stiffnessMatrix

using Quad4Pts
using multipleIntegral
using LinearAlgebra

testIntegralFunc(x, y) = x ^ 2 + y  # Should give 4/3 after integration

function F(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}, elasticityMatrix::AbstractArray)  # F = B^T * C * B * det(J)
    B = Quad4Pts.gradMatr(r, s, xCoords, yCoords)  # Gradient matrix
    BTransp = transpose(B)
    J = Quad4Pts.jacGlobToLoc(r, s, xCoords, yCoords)  # Jacobi's matrix
    #return B * elasticityMatrix * BTransp * det(J)
    return BTransp * elasticityMatrix * B * det(J)
end



function stiffnessMatrix(elasticityMatrix::AbstractArray, parameters::processPars, elementNum::Int)
    xCoords = [parameters.mesh.nodes[parameters.mesh.elements[elementNum][i]][1] for i in 1:4]
    yCoords = [parameters.mesh.nodes[parameters.mesh.elements[elementNum][i]][2] for i in 1:4]
    thickness = 1
    rLim = Dict{limits, Real}(lower => -1, upper => 1)
    sLim = Dict{limits, Real}(lower => -1, upper => 1)
    M = 500  # Number of sections by r
    N = 500  # Number of sections by s
    FIntegrate(r, s) = F(r, s, xCoords, yCoords, elasticityMatrix)  # F representation for integrating (depends only on r and s)
    K = multipleIntegral.cellMethodMatrix(FIntegrate, rLim, sLim, M, N)
    return K
end

