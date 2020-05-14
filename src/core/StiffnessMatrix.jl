export stiffnessMatrix

using Quad4Pts
using multipleIntegral
using LinearAlgebra

testIntegralFunc(x, y) = x ^ 2 + y  # Should give 4/3 after integration

function F(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}, elasticityMatrix::AbstractArray, elNum::Int)  # F = B^T * C * B * det(J)
    B = Quad4Pts.gradMatr(r, s, xCoords, yCoords)  # Gradient matrix
    # if elNum == 1 && xCoords == [0, 0.25, 0.25, 0] && yCoords == [0, 0, 0.25, 0.25]
    #     println(Quad4Pts.gradMatr(5, 10, xCoords, yCoords), "\n")
    # end
    BTransp = transpose(B)
    J = Quad4Pts.jacGlobToLoc(r, s, xCoords, yCoords)  # Jacobi's matrix
    #return B * elasticityMatrix * BTransp * det(J)
    return BTransp * elasticityMatrix * B * det(J)
end

function stiffnessMatrix(elasticityMatrix::AbstractArray, parameters::processPars, elementNum::Int)
    xCoords = [parameters.mesh.nodes[parameters.mesh.elements[elementNum][i]][1] for i in 1:4]
    yCoords = [parameters.mesh.nodes[parameters.mesh.elements[elementNum][i]][2] for i in 1:4]
    FIntegrate(r, s) = F(r, s, xCoords, yCoords, elasticityMatrix, elementNum)  # F representation for integrating (depends only on r and s)
    IntegrationOrder = 2
    K = multipleIntegral.gaussMethodMatrix(FIntegrate, IntegrationOrder)
    return K
end
