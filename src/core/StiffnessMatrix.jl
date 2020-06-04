# Calculating stiffness matrix

export stiffnessMatrix

using Quad4Pts
using multipleIntegral
using LinearAlgebra

"""
    testIntegralFunc(x, y)

Just random function depending on 2 variables: ``f(x, y) = x^2 + y``.
"""
testIntegralFunc(x, y) = x ^ 2 + y  # Should give 4/3 after integration

"""
    F(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}, elasticityMatrix::AbstractArray)

Function that should be integrated to calculate local stiffness matrix for element: ``F = B^T \\cdot C \\cdot B \\cdot det(J)``, 
where B is strain-displacement tramsformation matrix (in natural coordinate system), 
C is elasticity matrix and J is Jacobi's matrix.

# Arguments
- `r`: r-coordinate;
- `s`: s-coordinate;
- `xCoords::Array{Float64}`: x coordinates of each node in current element;
- `yCoords::Array{Float64}`: y coordinates of each node in current element;
- `elasticityMatrix::AbstractArray`: elasticity matrix of current element.
"""
function F(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}, elasticityMatrix::AbstractArray)  # F = B^T * C * B * det(J)
    B = Quad4Pts.gradMatr(r, s, xCoords, yCoords)  # Gradient matrix
    BTransp = transpose(B)
    J = Quad4Pts.jacGlobToLoc(r, s, xCoords, yCoords)  # Jacobi's matrix
    return BTransp * elasticityMatrix * B * det(J)
end  # F

"""
    stiffnessMatrix(elasticityMatrix::AbstractArray, parameters::processPars, elementNum::Int)

Calculate stiffness matrix of given element: ``K = \\int \\limits_S B^T \\cdot C \\cdot B \\cdot det(J) dS = 
\\int \\limits_S F dS``.

# Arguments
- `elasticityMatrix::AbstractArray`: elasticity matrix of current element;
- `parameters::processPars`: parameters of current model;
- `elementNum::Int`: number of given element.
"""
function stiffnessMatrix(elasticityMatrix::AbstractArray, parameters::processPars, elementNum::Int)
    xCoords = [parameters.mesh.nodes[parameters.mesh.elements[elementNum][i]][1] for i in 1:4]
    yCoords = [parameters.mesh.nodes[parameters.mesh.elements[elementNum][i]][2] for i in 1:4]
    FIntegrate(r, s) = F(r, s, xCoords, yCoords, elasticityMatrix)  # F representation for integrating (depends only on r and s)
    IntegrationOrder = 4
    K = multipleIntegral.gaussMethodMatrix(FIntegrate, IntegrationOrder)
    return K
end  # stiffnessMatrix
