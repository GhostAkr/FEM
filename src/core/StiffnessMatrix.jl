# Calculating stiffness matrix

export stiffnessMatrix

using Quad8Pts
using Quad4Pts
using MultipleIntegral
using LinearAlgebra

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
function F(r, s, xCoords::Array{Float64}, yCoords::Array{Float64}, elasticityMatrix::AbstractArray, elemTypeInd::FiniteElement)  # F = B^T * C * B * det(J)
    B = gradMatr(r, s, xCoords, yCoords, elemTypeInd)  # Gradient matrix
    BTransp = transpose(B)
    J = jacGlobToLoc(r, s, xCoords, yCoords, elemTypeInd)  # Jacobi's matrix
    
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
function stiffnessMatrix(elasticityMatrix::AbstractArray, parameters::processPars, elementNum::Int, intOrder::Int, elemTypeInd::FiniteElement)
    nodesPerElement = length(parameters.mesh.elements[elementNum])
    xCoords = [parameters.mesh.nodes[parameters.mesh.elements[elementNum][i]][1] for i in 1:nodesPerElement]
    yCoords = [parameters.mesh.nodes[parameters.mesh.elements[elementNum][i]][2] for i in 1:nodesPerElement]
    FIntegrate(r, s) = F(r, s, xCoords, yCoords, elasticityMatrix, elemTypeInd)  # F representation for integrating (depends only on r and s)
    K = MultipleIntegral.gaussmethod_matrix_2d(FIntegrate, intOrder)
    return K
end  # stiffnessMatrix

function F3D(r, s, t, xCoords::Array{Float64}, yCoords::Array{Float64}, zCoords::Array{Float64}, elasticityMatrix::AbstractArray, elemTypeInd::FiniteElement)
    B = gradMatr(r, s, t, xCoords, yCoords, zCoords, elemTypeInd)
    BTransp = transpose(B)
    J = jacGlobToLoc(r, s, t, xCoords, yCoords, zCoords, elemTypeInd)

    return BTransp * elasticityMatrix * B * det(J)
end

function stiffnessMatrix3D(elasticityMatrix::AbstractArray, parameters::processPars, elementNum::Int, intOrder::Int, elemTypeInd::FiniteElement)
    nodesPerElement = length(parameters.mesh.elements[elementNum])
    xCoords = [parameters.mesh.nodes[parameters.mesh.elements[elementNum][i]][1] for i in 1:nodesPerElement]
    yCoords = [parameters.mesh.nodes[parameters.mesh.elements[elementNum][i]][2] for i in 1:nodesPerElement]
    zCoords = [parameters.mesh.nodes[parameters.mesh.elements[elementNum][i]][3] for i in 1:nodesPerElement]
    FIntegrate(r, s, t) = F3D(r, s, t, xCoords, yCoords, zCoords, elasticityMatrix, elemTypeInd)
    K = MultipleIntegral.gaussmethod_matrix_3d(FIntegrate, intOrder)
    return K
end
