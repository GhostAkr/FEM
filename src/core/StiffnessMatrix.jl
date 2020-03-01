export stiffnessMatrix

using Quad4Pts

function stiffnessMatrix(elasticityMatrix::AbstractArray, parameters::processPars, elementNum::Int)
    xCoords = [parameters.mesh.nodes[parameters.mesh.elements[elementNum][i]][1] for i in 1:4]
    yCoords = [parameters.mesh.nodes[parameters.mesh.elements[elementNum][i]][2] for i in 1:4]
    x = 0
    for i in 1:length(Quad4Pts.interFunc)
        x[i] += xCoords[i] * Quad4Pts.interFunc[i]
    end
    x(r, s) = sum(interFunc[i] * )
end
