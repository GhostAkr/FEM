using Quad4Pts
using Quad8Pts
using multipleIntegral

function calculateDeformations(displacements::Array, pars::processPars, intOrder::Int, elemTypeInd::FiniteElement)
    nOfNodes = size(pars.mesh.nodes)[1]
    nOfDeformationsTypes = 3
    deformations = zeros(Float64, nOfNodes, nOfDeformationsTypes)
    averagingNums = zeros(Int, nOfNodes)
    mid = 0
    num = 0
    for element in pars.mesh.elements
        elementDisplacements = []
        for node in element
            push!(elementDisplacements, displacements[2 * node - 1])
            push!(elementDisplacements, displacements[2 * node])
        end
        # gaussPoints = multipleIntegral.getGaussPoints2D(intOrder)
        for nodeIndex in eachindex(element)
            # rCoord = gaussPoints[nodeIndex][1]
            # sCoord = gaussPoints[nodeIndex][2]
            rCoord = getRSFromNode(nodeIndex, elemTypeInd)[1]
            sCoord = getRSFromNode(nodeIndex, elemTypeInd)[2]
            nodesPerElement = length(element)
            xCoords = [pars.mesh.nodes[element[i]][1] for i in 1:nodesPerElement]
            yCoords = [pars.mesh.nodes[element[i]][2] for i in 1:nodesPerElement]
            gaussNodeDeformations = gradMatr(rCoord, sCoord, xCoords, yCoords, elemTypeInd::FiniteElement) * elementDisplacements
            for typeIndex in 1:nOfDeformationsTypes
                deformations[element[nodeIndex], typeIndex] += gaussNodeDeformations[typeIndex]
            end
            averagingNums[element[nodeIndex]] += 1
        end
    end
    for nodeIndex in 1:size(deformations)[1]
        for typeIndex in 1:size(deformations)[2]
            deformations[nodeIndex, typeIndex] /= averagingNums[nodeIndex]
        end
    end
    return deformations
end  # calculateDeformations