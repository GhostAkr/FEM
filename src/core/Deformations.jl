using Quad4Pts
using Quad8Pts

function calculateDeformations(displacements::Array, pars::processPars)
    calculatedNodes = []
    nOfNodes = size(pars.mesh.nodes)[1]
    deformations = Array{Array{Float64}}(undef, nOfNodes)
    for element in pars.mesh.elements
        elementDisplacements = []
        for node in element
            push!(elementDisplacements, displacements[2 * node - 1])
            push!(elementDisplacements, displacements[2 * node])
        end
        for nodeIndex in eachindex(element)
            if findfirst(isequal(element[nodeIndex]), calculatedNodes) !== nothing
                continue
            end
            rCoord = Quad8Pts.getRSFromNode(nodeIndex)[1]
            sCoord = Quad8Pts.getRSFromNode(nodeIndex)[2]
            nodesPerElement = length(element)
            xCoords = [pars.mesh.nodes[element[i]][1] for i in 1:nodesPerElement]
            yCoords = [pars.mesh.nodes[element[i]][2] for i in 1:nodesPerElement]
            nodeDeformations = Quad8Pts.gradMatr(rCoord, sCoord, xCoords, yCoords) * elementDisplacements
            deformations[element[nodeIndex]] = nodeDeformations
            push!(calculatedNodes, element[nodeIndex])
        end
    end
    return deformations
end  # calculateDeformations