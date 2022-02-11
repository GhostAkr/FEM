using Quad4Pts
using Quad8Pts
using Iso8Pts3D
using MultipleIntegral

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
        # gaussPoints = MultipleIntegral.getgausspoints_2d(intOrder)
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

"""
    calculate_deformations_3d(displacements::Array, pars::processPars, 
        elemtype::FiniteElement)

Calculate 3D deformations by given displacements vector.

# Arguments
- `displacements::Array`: displacements vector;
- `pars::processPars`: process parameters;
- `elemtype::FiniteElement`: type of finite element.
"""
function calculate_deformations_3d(displacements::Array, pars::processPars, 
    elemtype::FiniteElement
)
    n_of_nodes = size(pars.mesh.nodes)[1]
    n_of_deformations_types = 6
    deformations = zeros(Float64, n_of_nodes, n_of_deformations_types)
    averaging_nums = zeros(Int, n_of_nodes)
    for element in pars.mesh.elements
        element_displacements = []
        for node in element
            push!(element_displacements, displacements[3 * node - 2])
            push!(element_displacements, displacements[3 * node - 1])
            push!(element_displacements, displacements[3 * node])
        end
        # gaussPoints = MultipleIntegral.getgausspoints_2d(intOrder)
        for nodeidx in eachindex(element)
            # rCoord = gaussPoints[nodeidx][1]
            # sCoord = gaussPoints[nodeidx][2]
            loc_coords = getRSFromNode(nodeidx, elemtype)
            nodes_per_elem = length(element)
            x_coords = [pars.mesh.nodes[element[i]][1] for i in 1:nodes_per_elem]
            y_coords = [pars.mesh.nodes[element[i]][2] for i in 1:nodes_per_elem]
            z_coords = [pars.mesh.nodes[element[i]][3] for i in 1:nodes_per_elem]
            gauss_node_deformations = gradMatr(loc_coords[1], loc_coords[2], loc_coords[3], 
                x_coords, y_coords, z_coords, elemtype) * element_displacements
            for typeidx in 1:n_of_deformations_types
                deformations[element[nodeidx], typeidx] += gauss_node_deformations[typeidx]
            end
            averaging_nums[element[nodeidx]] += 1
        end
    end
    for nodeidx in 1:size(deformations)[1]
        for typeidx in 1:size(deformations)[2]
            deformations[nodeidx, typeidx] /= averaging_nums[nodeidx]
        end
    end
    return deformations
end  # calculateDeformations
