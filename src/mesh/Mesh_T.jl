export Mesh2D_T, generateTestMesh2D, printNodesMesh2D, printElementsMesh2D

struct Mesh2D_T
    nodes::Vector{Tuple{Vararg{Float64}}}  # Matrix of nodes (contains node's coordinates)
    elements::Vector{Tuple{Vararg{Int}}}  # Matrix of elements (contains elements's nodes)
end  # Mesh2D_T

function generateTestMesh2D()  # Simple mesh built on a square 1x1
    # Parameters
    nOfNodes = 9
    nOfElements = 4  # For square mesh it should be square of smth
    dimension = 2
    nodesPerElement = 4
    # Creation
    resultMesh = Mesh2D_T(Vector{Tuple{Vararg{Float64}}}(undef, nOfNodes), Vector{Tuple{Vararg{Int}}}(undef, nOfElements))
    nodeIndex = 1
    for i in 0:50:100
        for j in 0:50:100
            resultMesh.nodes[nodeIndex] = (j, i)
            nodeIndex += 1
        end
    end
    for i in eachindex(resultMesh.elements)
        blNum = i + div(i - 1, sqrt(nOfElements))  # Bottom left node of element
        element = Tuple{blNum, blNum + 1, blNum + 4, blNum + 3}
        resultMesh.elements[i] = (blNum, blNum + 1, blNum + 4, blNum + 3)
    end
    return resultMesh
end  # generateRandomMesh2D

function printNodesMesh2D(mesh::Mesh2D_T)
    for node in mesh.nodes
        print(node, "\n")
    end
    print("\n")
end  # printMesh2D

function printElementsMesh2D(mesh::Mesh2D_T)
    for element in mesh.elements
        print(element, "\n")
    end
    print("\n")
end  # printElementsMesh2D
