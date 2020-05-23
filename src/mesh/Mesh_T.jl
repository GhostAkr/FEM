export Mesh2D_T, generateTestMesh2D, printNodesMesh2D, printElementsMesh2D

"""
    Mesh2D_T

2D mesh main structure.
Basically contains 2 fields:
1. **nodes** - matrix of nodes (contains node's coordinates);
2. **elements** - matrix of elements (contains elements's nodes).
"""
struct Mesh2D_T
    "Matrix of nodes (contains node's coordinates)"
    nodes::Vector{Tuple{Vararg{Float64}}}
    "Matrix of elements (contains elements's nodes)"
    elements::Vector{Tuple{Vararg{Int}}}
end  # Mesh2D_T

"""
    generateTestMesh2D()

Build a simple 2D mesh on a square 100x100 with ``n^2`` elements on it.
"""
function generateTestMesh2D(n::Int)
    # Parameters
    nOfNodes = (n + 1)^2
    nOfElements = n^2  # For square mesh it should be square of smth
    dimension = 2
    nodesPerElement = 4
    # Creation
    resultMesh = Mesh2D_T(Vector{Tuple{Vararg{Float64}}}(undef, nOfNodes), Vector{Tuple{Vararg{Int}}}(undef, nOfElements))
    nodeIndex = 1
    step = 100 / n
    for i in 0:step:100
        for j in 0:step:100
            resultMesh.nodes[nodeIndex] = (j, i)
            nodeIndex += 1
        end
    end
    for i in eachindex(resultMesh.elements)
        blNum = i + div(i - 1, sqrt(nOfElements))  # Bottom left node of element
        element = Tuple{blNum, blNum + 1, blNum + (n + 1) + 1, blNum + (n + 1)}
        resultMesh.elements[i] = (blNum, blNum + 1, blNum + (n + 1) + 1, blNum + (n + 1))
    end
    return resultMesh
end  # generateRandomMesh2D


"""
    printNodesMesh2D(mesh::Mesh2D_T)

Print nodes of given mesh.
"""
function printNodesMesh2D(mesh::Mesh2D_T)
    for node in mesh.nodes
        print(node, "\n")
    end
    print("\n")
end  # printMesh2D


"""
    printElementsMesh2D(mesh::Mesh2D_T)

Print elements of given mesh.
"""
function printElementsMesh2D(mesh::Mesh2D_T)
    for element in mesh.elements
        print(element, "\n")
    end
    print("\n")
end  # printElementsMesh2D
