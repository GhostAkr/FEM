# Mesh type
module Mod_Mesh_T

    # Import block
    import Core

    struct Mesh2D_T
        nodes::AbstractVector{Tuple{Vararg{Real}}}  # Matrix of nodes (contains node's coordinates)
        elements::AbstractVector{Tuple{Vararg{Int}}}  # Matrix of elements (contains elements's nodes)
    end  # Mesh2D_T

    function generateTestMesh2D()  # Simple mesh built on a square 1x1
        # Parameters
        nOfNodes = 25
        nOfElements = 16
        dimension = 2
        nodesPerElement = 4
        # Creation
        resultMesh = Mesh2D_T(Vector{Tuple{Vararg{Real}}}(undef, nOfNodes), Vector{Tuple{Vararg{Int}}}(undef, nOfElements))
        nodeIndex = 1
        for i in 0:0.25:1
            for j in 0:0.25:1
                print("i = ", i, "; j = ", j, "\n")
                resultMesh.nodes[nodeIndex] = (j, i)
                nodeIndex += 1
            end
        end
        for i in eachindex(resultMesh.elements)
            blNum = i + div(i - 1, nodesPerElement)  # Bottom left node of element
            element = Tuple{blNum, blNum + 1, blNum + 6, blNum + 5}
            resultMesh.elements[i] = (blNum, blNum + 1, blNum + 6, blNum + 5)
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

    # Export block
    export Mesh2D_T, generateTestMesh2D, printElementsMesh2D, printNodesMesh2D

end  # Mod_Mesh_T
