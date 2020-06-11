include("MeshVars.jl")

export Mesh2D_T, generateTestMesh2D, printNodesMesh2D, printElementsMesh2D, readMeshFromSalomeDAT

"""
    Mesh2D_T

2D mesh main structure.
Basically contains 2 fields:
1. **nodes** - matrix of nodes (contains node's coordinates);
2. **elements** - matrix of elements (contains elements's nodes).
"""
mutable struct Mesh2D_T
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
        resultMesh.elements[i] = (blNum + (n + 1) + 1, blNum + (n + 1), blNum, blNum + 1)
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
    i = 1
    for element in mesh.elements
        print(i, " ")
        print(element, "\n")
        i += 1
    end
    print("\n")
end  # printElementsMesh2D

"""
    renumerateNodes(mesh::Mesh2D_T)

Renumerate nodes in each element according to FEM model.

# Arguments
- `mesh`: given model mesh.
"""
function renumerateNodes!(mesh::Mesh2D_T)
    for i in eachindex(mesh.elements)
        newNodes = (mesh.elements[i][2], mesh.elements[i][1], mesh.elements[i][4], mesh.elements[i][3])
        mesh.elements[i] = newNodes
    end
end

"""
    readMeshFromSalomeDAT(pathToFile::String)

Read mesh from .dat file generated via Salome platform.

# Arguments
- `pathToFile::String`: Path to given file;
- `type::meshType`: Type of given mesh.
"""
function readMeshFromSalomeDAT(pathToFile::String, type::meshType)
    mesh = Mesh2D_T(Vector{Tuple{Vararg{Float64}}}(undef, 0), Vector{Tuple{Vararg{Int}}}(undef, 0))
    open(pathToFile, "r") do file
        fileContents = split(read(file, String))  # Read once to provide more efficiency on big meshes
        nOfNodes = parse(Int, fileContents[1])
        nOfEntities = parse(Int, fileContents[2])
        nodes = Vector{Tuple{Vararg{Float64}}}(undef, nOfNodes)
        for i in 1:nOfNodes
            nodeIndex = parse(Int, fileContents[(i - 1) * 4 + 3])
            nodeX = parse(Float64, fileContents[(i - 1) * 4 + 4])
            nodeY = parse(Float64, fileContents[(i - 1) * 4 + 5])
            nodeZ = parse(Float64, fileContents[(i - 1) * 4 + 6])
            nodes[nodeIndex] = (nodeX, nodeY, nodeZ)
        end
        if type == Quad4Pts2D
            salomeTypeId = "204"
            nodesPerElement = 4
        else
            println("Given mesh type is not supported")
        end
        nOfElements = size(findall(isequal(salomeTypeId), fileContents))[1]
        elements = Vector{Tuple{Vararg{Float64}}}()
        elemPos = findfirst(isequal(salomeTypeId), fileContents) - 1
        while findnext(isequal(salomeTypeId), fileContents, elemPos) !== nothing
            elemNodes = Tuple{Vararg{Float64}}
            nextIdPos = findnext(isequal(salomeTypeId), fileContents, elemPos + 2)
            if nextIdPos !== nothing
                lastNode = nextIdPos - 2
            else
                lastNode = size(fileContents)[1]
            end
            nodesArray = Vector{Float64}()
            for i in (elemPos + 2):lastNode
                push!(nodesArray, parse(Float64, fileContents[i]))
            end
            elemNodes = Tuple(nodesArray)
            push!(elements, elemNodes)
            elemPos = lastNode + 1
        end
        mesh = Mesh2D_T(nodes, elements)
        renumerateNodes!(mesh)
        return mesh
    end
end  # readMeshFromSalomeDAT
