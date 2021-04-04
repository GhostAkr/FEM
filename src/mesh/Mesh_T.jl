include("MeshVars.jl")

using AsterReader

export Mesh2D_T, generateTestMesh2D, printNodesMesh2D, printElementsMesh2D, readMeshFromSalomeDAT
export generateTestMesh3D

# TODO: Rename mesh type and print functions as well as it can be used for 3D models too.

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
    "Groups of nodes"
    groups::Dict{String, Vector{Vector{Int}}}
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
    resultMesh = Mesh2D_T(Vector{Tuple{Vararg{Float64}}}(undef, nOfNodes), Vector{Tuple{Vararg{Int}}}(undef, nOfElements), Dict{String, Vector{Vector{Int}}}())
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

function generateTestMesh3D()
    nOfNodes = 12
    nOfElements = 2 
    resultMesh = Mesh2D_T(Vector{Tuple{Vararg{Float64}}}(undef, nOfNodes), Vector{Tuple{Vararg{Int}}}(undef, nOfElements))
    # Nodes
    resultMesh.nodes[1] = (0, 0, 0)
    resultMesh.nodes[2] = (50, 0, 0)
    resultMesh.nodes[3] = (100, 0, 0)
    resultMesh.nodes[4] = (0, 0, 100)
    resultMesh.nodes[5] = (50, 0, 100)
    resultMesh.nodes[6] = (100, 0, 100)
    resultMesh.nodes[7] = (0, 100, 0)
    resultMesh.nodes[8] = (50, 100, 0)
    resultMesh.nodes[9] = (100, 100, 0)
    resultMesh.nodes[10] = (0, 100, 100)
    resultMesh.nodes[11] = (50, 100, 100)
    resultMesh.nodes[12] = (100, 100, 100)
    # Elements
    resultMesh.elements[1] = (5, 11, 10, 4, 2, 8, 7, 1)
    resultMesh.elements[2] = (6, 12, 11, 5, 3, 9, 8, 2)
    return resultMesh
end  # generateTestMesh3D


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
function renumerateNodes!(mesh::Mesh2D_T, type::meshType)
    for i in eachindex(mesh.elements)
        x = [mesh.nodes[mesh.elements[i][1]][1], mesh.nodes[mesh.elements[i][2]][1], mesh.nodes[mesh.elements[i][3]][1], mesh.nodes[mesh.elements[i][4]][1]]
        y = [mesh.nodes[mesh.elements[i][1]][2], mesh.nodes[mesh.elements[i][2]][2], mesh.nodes[mesh.elements[i][3]][2], mesh.nodes[mesh.elements[i][4]][2]]
        newNodes = nothing
        if type === Quad4Pts2D
            # All possible renumerations
            newNodes = (mesh.elements[i][1], mesh.elements[i][4], mesh.elements[i][3], mesh.elements[i][2])
            # newNodes = (mesh.elements[i][4], mesh.elements[i][3], mesh.elements[i][2], mesh.elements[i][1])
            # newNodes = (mesh.elements[i][3], mesh.elements[i][2], mesh.elements[i][1], mesh.elements[i][4])
            # newNodes = (mesh.elements[i][2], mesh.elements[i][1], mesh.elements[i][4], mesh.elements[i][3])
        elseif type === Quad8Pts2D
            # TODO: This works not for every mesh generated in Salome.
            # Need to check how Salome actually numerates nodes for such finite element.

            # All possible renumerations
            newNodes = (mesh.elements[i][3], mesh.elements[i][2], mesh.elements[i][1], mesh.elements[i][4], mesh.elements[i][6], mesh.elements[i][5], mesh.elements[i][8], mesh.elements[i][7])
            # newNodes = (mesh.elements[i][2], mesh.elements[i][1], mesh.elements[i][4], mesh.elements[i][3], mesh.elements[i][5], mesh.elements[i][8], mesh.elements[i][7], mesh.elements[i][6])
            # newNodes = (mesh.elements[i][1], mesh.elements[i][4], mesh.elements[i][3], mesh.elements[i][2], mesh.elements[i][8], mesh.elements[i][7], mesh.elements[i][6], mesh.elements[i][5])
            # newNodes = (mesh.elements[i][4], mesh.elements[i][3], mesh.elements[i][2], mesh.elements[i][1], mesh.elements[i][7], mesh.elements[i][6], mesh.elements[i][5], mesh.elements[i][8])
        elseif type === Iso8Pts3DMeshType
            # TODO: Renumerate nodes in this case
        else
            println("Unknown mesh type while renumerating nodes")
        end
        if (newNodes !== nothing)
            mesh.elements[i] = newNodes
        end
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
    mesh = Mesh2D_T(Vector{Tuple{Vararg{Float64}}}(undef, 0), Vector{Tuple{Vararg{Int}}}(undef, 0), Dict{String, Vector{Vector{Int}}}())
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
        elseif type == Quad8Pts2D
            salomeTypeId = "208"
            nodesPerElement = 8
        elseif type == Iso8Pts3DMeshType
            # TODO: Handle this case
        else
            println("Given mesh type is not supported")
        end
        elements = Vector{Tuple{Vararg{Float64}}}()
        index = 4 * nOfNodes + 3
        while true
            if index > size(fileContents)[1]
                break
            end
            index += 1
            if fileContents[index] == "102"
                index += 2
            elseif fileContents[index] == "103"
                index += 3
            elseif fileContents[index] == "204"
                nodesArray = Vector{Float64}()
                for i in (index + 1):(index + 4)
                    push!(nodesArray, parse(Float64, fileContents[i]))
                end
                elemNodes = Tuple(nodesArray)
                push!(elements, elemNodes)
                index += 4
            elseif fileContents[index] == "208"
                nodesArray = Vector{Float64}()
                for i in (index + 1):(index + 8)
                    push!(nodesArray, parse(Float64, fileContents[i]))
                end
                elemNodes = Tuple(nodesArray)
                push!(elements, elemNodes)
                index += 8
            end
            index += 1
        end
        groups = Dict{String, Vector{Vector{Int}}}()
        mesh = Mesh2D_T(nodes, elements, groups)
        x0 = []
        y0 = []
        epsNull = 1e-10
        for i in eachindex(mesh.nodes)
            if abs(mesh.nodes[i][1]) < epsNull
                push!(x0, i)
            end
            if abs(mesh.nodes[i][2]) < epsNull
                push!(y0, i)
            end
        end
        println("X = 0 points: ", x0)
        println("Y = 0 points: ", y0)
        renumerateNodes!(mesh, type)
        return mesh
    end
end  # readMeshFromSalomeDAT

"""
    nodes_group_by_name(mesh::Mesh2D_T, name::String)

Get nodes group by name from given mesh.

# Arguments
- `mesh::Mesh2D_T`: Mesh which contains target groups;
- `name::String`: Name of target group.
"""
function nodes_group_by_name(mesh::Mesh2D_T, name::String)
    if !(haskey(mesh.groups, name))
        @warn("Nodes group with such key not found")
        return nothing
    end
    targetGroup = mesh.groups[name]
    return targetGroup
end  # nodes_group_by_name

"""
    add_nodes_group!(mesh::Mesh2D_T, name::String, nodes::Array{Tuple{Vararg{Int}}})

Add nodes group to given mesh.

# Arguments
- `mesh::Mesh2D_T`: Mesh to which group should be added;
- `name::String`: Name of target group;
- `nodes::Array{Tuple{Vararg{Int}}}`: Array of nodes in target group.
"""
function add_nodes_group!(mesh::Mesh2D_T, name::String, nodes::Vector{Vector{Int}})
    if haskey(mesh.groups, name)
        @warn("Group with given name already exists, replacing to the new one")
    end
    mesh.groups[name] = nodes
end  # add_nodes_group!

"""
    add_nodes_multiple_groups!(mesh::Mesh2D_T, groups::Dict{String, Vector{Tuple{Vararg{Int}}}})

Add multiple nodes groups to given mesh.

# Arguments
- `mesh::Mesh2D_T`: Mesh to which group should be added;
- `groups::Dict{String, Vector{Tuple{Vararg{Int}}}}`: Dictionary of groups which should be added to given mesh.
"""
function add_nodes_multiple_groups!(mesh::Mesh2D_T, groups::Dict{String, Vector{Vector{Int}}})
    for key in keys(groups)
        mesh.groups[key] = groups[key]
    end
end  # add_nodes_multiple_groups!

"""
    delete_nodes_group_by_name!(mesh::Mesh2D_T, name::String)

Delete nodes group from given mesh by name.

# Arguments
- `mesh::Mesh2D_T`: Mesh to which group should be added;
- `name::String`: Name of target group.
"""
function delete_nodes_group_by_name!(mesh::Mesh2D_T, name::String)
    delete!(mesh.groups, name)
end  # delete_nodes_group_by_name!

"""
    print_nodes_groups(mesh::Mesh2D_T)

Print all nodes groups of given mesh to the REPL.

# Arguments
- `mesh::Mesh2D_T`: Mesh to which group should be added.
"""
function print_nodes_groups(mesh::Mesh2D_T)
    println(mesh.groups)
end  # print_nodes_groups

"""
    read_nodes_groups_from_salome(salome_mesh::String)

Read groups of nodes from Salome MED file. This method uses
AsterReader module to get access to MED file.

# Arguments
- `salome_mesh::String`: Path to source MED file.
"""
function read_nodes_groups_from_salome(salome_mesh::String)
    mesh = AsterReader.aster_read_mesh(salome_mesh)

    # Getting elements from mesh
    elements = mesh["elements"]
    # Getting groups from mesh
    element_sets = mesh["element_sets"]

    res_dict = Dict{String, Vector{Vector{Int}}}()
    for name in keys(element_sets)
        set_nodes = Vector{Vector{Int}}()
        elem_set = element_sets[name]
        for elem_num in elem_set
            nodes = elements[elem_num]
            push!(set_nodes, nodes)
        end
        res_dict[name] = set_nodes
    end

    return res_dict
end  # read_nodes_groups_from_salome
