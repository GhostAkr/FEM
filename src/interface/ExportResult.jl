using CoreFEM
using MeshFEM

# Method assumes that result is for 2D object
"""
    exportToCSV(result::Array, pars::processPars)

Export given result to CSV file.

# Arguments
- `result::Array`: results vector (assuming that given result is for some 2D object);
- `pars::processPars`: parameters of current model.
"""
function exportToCSV(result::Array, pars::processPars)
    displacementHeader = "X, Y, Z, Displacement\n"
    open("output/displacementX.csv", "w") do file
        write(file, displacementHeader)
        for i in 1:size(pars.mesh.nodes)[1]
            write(file, string(pars.mesh.nodes[i][1]) * ", " * string(pars.mesh.nodes[i][2]) * ", " * "0, " * string(result[2 * i - 1]) * "\n")
        end
    end
    open("output/displacementY.csv", "w") do file
        write(file, displacementHeader)
        for i in 1:size(pars.mesh.nodes)[1]
            write(file, string(pars.mesh.nodes[i][1]) * ", " * string(pars.mesh.nodes[i][2]) * ", " * "0, " * string(result[2 * i]) * "\n")
        end
    end
end  # exportToCSV

function exportToDAT2D(result::Array, pars::processPars)
    open("output/displacementX.dat", "w") do file
        for i in 1:size(pars.mesh.nodes)[1]
            write(file, string(pars.mesh.nodes[i][1]) * " " * string(pars.mesh.nodes[i][2]) * " " * string(result[2 * i - 1]) * "\n")
        end
    end
    open("output/displacementY.dat", "w") do file
        for i in 1:size(pars.mesh.nodes)[1]
            write(file, string(pars.mesh.nodes[i][1]) * " " * string(pars.mesh.nodes[i][2]) * " " * string(result[2 * i]) * "\n")
        end
    end
end # exportToDAT2D

function vtkCellsListSize(pars::processPars)
    cellsListSize = 0
    for element in pars.mesh.elements
        cellsListSize += 1
        for node in element
            cellsListSize += 1
        end
    end
    return cellsListSize
end  # vtkCellsListSize

function exportToVTK(result::Array, deformations, stresses, vonMises, pars::processPars, type::meshType)
    vtkIdentifier = "# vtk DataFile Version 3.0\n"
    vtkFormat = "ASCII\n"
    vtkDataSetKeyword = "DATASET"
    nOfNodes = size(pars.mesh.nodes)[1]
    nOfElements = size(pars.mesh.elements)[1]
    cellsListSize = vtkCellsListSize(pars)
    vtkPointsKeyword = "POINTS"
    vtkPointsType = "float"
    vtkCellsKeyword = "CELLS"
    vtkCellTypesKeyword = "CELL_TYPES"
    if type === Quad4Pts2D
        cellType = "9"  # 4-nodes elements
    elseif type === Quad8Pts2D
        cellType = "23"  # 8-nodes elements
    elseif type === Iso8Pts3DMeshType
        cellType = "12"  # 8-nodes 3D elements
    else
        println("Unknown mesh type while exporting to VTK")
        return
    end
    vtkPointDataKeyword = "POINT_DATA"
    vtkScalarKeyword = "SCALARS"
    vtkTableKeyword = "LOOKUP_TABLE"
    vtkTableDefault = "default"
    vtkVectorsKeyword = "VECTORS"
    open("output/Results.vtk", "w") do file
        # Writing headers
        vtkHeader = "FEM Results\n"
        write(file, vtkIdentifier)
        write(file, vtkHeader)
        write(file, vtkFormat)
        # Writing datasets
        vtkDataSetType = "UNSTRUCTURED_GRID"
        write(file, vtkDataSetKeyword * " " * vtkDataSetType * "\n")
        write(file, vtkPointsKeyword * " " * string(nOfNodes) * " " * vtkPointsType * "\n")
        for node in pars.mesh.nodes
            write(file, string(node[1]) * " " * string(node[2]) * " " * string(node[3]) * "\n")
        end
        write(file, "\n")
        write(file, vtkCellsKeyword * " " * string(nOfElements) * " " * string(cellsListSize) * "\n")
        for element in pars.mesh.elements
            nOfElemNodes = length(element)
            write(file, string(nOfElemNodes))
            for node in element
                write(file, " " * string(node - 1))
            end
            write(file, "\n")
        end
        write(file, "\n")
        write(file, vtkCellTypesKeyword * " " * string(nOfElements) * "\n")
        for element in pars.mesh.elements
            write(file, cellType * "\n")
        end
        write(file, "\n")

        write(file, vtkPointDataKeyword * " " * string(nOfNodes) * "\n")
        displacementsKeyword = "Displacements"
        write(file, vtkVectorsKeyword * " " * displacementsKeyword * " " * vtkPointsType * "\n")  # By defalt it uses 1 vector per point
        if type === Iso8Pts3DMeshType  # TODO: Make this check more scalable
            for nodeIndex in 1:nOfNodes
                write(file, string(result[3 * nodeIndex - 2]) * " " * string(result[3 * nodeIndex - 1]) * " " * string(result[3 * nodeIndex]) * "\n")
            end
        else
            for nodeIndex in 1:nOfNodes
                write(file, string(result[2 * nodeIndex - 1]) * " " * string(result[2 * nodeIndex]) * " 0.0\n")
            end
        end
        write(file, "\n")

        if !isnothing(deformations)
            deformations_cnt = size(deformations)[2]
            deformationsKeyword = "Deformations"
            write(file, vtkScalarKeyword * " " * deformationsKeyword * " " * vtkPointsType * 
                " " * string(deformations_cnt) * "\n" * vtkTableKeyword * " " * 
                vtkTableDefault * "\n")
            for nodeIndex in 1:nOfNodes
                node_def = ""
                for defidx in 1:deformations_cnt
                    node_def *= string(deformations[nodeIndex, defidx])
                    if defidx == deformations_cnt
                        break
                    end
                    node_def *= " "
                end
                node_def *= "\n"
                write(file, node_def)
            end
            write(file, "\n")
        end

        if !isnothing(stresses)
            stresses_cnt = size(deformations)[2]
            stressesKeyword = "Stresses"
            write(file, vtkScalarKeyword * " " * stressesKeyword * " " * vtkPointsType * 
            " " * string(stresses_cnt) * "\n" * vtkTableKeyword * " " * vtkTableDefault * 
            "\n")
            for nodeIndex in 1:nOfNodes
                node_stress = ""
                for stressidx in 1:stresses_cnt
                    node_stress *= string(stresses[nodeIndex, stressidx])
                    if stressidx == stresses_cnt
                        break
                    end
                    node_stress *= " "
                end
                node_stress *= "\n"
                write(file, node_stress)
            end
            write(file, "\n")
        end

        if !isnothing(vonMises)
            vonMisesKeyword = "VonMises"
            write(file, vtkScalarKeyword * " " * vonMisesKeyword * " " * vtkPointsType * "\n")  # By defalt it uses 1 scalar per point
            write(file, vtkTableKeyword * " " * vtkTableDefault * "\n")
            for nodeIndex in 1:nOfNodes
                write(file, string(vonMises[nodeIndex]) * "\n")
            end
        end
    end
end  # exportToVTK