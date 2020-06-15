using CoreFEM

include("ExportVars.jl")

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

function exportToVTK(result::Array, deformations::Array, stresses::Array, vonMises::Array,  pars::processPars)
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
    cellType = "23"  # 8-nodes elements
    vtkPointDataKeyword = "POINT_DATA"
    vtkScalarKeyword = "SCALARS"
    vtkTableKeyword = "LOOKUP_TABLE"
    vtkTableDefault = "default"
    vtkVectorsKeyword = "VECTORS"
    open("output/displacementX.vtk", "w") do file
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
        for nodeIndex in 1:nOfNodes
            write(file, string(result[2 * nodeIndex - 1]) * " " * string(result[2 * nodeIndex]) * " 0.0\n")
        end
        write(file, "\n")
        deformationsKeyword = "Deformations"
        write(file, vtkVectorsKeyword * " " * deformationsKeyword * " " * vtkPointsType * "\n")  # By defalt it uses 1 scalar per point
        for nodeIndex in 1:nOfNodes
            # write(file, string(deformations[nodeIndex][1]) * "\n")
            write(file, string(deformations[nodeIndex][1]) * " " * string(deformations[nodeIndex][2]) * " " * string(deformations[nodeIndex][3]) * "\n")
        end
        write(file, "\n")
        stressesKeyword = "Stresses"
        write(file, vtkVectorsKeyword * " " * stressesKeyword * " " * vtkPointsType * "\n")  # By defalt it uses 1 scalar per point
        for nodeIndex in 1:nOfNodes
            # write(file, string(stresses[nodeIndex][1]) * "\n")
            write(file, string(stresses[nodeIndex][1]) * " " * string(stresses[nodeIndex][2]) * " " * string(stresses[nodeIndex][3]) * "\n")
        end
        write(file, "\n")
        vonMisesKeyword = "VonMises"
        write(file, vtkScalarKeyword * " " * vonMisesKeyword * " " * vtkPointsType * "\n")  # By defalt it uses 1 scalar per point
        write(file, vtkTableKeyword * " " * vtkTableDefault * "\n")
        for nodeIndex in 1:nOfNodes
            write(file, string(vonMises[nodeIndex]) * "\n")
        end
    end
end  # exportToVTK