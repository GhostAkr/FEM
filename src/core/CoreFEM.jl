"""
    CoreFEM
Module describing FEM core.
"""
module CoreFEM

include("MaterialVars.jl")
include("Load.jl")
include("Input.jl")
include("StiffnessMatrix.jl")
include("Deformations.jl")
include("Constants.jl")
include("Stresses.jl")

using MeshFEM
using DelimitedFiles
using BaseInterface
using IterativeSolvers
using ElementTypes
using TestFEM

using Quad4Pts
using Quad8Pts

export fem2D

"""
    assembly_left_part!(pars::processPars, targetMatrix::Array, currentElementMatrix::Array, elementNum::Number, freedom_deg::Int)

Assemble left part of linear system of equations. This method applies given local stiffness matrix to global ensemble.

# Arguments
- `pars::processPars`: parameters of current model;
- `targetMatrix::Array`: global stiffness matrix that should be updated;
- `currentElementMatrix::Array`: local stiffness matrix of given element that should be applid to global ensemble;
- `elementNum::Number`: number of given element in current mesh;
- `freedom_deg::Int`: degree of freedom.
"""
function assembly_left_part!(pars::processPars, targetMatrix::Array, currentElementMatrix::Array, elementNum::Number, freedom_deg::Int)
    elementNodes = pars.mesh.elements[elementNum]
    for i in eachindex(elementNodes)
        for j in eachindex(elementNodes)
            for first_offset in 1:freedom_deg
                for second_offset in 1:freedom_deg
                    targetMatrix[freedom_deg * elementNodes[i] - (first_offset - 1), freedom_deg * elementNodes[j] - (second_offset - 1)] += 
                                currentElementMatrix[freedom_deg * i - (first_offset - 1), freedom_deg * j - (second_offset - 1)]
                end
            end
        end
    end
end  # assemblyFEM2D

"""
    applyConstraints(pars::processPars, loads::Array, globalK::Array)

Applies constraints from model parameters to given ensemble.

# Arguments
- `pars::processPars`: parameters of current model;
- `loads::Array`: global loads vector (right part of equations system);
- `globalK::Array`: global stiffness matrix (left part of equations system).
"""
function applyConstraints(pars::processPars, loads::Array, globalK::Array)
    for (node, bc) in pars.bc
        if bc == fixedX
            applyFixedX(node, loads, globalK)
        elseif bc == fixedY
            applyFixedY(node, loads, globalK)
        elseif bc == fixedXY
            applyFixedXY(node, loads, globalK)
        else
            println("Unhandled boundary condition")
        end
    end
end  # applyConstraints

"""
    solve(globalK::Array, loadVector::Array)

Solve given equation system.

Since Julia provides native workaround with linear algebra elements the best solution in most cases
is to use standart syntax to solve linear equations system.
According to official documentation Julia will choose the best solving method by itself.
If it's not, there should be a way to control it.

# Arguments
- `globalK::Array`: global stiffness matrix (left part of equations system);
- `loadVector::Array`: global loads vector (right part of equations system).
"""
function solve(globalK::Array, loadVector::Array)
    initialVector = fill(1.0, size(loadVector)[1])
    # return minres!(initialVector, globalK, loadVector, tol = 1e-10)
    # println("In solver")
    return globalK \ loadVector
end

function defineElemType(elemTypeID::FETypes)
    resElement = nothing
    if elemTypeID === Quad4TypeID
        resElement = Quad4Type("Quad4Type")
    elseif (elemTypeID === Quad8TypeID)
        resElement = Quad8Type("Quad8Type")
    else
        println("Unknown finite element type")
    end
    return resElement
end

function typeMeshFromElement(elemTypeID::FETypes)
    resMeshType = nothing
    if elemTypeID === Quad4TypeID
        resMeshType = Quad4Pts2D
    elseif elemTypeID === Quad8TypeID
        resMeshType = Quad8Pts2D
    else
        println("Unknown element type while converting to mesh type")
    end
    return resMeshType
end

"""
    fem2D()

Start calculation with given model.

# Arguments
- `meshPath::String`: Path to given mesh;
- `dataPath::String`: Path to given initial data.
"""
function fem2D(meshPath::String, dataPath::String, elemTypeID::FETypes)
    freedom_deg = 2
    # Getting element type
    elementType = defineElemType(elemTypeID)
    if (elementType === nothing)
        println("Element type passed to fem2D() is unknown")
        return
    end

    # Getting mesh type
    meshType = typeMeshFromElement(elemTypeID)

    parameters = processPars(testMaterialProperties(), testBC(), testLoad(), generateTestMesh2D(2))
    read_params_JSON!(dataPath, parameters)
    # readParameters!(dataPath, parameters)
    printProcessPars(parameters)
    parameters.mesh = readMeshFromSalomeDAT(meshPath, meshType)
    intOrder = 3
    nu = parameters.materialProperties[poisC]
    E = parameters.materialProperties[youngMod]
    C = elasticityMatrix(E, nu, 2)  # 1 - plane strain; 2 - plain stress
    ensembleMatrix = zeros(Float64, 2 * size(parameters.mesh.nodes)[1], 2 * size(parameters.mesh.nodes)[1])
    for elementNum in eachindex(parameters.mesh.elements)
        K = stiffnessMatrix(C, parameters, elementNum, intOrder, elementType)
        assembly_left_part!(parameters, ensembleMatrix, K, elementNum, freedom_deg)
    end
    loadVector = assembly_loads!(parameters, intOrder, elementType, freedom_deg)
    applyConstraints(parameters, loadVector, ensembleMatrix)
    # Writing left part to file
    open("equation/K", "w") do file
        writedlm(file, ensembleMatrix)
    end
    # Writing right part to file
    open("equation/F", "w") do file
        writedlm(file, loadVector)
    end
    println("Solving...")
    result = solve(ensembleMatrix, loadVector)

    if TestFEM.verify_example(meshPath, dataPath, result)
        @info "Result is correct"
    else
        @info "Result is INcorrect"
    end

    deformations = calculateDeformations(result, parameters, intOrder, elementType)
    stresses = calculateStresses(deformations, C, parameters)
    vonMises = calculateVonMises(stresses)
    # Writing result to file
    open("equation/result", "w") do file
        writedlm(file, result)
    end
    # Exporting results to VTK file
    BaseInterface.exportToVTK(result, deformations, stresses, vonMises, parameters, meshType)
    return result
end  # fem2D

end  # Core