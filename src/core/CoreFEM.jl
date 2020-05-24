"""
    CoreFEM
Module describing FEM core.
"""
module CoreFEM

include("MaterialVars.jl")
include("Load.jl")
include("Input.jl")
include("StiffnessMatrix.jl")

using MeshFEM
using DelimitedFiles
using BaseInterface

export fem2D

"""
    assemblyFEM2D(pars::processPars, targetMatrix::Array, currentElementMatrix::Array, elementNum::Number)

Assemble left part of linear system of equations. This method applies given local stiffness matrix to global ensemble.

# Arguments
- `pars::processPars`: parameters of current model;
- `targetMatrix::Array`: global stiffness matrix that should be updated;
- `currentElementMatrix::Array`: local stiffness matrix of given element that should be applid to global ensemble;
- `elementNum::Number`: number of given element in current mesh.
"""
function assemblyFEM2D(pars::processPars, targetMatrix::Array, currentElementMatrix::Array, elementNum::Number)
    elementNodes = pars.mesh.elements[elementNum]
    for i in eachindex(elementNodes)
        for j in eachindex(elementNodes)
            targetMatrix[2 * elementNodes[i] - 1, 2 * elementNodes[j] - 1] += currentElementMatrix[2 * i - 1, 2 * j - 1]
            targetMatrix[2 * elementNodes[i] - 1, 2 * elementNodes[j]] += currentElementMatrix[2 * i - 1, 2 * j]
            targetMatrix[2 * elementNodes[i], 2 * elementNodes[j] - 1] += currentElementMatrix[2 * i, 2 * j - 1]
            targetMatrix[2 * elementNodes[i], 2 * elementNodes[j]] += currentElementMatrix[2 * i, 2 * j]
        end
    end
end  # assemblyFEM2D

"""
    assemblyLoads(pars::processPars)

Assemble right part of linear system of equations. This method applies given local load vector to global ensemble.

# Arguments
- `pars::processPars`: parameters of current model.
"""
function assemblyLoads(pars::processPars)
    loadsVector = zeros(Real, 2 * size(pars.mesh.nodes)[1])
    # TODO: Make this method universal for any mesh and load.
    # for node in [5, 10, 15, 20, 25]
        F = elementLoad(5, pars)
        loadsVector[2 * 6 - 1] += F[3]
        loadsVector[2 * 6] += F[4]

        loadsVector[2 * 12 - 1] += F[5]
        loadsVector[2 * 12] += F[6]

        F = elementLoad(10, pars)
        loadsVector[2 * 12 - 1] += F[3]
        loadsVector[2 * 12] += F[4]

        loadsVector[2 * 18 - 1] += F[5]
        loadsVector[2 * 18] += F[6]

        F = elementLoad(15, pars)
        loadsVector[2 * 18 - 1] += F[3]
        loadsVector[2 * 18] += F[4]

        loadsVector[2 * 24 - 1] += F[5]
        loadsVector[2 * 24] += F[6]

        F = elementLoad(20, pars)
        loadsVector[2 * 24 - 1] += F[3]
        loadsVector[2 * 24] += F[4]

        loadsVector[2 * 30 - 1] += F[5]
        loadsVector[2 * 30] += F[6]

        F = elementLoad(25, pars)
        loadsVector[2 * 30 - 1] += F[3]
        loadsVector[2 * 30] += F[4]

        loadsVector[2 * 36 - 1] += F[5]
        loadsVector[2 * 36] += F[6]

    # end
    # for elementNum in eachindex(pars.mesh.elements)
    #     element = pars.mesh.elements[elementNum]
    #     H = elementLoad(elementNum)
    # end
    return loadsVector
end  # constructLoads

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
solve(globalK::Array, loadVector::Array) = globalK \ loadVector

"""
    fem2D()

Start calculation with test model.
"""
function fem2D()
    parameters = processPars(testMaterialProperties(), testBC(), testLoad(), generateTestMesh2D(5))
    nu = parameters.materialProperties[poisC]
    E = parameters.materialProperties[youngMod]
    elasticityMatrix = [1 nu 0; nu 1 0; 0 0 ((1 - nu) / 2)]
    elasticityMatrix *= E / (1 - nu ^ 2)
    ensembleMatrix = zeros(Real, 2 * size(parameters.mesh.nodes)[1], 2 * size(parameters.mesh.nodes)[1])
    for elementNum in eachindex(parameters.mesh.elements)
        K = stiffnessMatrix(elasticityMatrix, parameters, elementNum)
        assemblyFEM2D(parameters, ensembleMatrix, K, elementNum)
    end
    loadVector = assemblyLoads(parameters)
    applyConstraints(parameters, loadVector, ensembleMatrix)
    # Writing left part to file
    open("equation/K", "w") do file
        writedlm(file, ensembleMatrix)
    end
    # Writing right part to file
    open("equation/F", "w") do file
        writedlm(file, loadVector)
    end
    result = solve(ensembleMatrix, loadVector)
    # Writing result to file
    open("equation/result", "w") do file
        writedlm(file, result)
    end
    # Exporting results to CSV file
    BaseInterface.exportToCSV(result, parameters)
    # Exporting results to DAT file
    BaseInterface.exportToDAT2D(result, parameters)
    return result
end  # fem2D

end  # Core