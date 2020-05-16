module CoreFEM

include("Material.jl")
include("Load.jl")
include("Input.jl")
include("StiffnessMatrix.jl")

using MeshFEM
using DelimitedFiles
using BaseInterface

export fem2D

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

function assemblyLoads(pars::processPars)
    loadsVector = zeros(Real, 2 * size(pars.mesh.nodes)[1])
    # for (node, load) in pars.load
    #     loadsVector[2 * node - 1] = load[1]
    #     loadsVector[2 * node] = load[2]
    # end
    # for node in [5, 10, 15, 20, 25]
        F = elementLoad(2, pars)  # TODO: In this simple case loads vectors will be the same for all elements, so we just use 1. Need to make it depending on element
        loadsVector[2 * 3 - 1] += F[7]
        loadsVector[2 * 3] += F[8]

        loadsVector[2 * 6 - 1] += F[1]
        loadsVector[2 * 6] += F[2]

        F = elementLoad(4, pars)
        loadsVector[2 * 6 - 1] += F[7]
        loadsVector[2 * 6] += F[8]

        loadsVector[2 * 9 - 1] += F[1]
        loadsVector[2 * 9] += F[2]

        # F = elementLoad(12, pars)
        # loadsVector[2 * 15 - 1] += F[7]
        # loadsVector[2 * 15] += F[8]

        # loadsVector[2 * 20 - 1] += F[1]
        # loadsVector[2 * 20] += F[2]

        # F = elementLoad(16, pars)
        # loadsVector[2 * 20 - 1] += F[7]
        # loadsVector[2 * 20] += F[8]

        # loadsVector[2 * 25 - 1] += F[1]
        # loadsVector[2 * 25] += F[2]
    # end
    # for elementNum in eachindex(pars.mesh.elements)
    #     element = pars.mesh.elements[elementNum]
    #     H = elementLoad(elementNum)
    # end
    return loadsVector
end  # constructLoads

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

# Since Julia provides native workaround with linear algebra elements the best solution in most cases
# is to use standart syntax to solve linear equations system.
# According to official documentation Julia will choose the best solving method by itself.
# If it's not, there should be a way to control it.
solve(globalK::Array, loadVector::Array) = globalK \ loadVector

function fem2D()
    parameters = processPars(testMaterialProperties(), testBC(), testLoad(), generateTestMesh2D())
    printProcessPars(parameters)
    nu = parameters.materialProperties[poisC]
    E = parameters.materialProperties[youngMod]
    println("Nu = ", nu)
    println("E = ", E)
    elasticityMatrix = [1 nu 0; nu 1 0; 0 0 ((1 - nu) / 2)]
    elasticityMatrix *= E / (1 - nu ^ 2)
    println("C:")
    println(elasticityMatrix)
    ensembleMatrix = zeros(Real, 2 * size(parameters.mesh.nodes)[1], 2 * size(parameters.mesh.nodes)[1])
    for elementNum in eachindex(parameters.mesh.elements)
        K = stiffnessMatrix(elasticityMatrix, parameters, elementNum)
        if elementNum == 1
            println("K:")
            println(K)
            println()
        end
        assemblyFEM2D(parameters, ensembleMatrix, K, elementNum)
    end
    loadVector = assemblyLoads(parameters)
    applyConstraints(parameters, loadVector, ensembleMatrix)
    open("equation/K", "w") do file
        writedlm(file, ensembleMatrix)
    end
    open("equation/F", "w") do file
        writedlm(file, loadVector)
    end
    println("Size of K is ", size(ensembleMatrix))
    println("Size of F is ", size(loadVector), "\n")
    println("11th row in K: ")
    for i in 1:size(ensembleMatrix)[1]
        print(ensembleMatrix[11, i], "; ")
    end
    println()
    println("11th element in F: ", loadVector[11], "\n")
    result = solve(ensembleMatrix, loadVector)
    open("equation/result", "w") do file
        writedlm(file, result)
    end
    BaseInterface.exportToCSV(result, parameters)
    return result
end  # fem2D

end  # Core