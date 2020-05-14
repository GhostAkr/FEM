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
    element = pars.mesh.elements[elementNum]
    for i in eachindex(element)
        for j in eachindex(element)
            targetMatrix[2 * element[i] - 1, 2 * element[j] - 1] += currentElementMatrix[2 * i - 1, 2 * j - 1]
            targetMatrix[2 * element[i] - 1, 2 * element[j]] += currentElementMatrix[2 * i - 1, 2 * j]
            targetMatrix[2 * element[i], 2 * element[j] - 1] += currentElementMatrix[2 * i, 2 * j - 1]
            targetMatrix[2 * element[i], 2 * element[j]] += currentElementMatrix[2 * i, 2 * j]
        end
    end
end  # assemblyFEM2D

function constructLoads(pars::processPars)
    loads = pars.load
    loadsVector = zeros(Real, 2 * size(pars.mesh.nodes)[1])
    for (node, load) in pars.load
        loadsVector[2 * node - 1] = load[1]
        loadsVector[2 * node] = load[2]
    end
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
            println("Element #", elementNum, "; Nodes: ", parameters.mesh.elements[elementNum])
            println(K)
        end
        assemblyFEM2D(parameters, ensembleMatrix, K, elementNum)
    end
    loadVector = constructLoads(parameters)
    applyConstraints(parameters, loadVector, ensembleMatrix)
    open("equation/K", "w") do file
        writedlm(file, ensembleMatrix)
    end
    open("equation/F", "w") do file
        writedlm(file, loadVector)
    end
    result = solve(ensembleMatrix, loadVector)
    open("equation/result", "w") do file
        writedlm(file, result)
    end
    BaseInterface.exportToCSV(result, parameters)
    return result
end  # fem2D

end  # Core