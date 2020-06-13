include("LoadVars.jl")
include("InputVars.jl")

export applyFixedX, applyFixedY, applyFixedXY, elementLoad

"""
    applyFixedX(node::Int, loads::Array, globalK::Array)

Apply fixing by OX axis.

# Arguments
- `node::Int`: node to which fixing should be applied;
- `loads::Array`: global loads vector;
- `globalK::Array`: global stiffness matrix.
"""
function applyFixedX(node::Int, loads::Array, globalK::Array)
    loads[2 * node - 1] = 0
    for col in 1:size(globalK)[2]
        globalK[2 * node - 1, col] = 0
    end
    for row in 1:size(globalK)[1]
        globalK[row, 2 * node - 1] = 0
    end
    globalK[2 * node - 1, 2 * node - 1] = 1
end  # applyFixedX

"""
    applyFixedY(node::Int, loads::Array, globalK::Array)

Apply fixing by OY axis.

# Arguments
- `node::Int`: node to which fixing should be applied;
- `loads::Array`: global loads vector;
- `globalK::Array`: global stiffness matrix.
"""
function applyFixedY(node::Int, loads::Array, globalK::Array)
    loads[2 * node] = 0
    for col in 1:size(globalK)[2]
        globalK[2 * node, col] = 0
    end
    for row in 1:size(globalK)[1]
        globalK[row, 2 * node] = 0
    end
    globalK[2 * node, 2 * node] = 1
end  # applyFixedY


"""
    applyFixedXY(node::Int, loads::Array, globalK::Array)

Apply fixing by both of OX and OY axis.

# Arguments
- `node::Int`: node to which fixing should be applied;
- `loads::Array`: global loads vector;
- `globalK::Array`: global stiffness matrix.
"""
function applyFixedXY(node::Int, loads::Array, globalK::Array)
    applyFixedX(node, loads, globalK)
    applyFixedY(node, loads, globalK)
end  # applyFixedXY

"""
    elementLoad(inputNodes::Array, pars::processPars, inputLoad::Array, loadDirect::loadDirection)

Return local load vector for given element.

# Arguments
- 
- `inputNodes::Array`: given nodes to which load should be applied;
- `pars::processPars`: parameters of current model;
- `inputLoad::Array`: given load;
- `loadDirect::loadDirection`: direction of given load.
"""
function elementLoad(elementNum::Int, pars::processPars, inputLoad::Array, loadDirect::loadDirection)
    nodesPerElement = length(pars.mesh.elements[elementNum])
    println("Nodes per element: ", nodesPerElement)
    xCoords = [pars.mesh.nodes[pars.mesh.elements[elementNum][i]][1] for i in 1:nodesPerElement]
    yCoords = [pars.mesh.nodes[pars.mesh.elements[elementNum][i]][2] for i in 1:nodesPerElement]
    load = inputLoad
    IntegrationOrder = 4
    FIntegrate(x) = 0
    if loadDirect == top
        FIntegrate(r) = transpose(Quad8Pts.displInterpMatr(r, 1)) * load * Quad8Pts.DetJs(r, 1, xCoords, yCoords)
    elseif loadDirect == left
        FIntegrate(s) = transpose(Quad8Pts.displInterpMatr(-1, s)) * load * Quad8Pts.DetJs(-1, s, xCoords, yCoords)
    elseif loadDirect == bottom
        FIntegrate(r) = transpose(Quad8Pts.displInterpMatr(r, -1)) * load * Quad8Pts.DetJs(r, -1, xCoords, yCoords)
    elseif loadDirect == right
        FIntegrate(s) = transpose(Quad8Pts.displInterpMatr(1, s)) * load * Quad8Pts.DetJs(1, s, xCoords, yCoords)
    else
        println("Given load direction is not supported")
        return nothing
    end
    F = multipleIntegral.gauss1DMethodMatrix(FIntegrate, IntegrationOrder)
    return F
end  # elementLoad

"""
    assemblyLoads(pars::processPars)

Assemble right part of linear system of equations. This method applies given local load vector to global ensemble.

# Arguments
- `pars::processPars`: parameters of current model.
"""
function assemblyLoads(pars::processPars)
    loadsVector = zeros(Float64, 2 * size(pars.mesh.nodes)[1])
    for (element, load) in pars.load
        elNum = element[1]
        direction = loadDirection(element[2])
        F = elementLoad(elNum, pars, load, direction)
        println("Element force: ", F)
        loadLocalNodes = Quad8Pts.nodesFromDirection(Int(direction))
        if (size(element)[1] - 2 != size(loadLocalNodes)[1])
            println("Incorrect input load")
            return nothing
        end
        loadGlobalNodes = [element[i] for i in 3:size(element)[1]]
        println("Load nodes for element ", elNum, ": ", loadGlobalNodes)
        for i in eachindex(loadLocalNodes)
            localIndex = loadLocalNodes[i]
            # globalIndex = pars.mesh.elements[elNum][i]
            globalIndex = loadGlobalNodes[i]
            println("Assigning load to ", globalIndex, " node (local: ", localIndex, ")")
            loadsVector[2 * globalIndex - 1] += F[2 * localIndex - 1]
            loadsVector[2 * globalIndex] += F[2 * localIndex]
        end
    end
    return loadsVector
end  # constructLoads
