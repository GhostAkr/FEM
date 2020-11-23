using Quad4Pts
using Quad8Pts

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

# TODO: Unify with `applyFixedX`
function applyFixedX3D(node::Int, loads::Array, globalK::Array)
    loads[3 * node - 2] = 0
    for col in 1:size(globalK)[2]
        globalK[3 * node - 2, col] = 0
    end
    for row in 1:size(globalK)[1]
        globalK[row, 3 * node - 2] = 0
    end
    globalK[3 * node - 2, 3 * node - 2] = 1
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

# TODO: Unify with `applyFixedY`
function applyFixedY3D(node::Int, loads::Array, globalK::Array)
    loads[3 * node - 1] = 0
    for col in 1:size(globalK)[2]
        globalK[3 * node - 1, col] = 0
    end
    for row in 1:size(globalK)[1]
        globalK[row, 3 * node - 1] = 0
    end
    globalK[3 * node - 1, 3 * node - 1] = 1
end  # applyFixedY

# TODO: Make an unparsed case for 2D model
function applyFixedZ3D(node::Int, loads::Array, globalK::Array)
    loads[3 * node] = 0
    for col in 1:size(globalK)[2]
        globalK[3 * node, col] = 0
    end
    for row in 1:size(globalK)[1]
        globalK[row, 3 * node] = 0
    end
    globalK[3 * node, 3 * node] = 1
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
function elementLoad(elementNum::Int, pars::processPars, inputLoad::Array, loadDirect::loadDirection, intOrder::Int, elemTypeInd::FiniteElement)
    nodesPerElement = length(pars.mesh.elements[elementNum])
    xCoords = [pars.mesh.nodes[pars.mesh.elements[elementNum][i]][1] for i in 1:nodesPerElement]
    yCoords = [pars.mesh.nodes[pars.mesh.elements[elementNum][i]][2] for i in 1:nodesPerElement]
    load = inputLoad
    # IntegrationOrder = 4
    FIntegrate(x) = 0
    if loadDirect == top
        FIntegrate(r) = transpose(displInterpMatr(r, 1, elemTypeInd)) * load * DetJs(r, 1, xCoords, yCoords, elemTypeInd)
    elseif loadDirect == left
        FIntegrate(s) = transpose(displInterpMatr(-1, s, elemTypeInd)) * load * DetJs(-1, s, xCoords, yCoords, elemTypeInd)
    elseif loadDirect == bottom
        FIntegrate(r) = transpose(displInterpMatr(r, -1, elemTypeInd)) * load * DetJs(r, -1, xCoords, yCoords, elemTypeInd)
    elseif loadDirect == right
        FIntegrate(s) = transpose(displInterpMatr(1, s, elemTypeInd)) * load * DetJs(1, s, xCoords, yCoords, elemTypeInd)
    else
        println("Given load direction is not supported")
        return nothing
    end
    F = multipleIntegral.gauss1DMethodMatrix(FIntegrate, intOrder)
    return F
end  # elementLoad

function elementLoad3D(elementNum::Int, pars::processPars, inputLoad::Array, loadDirect::loadDirection, intOrder::Int, elemTypeInd::FiniteElement)
    nodesPerElement = length(pars.mesh.elements[elementNum])
    xCoords = [pars.mesh.nodes[pars.mesh.elements[elementNum][i]][1] for i in 1:nodesPerElement]
    yCoords = [pars.mesh.nodes[pars.mesh.elements[elementNum][i]][2] for i in 1:nodesPerElement]
    yCoords = [pars.mesh.nodes[pars.mesh.elements[elementNum][i]][3] for i in 1:nodesPerElement]
    load = inputLoad
    # IntegrationOrder = 4
    FIntegrate(x) = 0
    if loadDirect == top
        FIntegrate(r, s) = transpose(displInterpMatr(r, s, 1 elemTypeInd)) * load * DetJs(r, s, 1, xCoords, yCoords, zCoords, elemTypeInd)
    elseif loadDirect == left
        FIntegrate(r, t) = transpose(displInterpMatr(r, -1, t, elemTypeInd)) * load * DetJs(r, -1, t, xCoords, yCoords, zCoords, elemTypeInd)
    elseif loadDirect == bottom
        FIntegrate(r, s) = transpose(displInterpMatr(r, s, 1, elemTypeInd)) * load * DetJs(r, s, -1, xCoords, yCoords, zCoords, elemTypeInd)
    elseif loadDirect == right
        FIntegrate(r, t) = transpose(displInterpMatr(r, 1, t, elemTypeInd)) * load * DetJs(r, 1, t, xCoords, yCoords, zCoords, elemTypeInd)
    elseif loadDirect == backwards  # To us
        FIntegrate(s, t) = transpose(displInterpMatr(1, s, t, elemTypeInd)) * load * DetJs(1, s, t, xCoords, yCoords, zCoords, elemTypeInd)
    elseif loadDirect == towards  # From us
        FIntegrate(s, t) = transpose(displInterpMatr(-1, s, t, elemTypeInd)) * load * DetJs(-1, s, t, xCoords, yCoords, zCoords, elemTypeInd)
    else
        println("Given load direction is not supported")
        return nothing
    end
    F = multipleIntegral.gaussMethodMatrix(FIntegrate, intOrder)

    return F
end  # elementLoad

"""
    assemblyLoads(pars::processPars)

Assemble right part of linear system of equations. This method applies given local load vector to global ensemble.

# Arguments
- `pars::processPars`: parameters of current model.
"""
function assemblyLoads(pars::processPars, intOrder::Int, elemTypeInd::FiniteElement)
    loadsVector = zeros(Float64, 2 * size(pars.mesh.nodes)[1])
    for (element, load) in pars.load
        elNum = element[1]
        direction = loadDirection(element[2])
        F = elementLoad(elNum, pars, load, direction, intOrder, elemTypeInd)
        loadLocalNodes = nodesFromDirection(Int(direction), elemTypeInd)
        # println("Size: ", size(element)[1])
        if (size(element)[1] - 2 != size(loadLocalNodes)[1])
            println("Incorrect input load")
            return nothing
        end
        loadGlobalNodes = [element[i] for i in 3:size(element)[1]]
        for i in eachindex(loadLocalNodes)
            localIndex = loadLocalNodes[i]
            globalIndex = loadGlobalNodes[i]
            loadsVector[2 * globalIndex - 1] += F[2 * localIndex - 1]
            loadsVector[2 * globalIndex] += F[2 * localIndex]
            # println("Loaded node: ", globalIndex)
        end
    end
    # println(loadsVector)
    return loadsVector
end  # assemblyLoads

function assemblyLoads3D(pars::processPars, intOrder::Int, elemTypeInd::FiniteElement)
    loadsVector = zeros(Float64, 3 * size(pars.mesh.nodes)[1])

    for (element, load) in pars.load
        elNum = element[1]
        direction = loadDirection(element[2])
        F = elementLoad3D(elNum, pars, load, direction, intOrder, elemTypeInd)
        loadLocalNodes = nodesFromDirection(Int(direction), elemTypeInd)

        if (size(element)[1] - 2 != size(loadLocalNodes)[1])
            println("Incorrect input load")
            return nothing
        end

        loadGlobalNodes = [element[i] for i in 3:size(element)[1]]

        for i in eachindex(loadLocalNodes)
            localIndex = loadLocalNodes[i]
            globalIndex = loadGlobalNodes[i]
            loadsVector[3 * globalIndex - 2] += F[2 * localIndex - 2]
            loadsVector[3 * globalIndex - 1] += F[3 * localIndex - 1]
            loadsVector[3 * globalIndex] += F[3 * localIndex]
        end
    end

    return loadsVector
end  # assemblyLoads
