include("LoadVars.jl")
include("InputVars.jl")

export applyFixedX, applyFixedY, applyFixedXY, elementLoad

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

function applyFixedXY(node::Int, loads::Array, globalK::Array)
    applyFixedX(node, loads, globalK)
    applyFixedY(node, loads, globalK)
end  # applyFixedXY

function elementLoad(elementNum::Int, pars::processPars)
    xCoords = [pars.mesh.nodes[pars.mesh.elements[elementNum][i]][1] for i in 1:4]
    yCoords = [pars.mesh.nodes[pars.mesh.elements[elementNum][i]][2] for i in 1:4]
    load = [100000; 0]  # TODO: Test load. Need to take it from input instead.
    FIntegrate(r, s) = transpose(Quad4Pts.displInterpMatr(r, s)) * load * Quad4Pts.DetJs(r, s, xCoords, yCoords)
    
    IntegrationOrder = 2
    F = multipleIntegral.gaussMethodMatrix(FIntegrate, IntegrationOrder)
    println("F: \n", F, "\n")
    return F
end  # elementLoad
