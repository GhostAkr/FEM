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
    elementLoad(elementNum::Int, pars::processPars)

Return local load vector for given element.

# Arguments
- `elementNum::Int`: number of given element;
- `pars::processPars`: parameters of current model.
"""
function elementLoad(elementNum::Int, pars::processPars)
    xCoords = [pars.mesh.nodes[pars.mesh.elements[elementNum][i]][1] for i in 1:4]
    yCoords = [pars.mesh.nodes[pars.mesh.elements[elementNum][i]][2] for i in 1:4]
    load = [10; 0]  # TODO: Test load. Need to take it from input instead.
    FIntegrate(s) = transpose(Quad4Pts.displInterpMatr(1, s)) * load * Quad4Pts.DetJs(1, s, xCoords, yCoords)
    IntegrationOrder = 4
    F = multipleIntegral.gauss1DMethodMatrix(FIntegrate, IntegrationOrder)
    return F
end  # elementLoad
