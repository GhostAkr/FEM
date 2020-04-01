export applyFixedX

@enum bc begin
    fixedX  # Dislpacement fixed by X
    fixedY  # Dislpacement fixed by Y
    fixedXY  # Dislpacement fixed by X and Y
end

function applyFixedX(node::Int, loads::Array, globalK::Array)
    loads[2 * node - 1] = 0
    for col in size(globalK)[2]
        globalK[2 * node - 1, col] = 0
    end
    for row in size(globalK)[1]
        globalK[row, 2 * node - 1] = 0
    end
    globalK[2 * node - 1, 2 * node - 1] = 1
end  # applyFixedX

function applyFixedY(node::Int, loads::Array, globalK::Array)
    loads[2 * node] = 0
    for col in size(globalK)[2]
        globalK[2 * node, col] = 0
    end
    for row in size(globalK)[1]
        globalK[row, 2 * node] = 0
    end
    globalK[2 * node, 2 * node] = 1
end  # applyFixedY

function applyFixedXY(node::Int, loads::Array, globalK::Array)
    applyFixedX(node, loads, globalK)
    applyFixedY(node, loads, globalK)
end  # applyFixedXY
