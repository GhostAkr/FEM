# Input parameters handling: material parameters, mesh

#module Input

export processPars, importMesh, testBC, testLoad, testMaterialProperties, printProcessPars

include("Material.jl")
include("InputVars.jl")

using MeshFEM

# Test input
testMaterialProperties() = Dict(poisC => 0.3, youngMod => 200 * 10^6)
# testLoad() = Dict(5 => [100, 0], 10 => [100, 0], 15 => [100, 0], 20 => [100, 0], 25 => [100, 0])
testLoad() = Dict(3 => [100, 0], 6 => [100, 0], 9 => [100, 0])
# testBC() = Dict(1 => fixedXY, 6 => fixedXY, 11 => fixedXY, 16 => fixedXY, 21 => fixedXY)
testBC() = Dict(1 => fixedXY, 4 => fixedXY, 7 => fixedXY)

function printProcessPars(processPars::processPars)
    print("Nodes:\n")
    printNodesMesh2D(processPars.mesh)
    print("Elements:\n")
    printElementsMesh2D(processPars.mesh)
    print("Material:\n")
    print(processPars.materialProperties, "\n")
    print("Boundary conditions:\n")
    print(processPars.bc, "\n")
    print("Loads:\n")
    print(processPars.load, "\n")
end

#end  # Input

