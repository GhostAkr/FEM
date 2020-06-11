# Input parameters handling: material parameters, mesh

export importMesh, testBC, testLoad, testMaterialProperties, printProcessPars

include("MaterialVars.jl")
include("InputVars.jl")

using MeshFEM

"""
    testMaterialProperties()

Return test material properties: Poisson's ratio ``\\nu = 0.3``, Young's modulus ``E = 200000``.
"""
testMaterialProperties() = Dict(poisC => 0.3, youngMod => 200000)

"""
    testLoad()

Return test load: surface force ``f^s = 10 \\text{МПа}`` on the right edge of test plate.
"""
testLoad() = Dict([2, 4] => [10., 0.], [4, 4] => [10., 0.])

"""
    testBC()

Return test boundary condition: fix plate on the left edge.
"""
testBC() = Dict(1 => fixedXY, 4 => fixedXY, 7 => fixedXY,)

"""
    printProcessPars(processPars::processPars)

Print given model parameters.

# Arguments
- `processPars::processPars`: parameters of current model.
"""
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
