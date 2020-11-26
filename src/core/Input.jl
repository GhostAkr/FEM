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
testLoad() = Dict([2, 4, 3, 6] => [10., 0.], [4, 4, 6, 9] => [10., 0.])

testLoad3D() = Dict([2, 4, 6, 12, 3, 9] => [10., 0., 0.])

"""
    testBC()

Return test boundary condition: fix plate on the left edge.
"""
testBC() = Dict(1 => fixedXY, 4 => fixedXY, 7 => fixedXY,)

testBC3D() = Dict(1 => fixedXYZ, 4 =>fixedXYZ, 10 => fixedXYZ, 7 => fixedXYZ, 2 => fixedZ, 5 => fixedZ, 11 => fixedZ, 8 => fixedZ, 3 => fixedZ, 6 => fixedZ, 12 => fixedZ, 9 => fixedZ)

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
