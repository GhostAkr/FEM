# Input parameters handling: material parameters, mesh

#module Input

export processPars, importMesh, testBC, testLoad, testMaterialProperties, printProcessPars

include("Material.jl")

#using Material
using MeshFEM

struct processPars
    materialProperties::Dict{materialProperty, Real}
    bc::Dict{Int, bc}  # {NumberOfPoint, TypeOfBC}
    load::Dict{Int, Vector{Real}}  # {NumberOfPoint, ForceVector}
    mesh::Mesh2D_T
end

# Test input
testMaterialProperties() = Dict(poisC => 0.3, youngMod => 2000)
testLoad() = Dict(5 => [1, 0], 10 => [1, 0], 15 => [1, 0], 20 => [1, 0], 25 => [1, 0])
testBC() = Dict(1 => fixedXY, 6 => fixedXY, 11 => fixedXY, 16 => fixedXY, 21 => fixedXY)

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

