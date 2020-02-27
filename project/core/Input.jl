# Input parameters handling: material parameters, mesh

module Input

export processPars, importMesh, testBC, testLoad, testMaterialProperties, printProcessPars

#include("Material.jl")
#include("Load.jl")
#include("../mesh/Mesh_T.jl")

using Material
using Mesh_T
using Load

struct processPars
    materialProperties::Dict{Material.materialProperty, Real}
    bc::Dict{Int, Load.bc}  # {NumberOfPoint, TypeOfBC}
    load::Dict{Int, Vector{Real}}  # {NumberOfPoint, ForceVector}
    mesh::Mesh_T.Mesh2D_T
end

function importMesh(meshTo::Mesh_T.Mesh2D_T, meshFrom::Mesh_T.Mesh2D_T)
    meshTo = meshFrom
end

function testMaterialProperties(targetProcessPars::processPars)  # Test material (steel)
    targetProcessPars = Dict(Material.poisC => 0.3, Material.youngMod => 2000)
end

function testLoad(targetProcessPars::processPars)
    force = Vector(1, 0)
    targetProcessPars.load = Dict(5 => force, 10 => force, 15 => force, 20 => force, 25 => force)
end

function testBC(targetProcessPars::processPars)
    bc = Material.fixedXY
    targetProcessPars.bc = Dict(1 => bc, 6 => bc, 11 => bc, 16 => bc, 21 => bc)
end

function printProcessPars(processPars::processPars)
    print("Nodes:\n")
    Mesh_T.printNodesMesh2D(processPars.mesh)
    print("Elements:\n")
    Mesh_T.printElementsMesh2D(processPars.mesh)
    print("Material:\n")
    print(processPars.materialProperties, "\n")
    print("Boundary conditions:\n")
    print(processPars.bc, "\n")
    print("Loads:\n")
    print(processPars.load, "\n")
end

end  # Input

