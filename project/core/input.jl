# Input parameters handling: material parameters, mesh

module Mod_Input

export processPars, importMesh, testBC, testLoad, testMaterialProperties, printProcessPars

include("material.jl")
include("load.jl")
include("../mesh/mesh_t.jl")

using Mod_Mat_Props
using Mod_Mesh_T
using Mod_Load

struct processPars
    materialProperties::Dict{Mod_Mat_Props.materialProperty, Real}
    bc::Dict{Int, Mod_Load.bc}  # {NumberOfPoint, TypeOfBC}
    load::Dict{Int, Vector{Real}}  # {NumberOfPoint, ForceVector}
    mesh::Mod_Mesh_T.Mesh2D_T
end

function importMesh(meshTo::Mod_Mesh_T.Mesh2D_T, meshFrom::Mod_Mesh_T.Mesh2D_T)
    meshTo = meshFrom
end

function testMaterialProperties(targetProcessPars::processPars)  # Test material (steel)
    targetProcessPars = Dict(Mod_Mat_Props.poisC => 0.3, Mod_Mat_Props.youngMod => 2000)
end

function testLoad(targetProcessPars::processPars)
    force = Vector(1, 0)
    targetProcessPars.load = Dict(5 => force, 10 => force, 15 => force, 20 => force, 25 => force)
end

function testBC(targetProcessPars::processPars)
    bc = Mod_Mat_Props.fixedXY
    targetProcessPars.bc = Dict(1 => bc, 6 => bc, 11 => bc, 16 => bc, 21 => bc)
end

function printProcessPars(processPars::processPars)
    print("Nodes:\n")
    Mod_Mesh_T.printNodesMesh2D(processPars.mesh)
    print("Elements:\n")
    Mod_Mesh_T.printElementsMesh2D(processPars.mesh)
    print("Material:\n")
    print(processPars.materialProperties, "\n")
    print("Boundary conditions:\n")
    print(processPars.bc, "\n")
    print("Loads:\n")
    print(processPars.load, "\n")
end

end  # Mod_Input

