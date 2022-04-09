# Base variables and types of input part

include("LoadVars.jl")
include("MaterialVars.jl")

using MeshFEM

export ProcessPars

"""
    ProcessPars

Parameters of current model. Contains 4 fields:
1. `materialProperties::Dict{materialProperty, Real}` - dictionary of material properties;
2. `bc::Dict{Int, bc}` - dictionary of boundary conditions;
3. `load::Dict{Int, Vector{Real}}` - dictionary of applied loads;
4. `mesh::Mesh2D_T` - model mesh.
"""
mutable struct ProcessPars
    materialProperties::Dict{materialProperty, Real}
    bc::Dict{Int, bc}  # {NumberOfPoint, TypeOfBC}
    load::Dict{Array{Int, 1}, Array{Float64, 1}}  # {NumberOfPoint, ForceVector}
    mesh::Mesh2D_T
end
