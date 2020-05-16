include("LoadVars.jl")

using MeshFEM

export processPars

struct processPars
    materialProperties::Dict{materialProperty, Real}
    bc::Dict{Int, bc}  # {NumberOfPoint, TypeOfBC}
    load::Dict{Int, Vector{Real}}  # {NumberOfPoint, ForceVector}
    mesh::Mesh2D_T
end
