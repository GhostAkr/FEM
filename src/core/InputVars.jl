# Base variables and types of input part

include("LoadVars.jl")
include("MaterialVars.jl")
include("NonLocVars.jl")
include("ModelVars.jl")

using MeshFEM

export ProcessPars

"""
    ProcessPars

Parameters of current model. Contains 4 fields:
1. `materialProperties::Dict{materialProperty, Real}` - dictionary of material properties;
2. `bc::Dict{Int, bc}` - dictionary of boundary conditions;
3. `load::Dict{Int, Vector{Real}}` - dictionary of applied loads;
4. `mesh::Mesh2D_T` - model mesh;
5. `nlpars::Union{NonLocPars, Nothing}` - optional parameters of non-local model.
"""
mutable struct ProcessPars
    # Default parameters
    materialProperties::Dict{materialProperty, Real}
    bc::Dict{Int, bc}  # {NumberOfPoint, TypeOfBC}
    load::Dict{Array{Int, 1}, Array{Float64, 1}}  # {NumberOfPoint, ForceVector}
    mesh::Mesh2D_T
    model::ModelPars
end

"""
    isnonloc(params::ProcessPars)

Check if given `params` can be applied to non-local simulation. Instance of `ProcessPars`
is compatible with non-local model when field `nlpars` in `model` is not `nothing`. Return 
true if given instance can be applied to non-local simulation and false - otherwise.

# Arguments
`params::ProcessPars`: simulation parameters.
"""
function isnonloc(params::ProcessPars)
    return !isnothing(params.model.nlpars)
end
