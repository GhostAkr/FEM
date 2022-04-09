# Base variables and types of input part

include("LoadVars.jl")
include("MaterialVars.jl")
include("NonLocVars.jl")

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

    # Non-local model parameters
    nlpars::Union{NonLocPars, Nothing}
end

# Constructors

"""
    ProcessPars(matprops::Dict{materialProperty, Real}, bc::Dict{Int, bc}, 
        oad::Dict{Array{Int, 1}, Array{Float64, 1}}, mesh::Mesh2D_T)

Construct `ProcessPars` instance where all optional arguments are passed as `nothing`.

# Arguments
- `matprops::Dict{materialProperty, Real}`: material properties;
- `bc::Dict{Int, bc}`: boundary conditions;
- `load::Dict{Array{Int, 1}, Array{Float64, 1}}`: dictionary of applied loads;
- `mesh::Mesh2D_T`: model mesh;
- `nlpars::Union{NonLocPars, Nothing}` - optional parameters of non-local model.
"""
function ProcessPars(matprops::Dict{materialProperty, Real}, bc::Dict{Int, bc}, 
    load::Dict{Array{Int, 1}, Array{Float64, 1}}, mesh::Mesh2D_T, 
    nlpars::Union{NonLocPars, Nothing} = nothing
)
    return ProcessPars(matprops, bc, load, mesh, nlpars)
end

"""
    isnonloc(params::ProcessPars)

Check if given `params` can be applied to non-local simulation. Instance of `ProcessPars`
is compatible with non-local model when it's field `nlpars` is not `nothing`. Return true if
given instance can be applied to non-local simulation and false - otherwise.

# Arguments
`params::ProcessPars`: simulation parameters.
"""
function isnonloc(params::ProcessPars)
    return isnothing(pars.nlpars)
end
