include("NonLocVars.jl")

"""
    ModelType

Possible types of simulation model:
- `default`: value which can be set if model is invalid;
- `gen_3d`: usual 3D model;
- `plainstress_2d`: plain stress 2D model;
- `plainstrain_2d`: plain strain 2D model.
"""
@enum ModelType begin
    default_modeltype
    gen_3d
    plainstress_2d
    plainstrain_2d
end

"""
    ModelPars
"""
struct ModelPars
    type::ModelType
    nlpars::Union{NonLocPars, Nothing}
end

"""
    ModelPars(type::ModelType, nlpars::Union{NonLocPars, Nothing} = nothing)

Construct `ModelPars` instance where all optional arguments are passed as `nothing`.

# Arguments
- `type::ModelType`: model type;
- `nlpars::Union{NonLocPars, Nothing}`: optional parameters of non-local model.
"""
function ModelPars(type::ModelType, nlpars::Union{NonLocPars, Nothing} = nothing)
    return ModelPars(type, nlpars)
end
