using Unicode
using CoreFEM
using Base.Iterators
import JSON

export read_params_JSON!

"""
    read_params_JSON!(file_path::String, params::CoreFEM.ProcessPars)

Read model parameters from given JSON file.

# Arguments
- `file_path::String`: path to given file;
- `params::CoreFEM.ProcessPars`: target process parameters.
"""
function read_params_JSON!(file_path::String, params::CoreFEM.ProcessPars)
    model_data = Dict()
    open(file_path, "r") do file
        model_data = JSON.parse(file)
    end
    
    for property in keys(model_data)
        parse_property_JSON!(property, model_data[property], params)
    end
end  # read_params_JSON!

"""
    parse_property_JSON!(property_name::String, property_data::Dict, params::CoreFEM.ProcessPars)

Process given property.

# Arguments
- `property_name::String`: name of given property;
- `property_data::Dict`: parameters of given property;
- `params::CoreFEM.ProcessPars`: target process parameters.
"""
function parse_property_JSON!(property_name::String, property_data::Dict, params::CoreFEM.ProcessPars)
    property_name = Unicode.normalize(property_name, casefold = true)

    if property_name == "material"
        params.materialProperties = parse_material_JSON(property_data)
    elseif property_name == "constraints"
        params.bc = parse_constraints_JSON(property_data, params.mesh)
    elseif property_name == "distributedload"
        params.load = parse_loads_JSON(property_data, params.mesh)
    elseif property_name == "model"
        params.model = parse_model_JSON(property_data)
    else
        @error("Given property is not supported")
    end
end  # parse_property_JSON!

"""
    parse_material_JSON(material_data::Dict)

Process given material. Return dictionary compatible with `ProcessPars` structure.

# Arguments
- `material_data::Dict`: parameters of given material.
"""
function parse_material_JSON(material_data::Dict)
    resDict = Dict{CoreFEM.materialProperty, Real}()

    for parameter in keys(material_data)
        value = material_data[parameter]
        parameter = Unicode.normalize(parameter, casefold = true)
        
        if parameter == "pratio"
            push!(resDict, CoreFEM.poisC => value)
        elseif parameter == "youngmod"
            push!(resDict, CoreFEM.youngMod => value)
        else
            @error("Given material parameter is not supported")
        end
    end

    return resDict
end  # parse_material_JSON

"""
    parse_constraints_JSON(constraints_data::Dict)

Process given constraints. Return dictionary compatible with `ProcessPars` structure.

# Arguments
- `constraints_data::Dict`: parameters of given constraints.
"""
function parse_constraints_JSON(constraints_data::Dict, mesh::Mesh2D_T)
    resDict = Dict{Int, CoreFEM.bc}()

    for parameter in keys(constraints_data)
        value = constraints_data[parameter]
        parameter = Unicode.normalize(parameter, casefold = true)

        # Converting parameter to appropriate type
        constraint_type = CoreFEM.bc
        if parameter == "fixedxy"
            constraint_type = CoreFEM.fixedXY
        elseif parameter == "fixedx"
            constraint_type = CoreFEM.fixedX
        elseif parameter == "fixedy"
            constraint_type = CoreFEM.fixedY
        elseif parameter == "fixedz"
            constraint_type = CoreFEM.fixedZ
        elseif parameter == "fixedxz"
            constraint_type = CoreFEM.fixedXZ
        elseif parameter == "fixedyz"
            constraint_type = CoreFEM.fixedYZ
        elseif parameter == "fixedxyz"
            constraint_type = CoreFEM.fixedXYZ
        else
            @error("Incorrect constraint while importing data")
        end

        # Here several types of value are supported:
        # String -- name of group of nodes;
        # Array -- group of nodes (explicit input).
        nodes = []
        if value isa String
            # Getting nodes from given mesh by group name
            if !haskey(mesh.groups, value)
                @error("Group of nodes not found (constraints)")
                return nothing
            end

            nodes = mesh.groups[value]
            # Flatting nodes array
            nodes = unique(collect(flatten(nodes)))
        elseif value isa Array
            nodes = value
        else
            @error("Incorrect constraint type while importing data")
            return nothing
        end

        for node in nodes
            push!(resDict, node => constraint_type)
        end
    end

    return resDict
end  # parse_constraints_JSON

"""
    parse_loads_JSON(loads_data::Dict)

Process given loads. Return dictionary compatible with `ProcessPars` structure.

# Arguments
- `loads_data::Dict`: parameters of given loads.
"""
function parse_loads_JSON(loads_data::Dict, mesh::Mesh2D_T)
    resDict = Dict{Array{Int, 1}, Array{Float64, 1}}()
    load = []

    for parameter in keys(loads_data)
        value = loads_data[parameter]
        parameter = Unicode.normalize(parameter, casefold = true)

        if parameter == "load"
            load = value
        elseif parameter == "elements"
            elements = []

            # Here several types of value are supported:
            # String -- name of group of nodes;
            # Array -- group of nodes (explicit input).
            if value isa String
                # Getting nodes from given mesh by group name
                if !haskey(mesh.groups, value)
                    @error("Group of nodes not found (loads)")
                    return nothing
                end

                elements = mesh.groups[value]
                # Pushing element index to the beggining of element nodes set
                for element in elements
                    for elementidx in eachindex(mesh.elements)
                        if issubset(element, mesh.elements[elementidx])
                            pushfirst!(element, elementidx)
                        end
                    end
                end
            elseif value isa Array
                elements = value
            else
                @error("Incorrect load type while importing data")
                return nothing
            end

            for nodes in elements
                push!(resDict, nodes => load)
            end
        else
            @error("Incorrect load paramter while importing data")
            return nothing
        end
    end

    return resDict
end

"""
    parse_model_JSON(model_data::Dict)

Process given model. Return parsed `ModelPars` structure.

# Arguments
`model_data::Dict`: parameters of given model.
"""
function parse_model_JSON(model_data::Dict)
    model_type::ModelType = default_modeltype
    nlpars::Union{NonLocPars, Nothing} = nothing

    for parameter in keys(model_data)
        value = model_data[parameter]
        parameter = Unicode.normalize(parameter, casefold = true)

        if parameter == "type"
            model_type = convert_model_type(value)
        elseif parameter == "nonlocparams"
            nlpars = NonLocPars(value["LocImpact"], value["NonLocImpact"], 
                value["ImpactDistance"])
        else
            @error("Given model parameter is not supported")
        end
    end

    return ModelPars(model_type, nlpars)
end

"""
    convert_model_type(type_str::String)

Convert string which denotes type of simulation model to `ModelType`. Return value of 
`ModelType` type.

# Arguments
- `type_str::String`: string denoting model type.
"""
function convert_model_type(type_str::String)
    type_str = Unicode.normalize(type_str, casefold = true)
    res_type::ModelType = default_modeltype

    if type_str == "default_modeltype"
        res_type = default_modeltype
    elseif type_str == "gen_3d"
        res_type = gen_3d
    elseif type_str == "plainstress_2d"
        res_type = plainstress_2d
    elseif type_str == "plainstrain_2d"
        res_type = plainstrain_2d
    else
        @error("Unknown type of model in convert_model_type()")
    end

    return res_type
end
