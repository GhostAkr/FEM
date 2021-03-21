using Unicode
using CoreFEM
import JSON

export read_params_JSON!

"""
    read_params_JSON!(file_path::String, params::CoreFEM.processPars)

Read model parameters from given JSON file.

# Arguments
- `file_path::String`: path to given file;
- `params::CoreFEM.processPars`: target process parameters.
"""
function read_params_JSON!(file_path::String, params::CoreFEM.processPars)
    model_data = Dict()
    open(file_path, "r") do file
        model_data = JSON.parse(file)
    end
    
    for property in keys(model_data)
        parse_property_JSON!(property, model_data[property], params)
    end
end  # read_params_JSON!

"""
    parse_property_JSON!(property_name::String, property_data::Dict, params::CoreFEM.processPars)

Process given property.

# Arguments
- `property_name::String`: name of given property;
- `property_data::Dict`: parameters of given property;
- `params::CoreFEM.processPars`: target process parameters.
"""
function parse_property_JSON!(property_name::String, property_data::Dict, params::CoreFEM.processPars)
    property_name = Unicode.normalize(property_name, casefold = true)

    if property_name == "material"
        params.materialProperties = parse_material_JSON(property_data)
    elseif property_name == "constraints"
        params.bc = parse_constraints_JSON(property_data)
    elseif property_name == "distributedload"
        params.load = parse_loads_JSON(property_data)
    else
        @error("Given property is not supported")
    end
end  # parse_property_JSON!

"""
    parse_material_JSON(material_data::Dict)

Process given material. Return dictionary compatible with `processPars` structure.

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

Process given constraints. Return dictionary compatible with `processPars` structure.

# Arguments
- `constraints_data::Dict`: parameters of given constraints.
"""
function parse_constraints_JSON(constraints_data::Dict)
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
        else
            @error("Incorrect constraint while importing data")
        end

        # Here several types of value are supported:
        # String -- name of group of nodes;
        # Array -- group of nodes (explicit input).
        nodes = []
        if value isa String
            # TODO: Implement this case (assign apropriate nodes array by the name of group)
            @warn("Input nodes by group is not fully supported yet")
            return nothing
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

Process given loads. Return dictionary compatible with `processPars` structure.

# Arguments
- `loads_data::Dict`: parameters of given loads.
"""
function parse_loads_JSON(loads_data::Dict)
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
                # TODO: Implement this case (assign apropriate nodes array by the name of group)
                @warn("Input nodes by group is not fully supported yet")
                return nothing
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