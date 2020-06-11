using Unicode

include("../core/MaterialVars.jl")
include("../core/LoadVars.jl")

export readParameters

# TODO: unify properties reading

"""
    parseMaterial(materialData::String)

Process given material. Return dictionary compatible with `processPars` structure.

# Arguments
- `materialData::String`: Parameters of given material.
"""
function parseMaterial(materialData::String)
    readPos = 1
    resDict = Dict{materialProperty, Real}()
    while findnext('=', materialData, readPos) !== nothing
        spacePos = findnext(' ', materialData, readPos)
        paramName = Unicode.normalize(String(SubString(materialData, readPos:(spacePos - 1))), casefold = true)
        if paramName == "pratio"
            eqSign = findnext('=', materialData, spacePos)
            lineEnd = findnext('\n', materialData, spacePos)
            value = parse(Float64, SubString(materialData, (eqSign + 1):(lineEnd - 2)))
            push!(resDict, poisC => value)
        elseif paramName == "youngmod"
            eqSign = findnext('=', materialData, spacePos)
            lineEnd = findnext('\n', materialData, spacePos)
            value = parse(Float64, SubString(materialData, (eqSign + 1):(lineEnd - 2)))
            push!(resDict, youngMod => value)
        else
            println("Given material parameter is not supported")
        end
        readPos = findnext('\n', materialData, readPos) + 1
    end
    return resDict
end  # parseMaterial


"""
    parseConstraints(constraintsData::String)

Process given constraints. Return dictionary compatible with `processPars` structure.

# Arguments
- `constraintsData::String`: Parameters of given constraints.
"""
function parseConstraints(constraintsData::String)
    readPos = 1
    resDict = Dict{Int, bc}()
    while findnext('=', constraintsData, readPos) !== nothing
        spacePos = findnext(' ', constraintsData, readPos)
        paramName = Unicode.normalize(String(SubString(constraintsData, readPos:(spacePos - 1))), casefold = true)
        if paramName == "fixedx"
            bktFirst = findnext('(', constraintsData, spacePos)
            bktLast = findnext(')', constraintsData, bktFirst + 1)
            valueStartPos = bktFirst + 1
            while true
                valueEndPos = findnext(',', constraintsData, valueStartPos)
                if valueEndPos === nothing
                    valueEndPos = bktLast
                elseif (valueEndPos > bktLast)
                    valueEndPos = bktLast
                end
                value = parse(Int, SubString(constraintsData, valueStartPos:(valueEndPos - 1)))
                valueStartPos = valueEndPos + 1
                push!(resDict, value => fixedX)
                nextComm = findnext(',', constraintsData, valueStartPos)
                if nextComm === nothing && valueStartPos > bktLast
                    break
                end
            end
        elseif paramName == "fixedy"
            bktFirst = findnext('(', constraintsData, spacePos)
            bktLast = findnext(')', constraintsData, bktFirst + 1)
            valueStartPos = bktFirst + 1
            while true
                valueEndPos = findnext(',', constraintsData, valueStartPos)
                if valueEndPos === nothing
                    valueEndPos = bktLast
                elseif (valueEndPos > bktLast)
                    valueEndPos = bktLast
                end
                value = parse(Int, SubString(constraintsData, valueStartPos:(valueEndPos - 1)))
                valueStartPos = valueEndPos + 1
                push!(resDict, value => fixedY)
                nextComm = findnext(',', constraintsData, valueStartPos)
                if nextComm === nothing && valueStartPos > bktLast
                    break
                end
            end
        elseif paramName == "fixedxy"
            bktFirst = findnext('(', constraintsData, spacePos)
            bktLast = findnext(')', constraintsData, bktFirst + 1)
            valueStartPos = bktFirst + 1
            while true
                valueEndPos = findnext(',', constraintsData, valueStartPos)
                if valueEndPos === nothing
                    valueEndPos = bktLast
                elseif (valueEndPos > bktLast)
                    valueEndPos = bktLast
                end
                value = parse(Int, SubString(constraintsData, valueStartPos:(valueEndPos - 1)))
                valueStartPos = valueEndPos + 1
                push!(resDict, value => fixedXY)
                nextComm = findnext(',', constraintsData, valueStartPos)
                if nextComm === nothing && valueStartPos > bktLast
                    break
                end
            end
        else
            println("Given constraints parameter is not supported")
        end
        readPos = findnext('\n', constraintsData, readPos) + 1
    end
    println(resDict)
    return resDict
end  # parseConstraints

function parseLoads(loadsData::String)
    readPos = 1
    resDict = Dict{materialProperty, Real}()
    while findnext('=', materialData, readPos) !== nothing
        
        readPos = findnext('\n', materialData, readPos) + 1
    end
    return resDict
end  # parseLoads

"""
    parseProperty(propertyName::String)

Process given property.

# Arguments
- `propertyName::String`: Name of given property;
- `propertyData::String`: Parameters of given property.
"""
function parseProperty(propertyName::String, propertyData::String)
    propertyName = Unicode.normalize(propertyName, casefold = true)
    if propertyName == "material"
        parseMaterial(propertyData)
    elseif propertyName == "constraints"
        parseConstraints(propertyData)
    else 
        println("Given property is not supported")
    end
end  # parseProperty

"""
    readParameters(filePath::String)

Read model parameters from given file.

# Arguments
- `filePath::String`: Path to given file.
"""
function readParameters(filePath::String)
    fileContents = []
    open(filePath, "r") do file
        fileContents = read(file, String)
    end
    readPos = 1
    while findnext("****", fileContents, readPos) !== nothing
        readPos = last(findnext("****", fileContents, readPos))
        readPos = findnext('\n', fileContents, readPos) + 1
        lineEnd = findnext('\n', fileContents, readPos)
        propertyName = String(SubString(fileContents, readPos:(lineEnd - 2)))
        readPos = lineEnd + 1
        nextProperty = findnext("****", fileContents, readPos)
        if nextProperty === nothing
            nextProperty = length(fileContents)[1]
            propertyEnd = length(fileContents)[1]
        else
            nextProperty = first(nextProperty)
            propertyEnd = findprev('\n', fileContents, nextProperty)
        end
        property = String(SubString(fileContents, readPos:(propertyEnd - 1)))
        parseProperty(propertyName, property)
    end
end  # readParameters



# function readParameters(filePath::String)
#     fileContents = []
#     open(filePath, "r") do file
#         fileContents = split(read(file, String))
#     end
#     println(fileContents)
#     currentParamPos = findnext(isequal("****"), fileContents, 1)
#     while currentParamPos !== nothing
#         if currentParamPos + 1 > size(fileContents)[1]
#             println("Wrong parameters file")
#         end
#         propertyName = fileContents[currentParamPos + 1]
#         lastPropertyIndex = findnext(isequal("****"), fileContents, currentParamPos + 1)
#         if lastPropertyIndex === nothing
#             lastPropertyIndex = size(fileContents)[1]
#         else
#         end
#         allProperties = [fileContents[i] for i in (currentParamPos + 2):(lastPropertyIndex - 1)]
#         println("Properties of ", propertyName)
#         println(allProperties)
#         currentParamPos = findnext(isequal("****"), fileContents, currentParamPos + 1)
#     end
# end  # readParameters
