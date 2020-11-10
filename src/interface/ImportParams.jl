using Unicode
using CoreFEM

# include("../core/MaterialVars.jl")
# include("../core/LoadVars.jl")

export readParameters!

# TODO: unify properties reading

"""
    parseMaterial(materialData::String)

Process given material. Return dictionary compatible with `processPars` structure.

# Arguments
- `materialData::String`: Parameters of given material.
"""
function parseMaterial(materialData::String)
    readPos = 1
    resDict = Dict{CoreFEM.materialProperty, Real}()
    while findnext('=', materialData, readPos) !== nothing
        eqSign = findnext('=', materialData, readPos)
        paramNameEndPos = eqSign - 1
        while materialData[paramNameEndPos] == ' ' || materialData[paramNameEndPos] == '\t'
            paramNameEndPos -= 1
        end
        paramName = Unicode.normalize(String(SubString(materialData, readPos:(paramNameEndPos))), casefold = true)
        if paramName == "pratio"
            # eqSign = findnext('=', materialData, paramNameEndPos)
            lineEnd = findnext('\n', materialData, paramNameEndPos)
            value = parse(Float64, SubString(materialData, (eqSign + 1):(lineEnd - 1)))
            push!(resDict, CoreFEM.poisC => value)
        elseif paramName == "youngmod"
            # eqSign = findnext('=', materialData, paramNameEndPos)
            lineEnd = findnext('\n', materialData, paramNameEndPos)
            value = parse(Float64, SubString(materialData, (eqSign + 1):(lineEnd - 1)))
            push!(resDict, CoreFEM.youngMod => value)
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
    resDict = Dict{Int, CoreFEM.bc}()
    while findnext('=', constraintsData, readPos) !== nothing
        eqSign = findnext('=', constraintsData, readPos)
        paramNameEndPos = eqSign - 1
        while constraintsData[paramNameEndPos] == ' ' || constraintsData[paramNameEndPos] == '\t'
            paramNameEndPos -= 1
        end
        paramName = Unicode.normalize(String(SubString(constraintsData, readPos:(paramNameEndPos))), casefold = true)
        if paramName == "fixedx"
            bktFirst = findnext('(', constraintsData, paramNameEndPos)
            bktLast = findnext(')', constraintsData, bktFirst + 1)
            valueStartPos = bktFirst + 1
            valueEndPos = 0
            while true
                valueEndPos = findnext(',', constraintsData, valueStartPos)
                if valueEndPos === nothing
                    valueEndPos = bktLast
                elseif (valueEndPos > bktLast)
                    valueEndPos = bktLast
                end
                value = parse(Int, SubString(constraintsData, valueStartPos:(valueEndPos - 1)))
                valueStartPos = valueEndPos + 1
                push!(resDict, value => CoreFEM.fixedX)
                nextComm = findnext(',', constraintsData, valueStartPos)
                if nextComm === nothing || valueStartPos > bktLast
                    break
                end
            end
        elseif paramName == "fixedy"
            bktFirst = findnext('(', constraintsData, paramNameEndPos)
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
                push!(resDict, value => CoreFEM.fixedY)
                nextComm = findnext(',', constraintsData, valueStartPos)
                if nextComm === nothing || valueStartPos > bktLast
                    break
                end
            end
        elseif paramName == "fixedxy"
            bktFirst = findnext('(', constraintsData, paramNameEndPos)
            bktLast = findnext(')', constraintsData, bktFirst + 1)
            valueStartPos = bktFirst + 1
            valueEndPos = 0
            while true
                valueEndPos = findnext(',', constraintsData, valueStartPos)
                if valueEndPos === nothing
                    valueEndPos = bktLast
                elseif (valueEndPos > bktLast)
                    valueEndPos = bktLast
                end
                value = parse(Int, SubString(constraintsData, valueStartPos:(valueEndPos - 1)))
                valueStartPos = valueEndPos + 1
                push!(resDict, value => CoreFEM.fixedXY)
                nextComm = findnext(',', constraintsData, valueStartPos)
                if valueEndPos == bktLast
                    break
                end
            end
        else
            println("Given constraints parameter is not supported")
        end
        readPos = findnext('\n', constraintsData, readPos) + 1
    end
    return resDict
end  # parseConstraints

function parseLoads(loadsData::String)
    readPos = 1
    resDict = Dict{Array{Int, 1}, Array{Float64, 1}}()
    elements = []
    load = []
    while findnext('=', loadsData, readPos) !== nothing
        eqSign = findnext('=', loadsData, readPos)
        paramNameEndPos = eqSign - 1
        while loadsData[paramNameEndPos] == ' ' || loadsData[paramNameEndPos] == '\t'
            paramNameEndPos -= 1
        end
        paramName = Unicode.normalize(String(SubString(loadsData, readPos:(paramNameEndPos))), casefold = true)
        if paramName == "load"
            bktFirst = findnext('(', loadsData, paramNameEndPos)
            bktLast = findnext(')', loadsData, bktFirst + 1)
            nextComm = findnext(',', loadsData, bktFirst)
            xForce = parse(Float64, SubString(loadsData, (bktFirst + 1):(nextComm - 1)))
            yForce = parse(Float64, SubString(loadsData, (nextComm + 1):(bktLast - 1)))
            load = [xForce, yForce]
        elseif paramName == "elements"
            bktFirst = findnext('(', loadsData, paramNameEndPos)
            bktLast = findnext(')', loadsData, bktFirst + 1)
            currElem = bktFirst
            while true
                elemBegin = findnext('{', loadsData, currElem)
                if elemBegin === nothing
                    break
                end
                elemEnd = findnext('}', loadsData, elemBegin)
                valueStartPos = elemBegin + 1
                valueEndPos = 0
                elemData = []
                while true
                    valueEndPos = findnext(',', loadsData, valueStartPos)
                    if valueEndPos === nothing
                        valueEndPos = elemEnd
                    elseif (valueEndPos > elemEnd)
                        valueEndPos = elemEnd
                    end
                    value = parse(Int, SubString(loadsData, valueStartPos:(valueEndPos - 1)))
                    valueStartPos = valueEndPos + 1
                    push!(elemData, value)
                    nextComm = findnext(',', loadsData, valueStartPos)
                    if valueEndPos == elemEnd
                        break
                    end
                end
                currElem = elemEnd
                push!(elements, elemData)
            end
        else
            println("Given material parameter is not supported")
        end
        readPos = findnext('\n', loadsData, readPos) + 1
    end
    for i in elements
        push!(resDict, i => load)
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
function parseProperty!(propertyName::String, propertyData::String, params::CoreFEM.processPars)
    propertyName = Unicode.normalize(propertyName, casefold = true)
    println(propertyName)
    if propertyName == "material"
        params.materialProperties = parseMaterial(propertyData)
    elseif propertyName == "constraints"
        params.bc = parseConstraints(propertyData)
    elseif propertyName == "distributedload"
        params.load = parseLoads(propertyData)
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
function readParameters!(filePath::String, params::CoreFEM.processPars)
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
            propertyEnd = findprev('\n', fileContents, propertyEnd)
        else
            nextProperty = first(nextProperty)
            propertyEnd = findprev('\n', fileContents, nextProperty)
        end
        property = String(SubString(fileContents, readPos:(propertyEnd)))
        parseProperty!(propertyName, property, params)
    end
end  # readParameters
