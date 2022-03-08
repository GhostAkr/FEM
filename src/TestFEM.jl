module TestFEM

using DelimitedFiles

export verify_example

function verify_example(meshPath::String, dataPath::String, result::Array)
    answer = nothing
    if meshPath == "examples/SmallPlate/Mesh.med" &&
        dataPath == "examples/SmallPlate/Task.json"
        answer = small_example()
    elseif meshPath == "examples/Beam/Mesh.med" &&
        dataPath == "examples/Beam/Task.json"
        answer = beam_example()
    elseif meshPath == "examples/Beam3D/SmallTask/Mesh.med" &&
        dataPath == "examples/Beam3D/SmallTask/TaskStretch.json"
        answer = beam3dsmall_example()
    elseif meshPath == "examples/Beam3D/BigTask/Mesh.med" &&
        dataPath == "examples/Beam3D/BigTask/TaskBind.json"
        answer = beam3dbig_example()
    else
        @warn "Example answer wasn't found"
        return false
    end

    if result ≈ answer
        return true
    else
        return false
    end
end

# Tests from examples

function small_example()
    answer_path = "examples/SmallPlate/Answer"
    answer = readdlm(answer_path, '\t', Float64, '\n')
    return answer
end

function beam_example()
    answer_path = "examples/Beam/Answer"
    answer = readdlm(answer_path, '\t', Float64, '\n')
    return answer
end

function beam3dsmall_example()
    answer_path = "examples/Beam3D/SmallTask/AnswerStretch"
    answer = readdlm(answer_path, '\t', Float64, '\n')
    return answer
end

function beam3dbig_example()
    answer_path = "examples/Beam3D/BigTask/AnswerBind"
    answer = readdlm(answer_path, '\t', Float64, '\n')
    return answer
end

end  # TestFEM
