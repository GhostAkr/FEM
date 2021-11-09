module TestFEM

using DelimitedFiles

export verify_example

function verify_example(meshPath::String, dataPath::String, result::Array)
    answer = nothing
    if meshPath == "examples/QuarterPlate/Quarter_Mesh.dat" &&
        dataPath == "examples/QuarterPlate/Quarter_Data.json" 
        answer = quarter_example()
    elseif meshPath == "examples/QuarterPlate8N/QuadrPlate8N_MeshNew.dat" &&
        dataPath == "examples/QuarterPlate8N/QuadrPlate8N_DataNew"
        answer = quarter8N_example()
    elseif meshPath == "examples/SmallPlate/PlateMeshSmall.med" &&
        dataPath == "examples/SmallPlate/SmallPlate_Data.json"
        answer = small_example()
    elseif meshPath == "examples/Beam/BeamMesh.med" &&
        dataPath == "examples/Beam/BeamData.json"
        answer = beam_example()
    elseif meshPath == "examples/Beam3DBindAnsys/Beam3DBindAnsys.med" &&
        dataPath == "examples/Beam3DBindAnsys/Beam3DBindAnsys.json"
        answer = beam3DBindAnsys_example()
    elseif meshPath == "examples/Beam3D/Beam3D.med" &&
        dataPath == "examples/Beam3D/Beam3D.json"
        answer = beam3d_example()
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
    answer_path = "examples/SmallPlate/SmallPlate_Answer"
    answer = readdlm(answer_path, '\t', Float64, '\n')
    return answer
end

function quarter_example()
    answer_path = "examples/QuarterPlate/Quarter_Answer"
    answer = readdlm(answer_path, '\t', Float64, '\n')
    return answer
end

function quarter8N_example()
    answer_path = "examples/QuarterPlate8N/QuarterPlate8N_Answer"
    answer = readdlm(answer_path, '\t', Float64, '\n')
    return answer
end

function beam_example()
    answer_path = "examples/Beam/BeamAnswer"
    answer = readdlm(answer_path, '\t', Float64, '\n')
    return answer
end

function beam3d_example()
    answer_path = "examples/Beam3D/Beam3D_Answer"
    answer = readdlm(answer_path, '\t', Float64, '\n')
    return answer
end

function beam3DBindAnsys_example()
    answer_path = "examples/Beam3DBindAnsys/Beam3DBindAnsys_Answer"
    answer = readdlm(answer_path, '\t', Float64, '\n')
    return answer
end

end  # TestFEM