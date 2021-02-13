module TestFEM

using DelimitedFiles

export verify_example

function verify_example(meshPath::String, dataPath::String, result::Array)
    answer = nothing
    if meshPath == "examples/QuarterPlate/Quarter_Mesh.dat" &&
        dataPath == "examples/QuarterPlate/Quarter_Data" 
        answer = quarter_example()
    elseif meshPath == "examples/QuarterPlate8N/QuadrPlate8N_MeshNew.dat" &&
        dataPath == "examples/QuarterPlate8N/QuadrPlate8N_DataNew"
        answer = quarter8N_example()
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

function quarter_example()
    quarter_example_path = "examples/QuarterPlate/Quarter_Answer"
    answer = readdlm(quarter_example_path, '\t', Float64, '\n')
    return answer
end

function quarter8N_example()
    quarter_example_path = "examples/QuarterPlate8N/QuarterPlate8N_Answer"
    answer = readdlm(quarter_example_path, '\t', Float64, '\n')
    return answer
end

end  # TestFEM