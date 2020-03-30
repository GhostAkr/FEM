module CoreFEM

include("Material.jl")
include("Load.jl")
include("Input.jl")
include("StiffnessMatrix.jl")

using MeshFEM
using DelimitedFiles

export fem2D

function fem2D()
    parameters = processPars(testMaterialProperties(), testBC(), testLoad(), generateTestMesh2D())
    nu = parameters.materialProperties[poisC]
    E = parameters.materialProperties[youngMod]
    elasticityMatrix = [1 nu 0; nu 1 0; 0 0 ((1 - nu) / 2)]
    elasticityMatrix *= E / (1 - nu ^ 2)
    names = Array{String}(undef, 16)
    for elementNum in eachindex(parameters.mesh.elements)
        println("Element #", elementNum)
        K = stiffnessMatrix(elasticityMatrix, parameters, elementNum)
        fileName = "K" * string(elementNum)
        file = open(fileName, "w")
        close(file)
        open(fileName, "a") do file
            writedlm(file, K)
        end
    end
end

end  # Core