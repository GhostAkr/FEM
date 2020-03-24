module CoreFEM

include("Material.jl")
include("Load.jl")
include("Input.jl")
include("StiffnessMatrix.jl")

using MeshFEM

export fem2D

function fem2D()
    parameters = processPars(testMaterialProperties(), testBC(), testLoad(), generateTestMesh2D())
    nu = parameters.materialProperties[poisC]
    E = parameters.materialProperties[youngMod]
    elasticityMatrix = [1 nu 0; nu 1 0; 0 0 ((1 - nu) / 2)]
    elasticityMatrix *= E / (1 - nu ^ 2)
    #for elementNum in eachindex(parameters.mesh.elements)
    #    println("Element #", elementNum)
    #    stiffnessMatrix(elasticityMatrix, parameters, elementNum)
    #end
    stiffnessMatrix(elasticityMatrix, parameters, 1)
end

end  # Core