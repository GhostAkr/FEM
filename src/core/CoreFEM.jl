module CoreFEM

include("Material.jl")
include("Load.jl")
include("Input.jl")

using MeshFEM

function fem2D()
    parameters = processPars(testMaterialProperties(), testBC(), testLoad(), generateTestMesh2D())
    printProcessPars(parameters)
end

end  # Core