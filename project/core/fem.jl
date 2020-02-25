# File with concrete FEM algorithms

module Mod_FEM
    # Include block
    include("material.jl")
    include("input.jl")
    include("../mesh/mesh_t.jl")

    # Import block
    import Mod_Input
    import Mod_Mesh_T

    function fem2D()
        parameters = Mod_Input.processPars(undef, undef, undef, undef)
        Mod_Input.importMesh(parameters.mesh, Mod_Mesh_T.generateTestMesh2D())
        Mod_Input.testMaterialProperties(parameters)
        Mod_Input.testLoad(parameters)
        Mod_Input.testBC(parameters)
        Mod_Input.printProcessPars(parameters)
    end
end  # Mod_FEM
