# File with concrete FEM algorithms

module Mod_FEM
    # Include block
    include("input.jl")
    include("../mesh/mesh_t.jl")

    # Import block
    import Mod_Input
    import Mod_Mesh_T

    function fem2D()
        parameters = Mod_Input.processPars(undef, undef, undef, undef)
        Mod_Input.importMesh(parameters.mesh, Mod_Mesh_T.generateTestMesh2D())
    end
end  # Mod_FEM
