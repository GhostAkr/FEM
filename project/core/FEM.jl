# File with concrete FEM algorithms

module FEM
    # Include block
    #include("Material.jl")
    #include("Input.jl")
    #include("../mesh/Mesh_T.jl")

    # Import block
    using Input
    using Mesh_T

    function fem2D()
        #parameters = Input.processPars
        #parameters = Input.processPars()
        importMesh(parameters.mesh, Mesh_T.generateTestMesh2D())
        Input.testMaterialProperties(parameters)
        Input.testLoad(parameters)
        Input.testBC(parameters)
        Input.printProcessPars(parameters)
    end
end  # FEM
