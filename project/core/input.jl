# Input parameters handling: material parameters, mesh

module Mod_Input

    # Include block
    include("material.jl")
    include("load.jl")
    include("../mesh/mesh_t.jl")

    # Import block
    import Mod_Mat_Props
    import Mod_Mesh_T
    import Mod_Load

    struct processPars
        materialProperties::Dict{Mod_Mat_Props.materialProperty, Real}
        bc::Dict{Int, Mod_Load.bc}  # {NumberOfPoint, TypeOfBC}
        load::Dict{Int, Vector{Real}}  # {NumberOfPoint, ForceVector}
        mesh::Mod_Mesh_T.Mesh2D_T
    end

    function importMesh(meshTo::Mod_Mesh_T.Mesh2D_T, meshFrom::Mod_Mesh_T.Mesh2D_T)
        meshTo = meshFrom
    end

    function testMateriaProperties(targetProcessPars::processPars)  # Test material (steel)
        targetProcessPars = Dict(Mod_Mat_Props.poisC => 0.3, Mod_Mat_Props.youngMod => 2000)
    end

end  # Mod_Input

