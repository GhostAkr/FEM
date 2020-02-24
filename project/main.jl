# Main file. Used as application start point

# Main module
module Mod_Main

    # Import block
    import Base

    # Include block
    include("include.jl")

    testMesh = Mod_Mesh_T.generateTestMesh2D()
    print("Elements in test mesh:\n")
    Mod_Mesh_T.printElementsMesh2D(testMesh)
    print("\n")
    print("Nodes in test mesh\n")
    Mod_Mesh_T.printNodesMesh2D(testMesh)
end  # Mod_Main
