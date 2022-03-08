# File contains commands which correspond to preloaded examples. To run example you need to
# load all necessary modules by launching src/Main.jl file in Julia REPL. That will create
# all necessary variables in current Julia session. 

# Each example below is incapsulated in code cell according to VS Code Julia extension 
# syntax ('##' separator).

# 2D models

# SmallPlate example
##
CoreFEM.fem2D("examples/SmallPlate/Mesh.med", "examples/SmallPlate/Task.json", Quad4TypeID)
##

# Beam example
##
CoreFEM.fem2D("examples/Beam/Mesh.med", "examples/Beam/Task.json", Quad4TypeID)
##

# 3D models

# Beam3D/SmallTask
##
CoreFEM.elasmech_3d("examples/Beam3D/SmallTask/Mesh.med",
    "examples/Beam3D/SmallTask/TaskStretch.json", Iso8Pts3DTypeID)
##

# Beam3D/BigTask
##
CoreFEM.elasmech_3d("examples/Beam3D/BigTask/Mesh.med", 
    "examples/Beam3D/BigTask/TaskBind.json", Iso8Pts3DTypeID)
##

# TODO: Add examples of non-local model usage.
# impact_distance = 5
# beta_loc = 0.8
# beta_nonloc = 0.2
# CoreFEM.elasmech_3d_nonloc("examples/Beam3D/SmallTask/Mesh.med", 
    # "examples/Beam3D/SmallTask/TaskStretch.json", impact_distance, beta_loc, beta_nonloc, 
    # Iso8Pts3DTypeID)
