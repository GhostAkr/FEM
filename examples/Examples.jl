# File contains commands which correspond to preloaded examples. To run example you need to
# load all necessary modules by launching src/Main.jl file in Julia REPL. That will create
# all necessary variables in current Julia session. 

# Each example below is incapsulated in code cell according to VS Code Julia extension 
# syntax ('##' separator).

# TODO: Check all example models for it's relevance.

# FIXME: examples folder contains some unused examples. Need to check each one for it's
# relevance.

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

# Beam3D example
##
CoreFEM.elasmech_3d("examples/Beam3D/Beam3D.med", "examples/Beam3D/Beam3D.json", Iso8Pts3DTypeID)
##

# Beam3D example (non-local model)
##
impact_distance = 5
beta_loc = 0.8
beta_nonloc = 0.2
CoreFEM.elasmech_3d_nonloc("examples/Beam3D/SmallTask/Mesh.med", 
    "examples/Beam3D/SmallTask/TaskStretch.json", impact_distance, beta_loc, beta_nonloc, 
    Iso8Pts3DTypeID)
##

# Beam3D/MidTask example (non-local model)
##
impact_distance = 5
beta_loc = 1
beta_nonloc = 0
CoreFEM.elasmech_3d_nonloc("examples/Beam3D/MidTask/Mesh.med", 
    "examples/Beam3D/MidTask/TaskBind.json", impact_distance, beta_loc, beta_nonloc, 
    Iso8Pts3DTypeID)
##

# Beam3D example (non-local model, 2D analogue)
##
impact_distance = 5
beta_loc = 0.8
beta_nonloc = 0.2
CoreFEM.elasmech_3d_nonloc("examples/Beam3D/2DAnalogue/Mesh.med", 
    "examples/Beam3D/2DAnalogue/TaskStretch.json", impact_distance, beta_loc, beta_nonloc, 
    Iso8Pts3DTypeID)
##
