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
CoreFEM.fem2D("examples/SmallPlate/PlateMeshSmall.med", 
    "examples/SmallPlate/SmallPlate_Data.json", Quad4TypeID)
##

# Beam example
##
CoreFEM.fem2D("examples/Beam/BeamMesh.med", "examples/Beam/BeamData.json", Quad4TypeID)
##

# Beam2DBindAnsys example
##
CoreFEM.fem2D("examples/Beam2DBindAnsys/Beam2DAnsys.med", 
    "examples/Beam2DBindAnsys/Beam2DBindAnsys.json", Quad4TypeID)
##

# 3D models

# Beam3D example
##
CoreFEM.elasmech_3d("examples/Beam3D/Beam3D.med", "examples/Beam3D/Beam3D.json", Iso8Pts3DTypeID)
##

# Beam3DBindSimple example
##
CoreFEM.elasmech_3d("examples/Beam3DBindSimple/Beam3DBindSimple.med", 
    "examples/Beam3DBindSimple/Beam3DBindSimple.json", Iso8Pts3DTypeID)
##

# Beam3DBindAnsys example
##
CoreFEM.elasmech_3d("examples/Beam3DBindAnsys/Beam3DBindAnsys.med", 
    "examples/Beam3DBindAnsys/Beam3DBindAnsys.json", Iso8Pts3DTypeID)
##

# NonLoc2DPlate example (local model)
##
CoreFEM.fem2D("examples/NonLoc2DPlate/Mesh.med", 
    "examples/NonLoc2DPlate/Data.json", Quad4TypeID)
##

# NonLoc2DPlate example (non-local model)
##
impact_distance = 80
beta_loc = 0.8
beta_nonloc = 0.2
CoreFEM.elasmech_2d_nonloc("examples/NonLoc2DPlate/Mesh.med", 
    "examples/NonLoc2DPlate/Data.json", impact_distance, beta_loc, beta_nonloc, Quad4TypeID)
##
