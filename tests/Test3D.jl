# Application entry point
include("../src/Main.jl")

# Simple 3D example
CoreFEM.fem3D("examples/Beam3DBindAnsys/Beam3DBindAnsys.med", 
    "examples/Beam3DBindAnsys/Beam3DBindAnsys.json", Iso8Pts3DTypeID)