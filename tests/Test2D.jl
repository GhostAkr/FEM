# Application entry point
include("../src/Main.jl")

# Simple 2D example
CoreFEM.fem2D("examples/Beam/BeamMesh.med", "examples/Beam/BeamData.json", Quad4TypeID)
