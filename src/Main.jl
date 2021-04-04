# To provide correct modules loading it's important to execeute Set_Path.jl first.

using CoreFEM
using ElementTypes

# Example tests
CoreFEM.fem2D("examples/SmallPlate/PlateMeshSmall.med", "examples/SmallPlate/SmallPlate_Data.json", Quad4TypeID)
CoreFEM.fem2D("examples/Beam/BeamMesh.med", "examples/Beam/BeamData.json", Quad4TypeID)

CoreFEM.fem3D("unparsed", "unparsed", Iso8Pts3DTypeID)
