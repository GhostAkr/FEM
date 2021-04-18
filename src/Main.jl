# To provide correct modules loading it's important to execeute Set_Path.jl first.

using CoreFEM
using ElementTypes

# Example tests
CoreFEM.fem2D("examples/SmallPlate/PlateMeshSmall.med", "examples/SmallPlate/SmallPlate_Data.json", Quad4TypeID)
CoreFEM.fem2D("examples/Beam/BeamMesh.med", "examples/Beam/BeamData.json", Quad4TypeID)

# CoreFEM.fem3D("unparsed", "unparsed", Iso8Pts3DTypeID)
CoreFEM.fem3D("examples/Beam3D/Beam3DMesh.med", "examples/Beam3D/Beam3DData.json", Iso8Pts3DTypeID)
CoreFEM.fem3D("examples/SimpleBeam3D/SimpleBeam3D.med", "examples/SimpleBeam3D/SimpleBeam3D.json", Iso8Pts3DTypeID)
