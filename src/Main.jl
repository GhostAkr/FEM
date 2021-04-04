# To provide correct modules loading it's important to execeute Set_Path.jl first.

using CoreFEM
using ElementTypes

# Example tests
CoreFEM.fem2D("examples/QuarterPlate8N/QuadrPlate8N_MeshNew.dat", "examples/QuarterPlate8N/QuadrPlate8N_DataNew", Quad8TypeID)
CoreFEM.fem2D("examples/QuarterPlate/Quarter_Mesh.dat", "examples/QuarterPlate/Quarter_Data.json", Quad4TypeID)
CoreFEM.fem2D("examples/SmallPlate/PlateMeshSmall.med", "examples/SmallPlate/SmallPlate_Data.json", Quad4TypeID)
CoreFEM.fem2D("examples/Beam/BeamMesh.med", "examples/Beam/BeamData.json", Quad4TypeID)

CoreFEM.fem3D("unparsed", "unparsed", Iso8Pts3DTypeID)
