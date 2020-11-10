# To provide correct modules loading it's important to execeute Set_Path.jl first.

using CoreFEM
using multipleIntegral
using MeshFEM
using BaseInterface
using Quad8Pts

# Res = CoreFEM.fem2D("examples/QuarterPlate/Quarter_Mesh.dat", "examples/QuarterPlate/Quarter_Data")
# Res = CoreFEM.fem2D("examples/QuarterPlate/Quarter_Mesh35mm.dat", "examples/QuarterPlate/Quarter_Data35mm")
# Res = CoreFEM.fem2D("examples/QuarterPlate/Small_Mesh.dat", "examples/QuarterPlate/Small_Input")

# Res = CoreFEM.fem2D("examples/SimplePlate8N/SimplePlate8N_Mesh.dat", "examples/SimplePlate8N/SimplePlate8N_Data")
# Res = CoreFEM.fem2D("examples/SimplePlate/Squared2x2_Mesh.dat", "examples/SimplePlate/Squared2x2_Data")
# Res = CoreFEM.fem2D("examples/QuarterPlate8N/QuadrPlate8N_Mesh.dat", "examples/QuarterPlate8N/QuadrPlate8N_Data")
Res = CoreFEM.fem2D("examples/QuarterPlate8N/QuadrPlate8N_MeshNew.dat", "examples/QuarterPlate8N/QuadrPlate8N_DataNew")
